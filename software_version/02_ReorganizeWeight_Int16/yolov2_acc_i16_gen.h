
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>
#include <sys/stat.h>
#include <unistd.h>

unsigned long get_file_size(const char *filename)  
{  
    struct stat buf;  
    if(stat(filename, &buf)<0)  
    {  
        return 0;  
    }  
    return (unsigned long)buf.st_size;  
}

#include "acc.h"
#include "reorder_quanti_int16.h"

#define MEM_LEN (416*416*32+208*208*32)
void generate_iofm_offset(float* in_ptr[32], float* out_ptr[32], float *Memory_buf, network *net)
{
#define ROUTE16_LEN (26*26*512)
#define CONV27_LEN (13*13*256)
#define CONV24_LEN (13*13*1024)

	float *Memory_top = Memory_buf+512;
	float *Memory_bottom = Memory_top + MEM_LEN;
	int x;
	for(x=0;x<18;x++)
	{
		if(x%2==0)
		{
			in_ptr[x] = Memory_top;
			out_ptr[x] = Memory_bottom - net->layers[x].outputs;
		}
		else
		{
			in_ptr[x] = out_ptr[x-1];
			out_ptr[x] = Memory_top;
		}
	}

	for(x=18;x<25;x++)
	{
		if(x%2==0)
		{
			in_ptr[x] = Memory_top;
			out_ptr[x] = Memory_bottom - ROUTE16_LEN - net->layers[x].outputs;
		}else
		{
			in_ptr[x] = out_ptr[x-1];
			out_ptr[x] = Memory_top;
		}
	}

	in_ptr[26] = Memory_bottom - ROUTE16_LEN;
	out_ptr[26] = Memory_top;

	in_ptr[27] = Memory_top;
	out_ptr[27] = Memory_bottom - ROUTE16_LEN - CONV24_LEN - CONV27_LEN;

	in_ptr[29] = out_ptr[27];
	out_ptr[29] = Memory_top;

	in_ptr[30] = Memory_top;
	out_ptr[30] = Memory_bottom - net->layers[30].outputs;

	in_ptr[31] = out_ptr[30];
}

void reorg_cpu(float *x, int w, int h, int c, int stride, float *out)
{
    int i,j,k;
    int out_c = c/(stride*stride);

	for(k = 0; k < c; ++k){
	    for(j = 0; j < h; ++j){
		for(i = 0; i < w; ++i){
		    int in_index  = i + w*(j + h*k);
		    int c2 = k % out_c;
		    int offset = k / out_c;
		    int w2 = i*stride + offset % stride;
		    int h2 = j*stride + offset / stride;
		    int out_index = w2 + w*stride*(h2 + h*stride*c2);
		    out[in_index] = x[out_index];
		}
	    }
	}
}

void yolov2_hls_ps(network *net, float *input)
{
	int weight_offset[32] = {864, 18432, 73728, 8192, 73728,
		294912, 32768, 294912, 1179648, 131072, 1179648, 131072,
		1179648, 4718592, 524288, 4718592, 524288, 4718592, 9437184,
		9437184, 32768, 11796480, 435200, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	int beta_offset[32] = {32, 64, 128, 64, 128, 256, 128, 256, 512, 256, 512, 256, 512, 1024,
		512, 1024, 512, 1024, 1024, 1024, 64, 1024, 425, 0, 0, 0, 0, 0, 0, 0, 0, 0};

//////////////////////////////////////////////////////////////
	short ap16_min = 0x8000;
	short ap16_max = 0x7fff;
	float ap16_range[16*2];
	for(int i=0;i<16;i++)
	{
		//printf("Q%2d:",i);
		ap16_range[2*i]   = (float)ap16_min*pow((float)2,-i);//min
		ap16_range[2*i+1] = (float)ap16_max*pow((float)2,-i);//max
		// printf("min=%.7lf,max=%.7lf\n",ap16_range[2*i],ap16_range[2*i+1]);
	}
	int maxQ_array[64];
	int add_offset[64];
	int io_maxQarray[64];
	memset(io_maxQarray, 0, sizeof(io_maxQarray));
	FILE *fp_ioq;

/*	if(access("iofm_Q.bin", F_OK)==0)//iofm_Q existed.
	{
		int file_sz = get_file_size("iofm_Q.bin");
		fp_ioq = fopen("iofm_Q.bin", "rb");
	    	if(!fp_ioq) file_error("iofm_Q.bin");
		fread(io_maxQarray, sizeof(int), file_sz/4, fp_ioq);
		fclose(fp_ioq);//iofm quanti
		printf("iofm_Q.bin existed. sz=%d\n", file_sz/4);
	}*/

	FILE *fout = fopen("weight_reorg_ap16.bin","wb");
    	if(!fout) printf("fopen weight_reorg_ap16.bin error\n");
/////////////////////////////////////////////////////////////////////////

	float *Weight_buf = (float *)calloc(203767168/4,sizeof(float));
	float *Beta_buf   = (float *)calloc(43044/4,sizeof(float));

	FILE *fp_w = fopen("weights.bin", "rb");
    	if(!fp_w) file_error("weights.bin");

	FILE *fp_b = fopen("bias.bin", "rb");
    	if(!fp_b) file_error("bias.bin");

	fread(Weight_buf, sizeof(float), 203767168/4, fp_w);
	fread(Beta_buf, sizeof(float), 43044/4, fp_b);
	
	fclose(fp_w);
	fclose(fp_b);

//leave some memories for overflow, because the load_module will load extra pixels near boundary for padding
	float *Memory_buf = (float*)calloc(MEM_LEN+512*2,sizeof(float));
	float* in_ptr[32];
	float* out_ptr[32];
	generate_iofm_offset( in_ptr, out_ptr, Memory_buf, net);

	//memcpy(in_ptr[0], input, 416*416*3*sizeof(float));//416x416x3 input_pic
////////////////////////////////////////////////////////////////////
	float* tmp_ptr_ft32 = in_ptr[0];
	for(int i=0;i<416*416*3;i++)//1st Layer input Q14
	{
		tmp_ptr_ft32[i] = ((short)(input[i]*pow(2.0,14)))*pow(2.0,-14);
	}
	float *inout_fixed_buf = (float *)calloc(sizeof(float),416*416*32); 
	io_maxQarray[0] = 14;//1st layer input Q14
//////////////////////////////////////////////////////////////

	int offset_index = 0;
	int woffset = 0;
	int boffset = 0;
	int mLoops;
	double sum_gop = 0.0;
	uint16_t ofm_w;
	uint16_t ofm_h;

    	for(int i = 0; i < net->n; ++i)
	{
        	layer l = net->layers[i];
		printf("Layer[%2d]: ",i);

		if((l.type == CONVOLUTIONAL) || (l.type == MAXPOOL))
		{
			uint8_t ltype = 255;
			if(l.type == MAXPOOL)
				ltype = LT_MAXPOOL;
			else if(l.type == CONVOLUTIONAL)
				ltype = LT_CONV;
			uint8_t kernel_size = l.size;
			uint8_t kernel_stride = l.stride;
			uint8_t pad_int = l.pad;
			if(l.type == MAXPOOL)
				pad_int = 0;		

			uint16_t ifm_num = l.c;
			uint16_t ifm_w = l.w;
			uint16_t ifm_h = l.h;
			uint16_t ofm_num = l.n;
			if(l.type == MAXPOOL)
				ofm_num = l.c;
			uint16_t ofm_w = l.out_w;
			uint16_t ofm_h = l.out_h;

			uint16_t TR,TC,TM,TN;
			TR = MIN_diy(((OnChipIB_Height-kernel_size)/kernel_stride+1),Tr);//keep Kernel_stride>=1
			TR = MIN_diy(ofm_h,TR);
			TC = MIN_diy(((OnChipIB_Width-kernel_size)/kernel_stride+1),Tc);
			TC = MIN_diy(ofm_w,TC);
			TM = MIN_diy(ofm_num,Tm);
			TN = MIN_diy(ifm_num,Tn);

			if(ltype != LT_CONV)
			{
				TM = MIN_diy(TM,TN);
				TN = 0;
			}

			bool IsNotConv = (ltype != LT_CONV);
			int mLoops = (int)ceil(((float)ofm_num)/TM);
			int OFM_num_bound = IsNotConv ? (mLoops + 2)*TM : (mLoops + 1)*TM;
			int mLoopsxTM = mLoops*TM;
			int mLoops_a1xTM = (mLoops+1)*TM;

			uint16_t TRow = (TR-1)*kernel_stride + kernel_size;
			uint16_t TCol = (TC-1)*kernel_stride + kernel_size;
			uint8_t  KK = kernel_size*kernel_size;
			uint32_t INumxKK = ifm_num*KK;//24bit

			uint32_t TRTC = (TR << 16) | TC;
			uint32_t TMTN = (TM << 16) | TN;
			uint32_t ofm_w_h = (ofm_w << 16) | ofm_h;
			uint32_t ifm_w_h = (ifm_w << 16) | ifm_h;
			uint32_t iofm_num = (ifm_num << 16) | ofm_num;
			uint32_t k_s_pad_ltype = (kernel_size << 24) | (kernel_stride << 16) | (pad_int << 8) | (ltype);
			uint32_t TRowTCol = (TRow << 16) | TCol;
			uint32_t IHW = ifm_h*ifm_w;
			uint32_t OHW = ofm_h*ofm_w;
			uint32_t KK_INumxKK = (KK << 24) | INumxKK;

			uint32_t en_bits = 0x0;//enable_bits[2:0]={IsReLU, LoadBias, IsNotConv}
			if(IsNotConv)
				en_bits |= 0x1;
			if(ltype == LT_CONV)
				en_bits |= (0x1 << 1);
			switch(l.type)
			{
				case CONVOLUTIONAL:
					printf("outputMemory:%8d;BN=%d;Activation=%d;conv  %5d %2d x%2d /%2d  %4d x%4d x%4d   ->  %4d x%4d x%4d  %5.3f BFLOPs\n",l.outputs,l.batch_normalize,l.activation,
					 l.n, l.size, l.size, l.stride, l.w, l.h, l.c, l.out_w, l.out_h, l.out_c, (2.0 * l.n * l.size*l.size*l.c/l.groups * l.out_h*l.out_w)/1000000000.);
					sum_gop += (2.0 * l.n * l.size*l.size*l.c/l.groups * l.out_h*l.out_w)/1000000000.;
					//output_w = (l.w - l.size + 2*l.pad)/l.stride + 1;
					//output_h = (l.h - l.size + 2*l.pad)/l.stride + 1;
					if(l.activation==LEAKY)
						en_bits |= (0x1 << 2);

					FPGA_Acc(in_ptr[i], out_ptr[i], Weight_buf+woffset, Beta_buf+boffset,
						k_s_pad_ltype, iofm_num, ifm_w_h, ofm_w_h, TRTC, TMTN, OFM_num_bound, mLoopsxTM, mLoops_a1xTM, 0.0f,
						TRowTCol, IHW, OHW, KK_INumxKK, en_bits);
//////////////////////////////////////////////////////#ifdef REORDER_GEN
					//Reorder_weight(Weight_buf + woffset, Weight_reorg_buf + woffset, l.size, l.c, l.n, LT_CONV, TM, TN);
					quantize_ifm_int16(out_ptr[i],inout_fixed_buf, l.outputs, ap16_range, &io_maxQarray[offset_index+1]);
					memcpy(out_ptr[i],inout_fixed_buf,l.outputs*sizeof(float));
					reorg_quantize_weight_int16(Weight_buf + woffset, l.c,l.n,l.size,TM,TN, ap16_range, &maxQ_array[offset_index], &add_offset[offset_index], fout);
/////////////////////////////////////////////////////#endif
					woffset += weight_offset[offset_index];
					boffset += beta_offset[offset_index];
					offset_index++;

					break;
				case MAXPOOL:
					printf("outputMemory:%8d;max          %d x %d / %d  %4d x%4d x%4d   ->  %4d x%4d x%4d\n",l.outputs, l.size, l.size, l.stride, l.w, l.h, l.c,
						 l.out_w, l.out_h, l.out_c);
					sum_gop += (4.0 * l.c * l.out_h*l.out_w)/1000000000.;
					//output_w = (l.w - l.size)/l.stride + 1 ;
					//output_h = (l.h - l.size)/l.stride + 1 ;
					//very strange that yolo's maxpool padding is all 0?!!
					FPGA_Acc(in_ptr[i], out_ptr[i], NULL, NULL,
						k_s_pad_ltype, iofm_num, ifm_w_h, ofm_w_h, TRTC, TMTN, OFM_num_bound, mLoopsxTM, mLoops_a1xTM, MIN_NEG,
						TRowTCol, IHW, OHW, KK_INumxKK, en_bits);

					break;
			}
		}else
		{
			switch(l.type)
			{
				case REORG:
					printf("outputMemory:%8d;reorg              /%2d  %4d x%4d x%4d   ->  %4d x%4d x%4d\n",l.outputs,  l.stride, l.w, l.h, l.c,
						 l.out_w, l.out_h, l.out_c);			
					ofm_w = 26;
					ofm_h = 32*13;

					reorg_cpu(in_ptr[i], ofm_w, ofm_h, 4, 2, out_ptr[i]);
					break;
				case ROUTE:
					printf("outputMemory:%8d;route ",l.outputs);
					for(int j = 0; j < l.n; ++j){
						printf(" %d", l.input_layers[j]);
					}
					printf("\n");
					break;
				case REGION:
					printf("outputMemory:%8d;Detection\n",l.outputs);
					forward_region_layer(l, in_ptr[i]);
					break;
			}
		}
    	}
	printf("SUM_GOP=%.7lf\n",sum_gop);

////////////////////////////////////////////////#ifdef REORDER_GEN
	fclose(fout);//weight reorg quanti

	int conv_layer_num = offset_index;

	fp_ioq = fopen("iofm_Q.bin", "wb");
    	if(!fp_ioq) file_error("iofm_Q.bin");
	fwrite(io_maxQarray, sizeof(int), conv_layer_num + 1, fp_ioq);
	fclose(fp_ioq);//iofm quanti
	
	fout = fopen("weights_ap16_offset_add.bin","wb");
    	if(!fout) printf("fopen weights_ap16_offset_add.bin error\n");
	fwrite(add_offset, sizeof(int), conv_layer_num, fout);
	fclose(fout);

	fout = fopen("weights_ap16_Q.bin","wb");
    	if(!fout) printf("fopen weights_ap16_Q.bin error\n");
	fwrite(maxQ_array, sizeof(int), conv_layer_num, fout);
	fclose(fout);	

	short *tmp_buf_reorg = (short *)calloc(11796480+65536,sizeof(short));
	if(!tmp_buf_reorg) printf("tmp_buf_reorg alloc failed.\n");

	fout = fopen("bias_ap16.bin","wb");
    	if(!fout) printf("fopen bias_ap16.bin error\n");

	quantize_bias_int16(Beta_buf, tmp_buf_reorg, beta_offset, conv_layer_num, ap16_range, maxQ_array, add_offset, fout);

	fclose(fout);
	free(tmp_buf_reorg);

	fout = fopen("bias_ap16_offset_add.bin","wb");
    	if(!fout) printf("fopen bias_ap16_offset_add.bin error\n");
	fwrite(add_offset, sizeof(int), conv_layer_num, fout);
	fclose(fout);

	fout = fopen("bias_ap16_Q.bin","wb");
    	if(!fout) printf("fopen bias_ap16_Q.bin error\n");
	fwrite(maxQ_array, sizeof(int), conv_layer_num, fout);
	fclose(fout);		

	free(inout_fixed_buf);
////////////////////////////////////////////////////#endif
	free(Memory_buf);
	free(Weight_buf);
	free(Beta_buf);

}
///////////////////////////////////////////////////////////////////////20181229 anti-reorg ok end n4m32
