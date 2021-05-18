
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/time.h>

double what_time_is_it_now()
{
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

unsigned long get_file_size(const char *filename)  
{  
    struct stat buf;  
    if(stat(filename, &buf)<0)  
    {  
        return 0;  
    }  
    return (unsigned long)buf.st_size;  
}

//#define INTERWIDTH 19
#include "cnn.h"

//256b
#define MEM_LEN (416*416*32/2*4+208*208*32/2*4)
void generate_iofm_offset(uint32_t in_ptr[32], uint32_t out_ptr[32], uint32_t Memory_buf, network *net)
{
#define ROUTE16_LEN (26*26*512/2*4)
#define CONV27_LEN (13*14*256/2*4)
#define CONV24_LEN (13*14*1024/2*4)

	uint32_t Memory_top = Memory_buf;
	uint32_t Memory_bottom = Memory_top + MEM_LEN;
	int x;
	for(x=0;x<18;x++)
	{
		int out_w = net->layers[x].out_w;
		int out_w_align_256b = (out_w >> 1) << 1;
		if(out_w & 0x1)
			out_w_align_256b += 2;

		if(x%2==0)
		{
			in_ptr[x] = Memory_top;
			out_ptr[x] = Memory_bottom - net->layers[x].out_c * net->layers[x].out_h * out_w_align_256b/2*4;
		}
		else
		{
			in_ptr[x] = out_ptr[x-1];
			out_ptr[x] = Memory_top;
		}
	}

	for(x=18;x<25;x++)
	{
		int out_w = net->layers[x].out_w;
		int out_w_align_256b = (out_w >> 1) << 1;
		if(out_w & 0x1)
			out_w_align_256b += 2;

		if(x%2==0)
		{
			in_ptr[x] = Memory_top;
			out_ptr[x] = Memory_bottom - ROUTE16_LEN - net->layers[x].out_c * net->layers[x].out_h * out_w_align_256b/2*4;
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
	out_ptr[30] = Memory_bottom - (net->layers[30].outputs + 13*1*425)/2*4;//13*14*425/2

	in_ptr[31] = out_ptr[30];
}

void reorg_cpu(int16_t *x, int w, int h, int c, int stride, int16_t *out)
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

	uint32_t file_size = 0;

	file_size = get_file_size("iofm_Q.bin");
	int32_t *iofmQ = (int32_t *)malloc(file_size);
	if(!iofmQ) printf("iofmQ alloc failed.\n");
	printf("iofm_Q's size = %d\n", file_size);
	FILE *fp_ioq = fopen("iofm_Q.bin", "rb");
    	if(!fp_ioq) file_error("iofm_Q.bin");
	fread(iofmQ, sizeof(int32_t), file_size/4, fp_ioq);
	fclose(fp_ioq);
	if(iofmQ[20] < iofmQ[21])
		iofmQ[21] = iofmQ[20];
	else
		iofmQ[20] = iofmQ[21];

	FILE *fp_w, *fp_b;
	int weights_ap16_offset_add[32];
	int WeightQ[32];

	file_size = get_file_size("weight_reorg_ap16.bin");
	copy_file2mem("weight_reorg_ap16.bin", file_size, WEIGHT_BASEADDR);

	file_size = get_file_size("bias_ap16.bin");
	copy_file2mem("bias_ap16.bin", file_size, BETA_BASEADDR);

	file_size = get_file_size("weights_ap16_offset_add.bin");
	printf("weights_ap16_offset_add's size = %d\n", file_size);
	fp_w = fopen("weights_ap16_offset_add.bin", "rb");
    	if(!fp_w) file_error("weights_ap16_offset_add.bin");
	fread(weights_ap16_offset_add, sizeof(int32_t), file_size/4, fp_w);
	fclose(fp_w);

	file_size = get_file_size("weights_ap16_Q.bin");
	printf("weights_ap16_Q's size = %d\n", file_size);
	fp_w = fopen("weights_ap16_Q.bin", "rb");
    	if(!fp_w) file_error("weights_ap16_Q.bin");
	fread(WeightQ, sizeof(int32_t), file_size/4, fp_w);
	fclose(fp_w);

	int bias_ap16_offset_add[32];
	int BetaQ[32];

	file_size = get_file_size("bias_ap16_offset_add.bin");
	printf("bias_ap16_offset_add's size = %d\n", file_size);
	fp_b = fopen("bias_ap16_offset_add.bin", "rb");
    	if(!fp_b) file_error("bias_ap16_offset_add.bin");
	fread(bias_ap16_offset_add, sizeof(int32_t), file_size/4, fp_b);
	fclose(fp_b);

	file_size = get_file_size("bias_ap16_Q.bin");
	printf("bias_ap16_Q's size = %d\n", file_size);
	fp_b = fopen("bias_ap16_Q.bin", "rb");
    	if(!fp_b) file_error("bias_ap16_Q.bin");
	fread(BetaQ, sizeof(int32_t), file_size/4, fp_b);
	fclose(fp_b);

//	int32_t *Memory_buf = (int32_t*)calloc(MEM_LEN+512*2,sizeof(int32_t));
	uint32_t in_ptr[32];
	uint32_t out_ptr[32];
	generate_iofm_offset( in_ptr, out_ptr, MEM_BASEADDR, net);

	int16_t *region_buf = (int16_t *)calloc(416*416*3,sizeof(int16_t));
	float *region_buf2 = (float *)calloc(13*16*425,sizeof(float));

	double tmp_iofmQ_pow = pow(2.0,iofmQ[0]);
	int32_t* tmp_ptr_int32;
	int16_t* tmp_ptr_int16;
	//memcpy(in_ptr[0], input, 416*416*3*sizeof(float));//416x416x3 input_pic
	for(int i=0;i<416*416*3;i++)//1st Layer input Q14
	{
		region_buf[i] = (int16_t)(input[i]*tmp_iofmQ_pow);
	}
	copy_mem2dev((uint8_t *)region_buf, 2*416*416*3, in_ptr[0]);

	int offset_index = 0;
	int woffset = 0;
	int boffset = 0;
	double sum_gop = 0.0;
	uint16_t ofm_w;
	uint16_t ofm_h;
	int inputQ_idx;
	double first, second;
	uint32_t out_align;
	uint32_t out_left_num;

	first = what_time_is_it_now();
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

					inputQ_idx = offset_index;
					if(i==26)
					{
						inputQ_idx = 13;
					}

					FPGA_Acc(in_ptr[i], out_ptr[i], woffset, boffset,
						k_s_pad_ltype, iofm_num, ifm_w_h, ofm_w_h, TRTC, TMTN, OFM_num_bound, mLoopsxTM, mLoops_a1xTM, 0,
						TRowTCol, IHW, OHW, KK_INumxKK, en_bits, WeightQ[offset_index], BetaQ[offset_index], iofmQ[inputQ_idx], iofmQ[offset_index+1]);

					woffset += (weight_offset[offset_index] + weights_ap16_offset_add[offset_index])/2;
					boffset += (beta_offset[offset_index] + bias_ap16_offset_add[offset_index])/2;
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
						TRowTCol, IHW, OHW, KK_INumxKK, en_bits, INTERWIDTH, INTERWIDTH, iofmQ[offset_index], iofmQ[offset_index]);

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
					//reorg_cpu(in_ptr[i], ofm_w, ofm_h, 4, 2, out_ptr[i]);

					copy_dev2mem((uint8_t *)region_buf, 26*26*64*sizeof(int16_t), in_ptr[i]);
//					tmp_ptr_int32 = in_ptr[i];
//					memcpy((int16_t *)(region_buf), (int16_t *)(tmp_ptr_int32), 26*26*64*sizeof(int16_t));

				//	for(int k = 0; k<26*64; k++)
				//		memcpy((int16_t *)(region_buf) + k*26, (int16_t *)(tmp_ptr_int32) + k*32, 26*sizeof(int16_t));
					out_align = out_ptr[i] & 0xFFFFF000;
					out_left_num = out_ptr[i] & 0xFFF;

					tmp_ptr_int16 = (int16_t *)(region_buf2);
					reorg_cpu(region_buf, ofm_w, ofm_h, 4, 2, tmp_ptr_int16);
					for(int k = 0; k<13*256; k++)
						memcpy(region_buf + k*14 + (out_left_num/2), tmp_ptr_int16 + k*13, 13*sizeof(int16_t));
					copy_mem2dev((uint8_t *)region_buf, out_left_num + 13*14*256*sizeof(int16_t), out_align);
//					memcpy(out_ptr[i], (int32_t *)region_buf, 13*14*256*sizeof(int32_t)/2);
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

					out_align = in_ptr[i] & 0xFFFFF000;
					out_left_num = in_ptr[i] & 0xFFF;

					copy_dev2mem((uint8_t *)region_buf, out_left_num + 13*14*425*sizeof(int16_t), out_align);
					tmp_iofmQ_pow = pow(2.0, -iofmQ[offset_index]);
					for(int k = 0; k<13*425; k++)
						for(int j = 0; j < 13; j++)
						{
							region_buf2[k*13 + j] = region_buf[k*14 + j + (out_left_num/2)]*tmp_iofmQ_pow;
						}
					forward_region_layer(l, region_buf2);
					//forward_region_layer(l, in_ptr[i]);
					break;
			}
		}
    }
    second = what_time_is_it_now();
    printf("Forward time: %5.5lf seconds.\n", second - first);
	printf("SUM_GOP=%.7lf\n",sum_gop);

	free(region_buf);
	free(region_buf2);

//	free(Memory_buf);
//	free(Weight_buf);
//	free(Beta_buf);
}
///////////////////////////////////////////////////////////////////////20181229 anti-reorg ok end n4m32
