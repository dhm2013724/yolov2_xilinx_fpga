
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>
#include "yolov2.h"

#include "acc.h"
#include "reorder_f32.h"

uint32_t generate_iofm_offset(uint32_t in_ptr[32], uint32_t out_ptr[32], network *net, uint32_t *iofm_size_max)
{
	const uint32_t MEM_LEN = (416*416*32+208*208*32);
	const uint32_t ROUTE16_LEN = (26*26*512);
	const uint32_t CONV27_LEN = (13*13*256);
	const uint32_t CONV24_LEN = (13*13*1024);
	const uint32_t OUT_LEN = (13*13*425);

	uint32_t iofm_max = ROUTE16_LEN;
	uint32_t top = 0;
	uint32_t bottom = top + MEM_LEN;
	for(int x=0;x<18;x++)
	{
		layer l = net->layers[x];
		uint16_t ifm_num = l.c;
		uint16_t ifm_w = l.w;
		uint16_t ifm_h = l.h;
		uint16_t ofm_num = l.n;
		if(l.type == MAXPOOL)
			ofm_num = l.c;
		uint16_t ofm_w = l.out_w;
		uint16_t ofm_h = l.out_h;

		uint32_t ifm_offset = ifm_num*ifm_h*ifm_w;
		uint32_t ofm_offset = ofm_num*ofm_h*ofm_w;
		if(ofm_offset > iofm_max)
			iofm_max = ofm_offset;
		if(ifm_offset > iofm_max)
			iofm_max = ifm_offset;

		if(x%2==0)
		{
			in_ptr[x] = top;
			out_ptr[x] = bottom - ofm_offset;
		}
		else
		{
			in_ptr[x] = out_ptr[x-1];
			out_ptr[x] = top;
		}
	}

	for(int x=18;x<25;x++)
	{
		layer l = net->layers[x];
		uint16_t ofm_num = l.n;
		if(l.type == MAXPOOL)
			ofm_num = l.c;
		uint16_t ofm_w = l.out_w;
		uint16_t ofm_h = l.out_h;
		uint16_t ifm_num = l.c;
		uint16_t ifm_w = l.w;
		uint16_t ifm_h = l.h;

		uint32_t ifm_offset = ifm_num*ifm_h*ifm_w;
		uint32_t ofm_offset = ofm_num*ofm_h*ofm_w;
		if(ofm_offset > iofm_max)
			iofm_max = ofm_offset;
		if(ifm_offset > iofm_max)
			iofm_max = ifm_offset;

		if(x%2==0)
		{
			in_ptr[x] = top;
			out_ptr[x] = bottom - ROUTE16_LEN - ofm_offset;
		}else
		{
			in_ptr[x] = out_ptr[x-1];
			out_ptr[x] = top;
		}
	}

	in_ptr[26] = bottom - ROUTE16_LEN;
	out_ptr[26] = top;

	in_ptr[27] = top;
	out_ptr[27] = bottom - ROUTE16_LEN - CONV24_LEN - CONV27_LEN;

	in_ptr[29] = out_ptr[27];
	out_ptr[29] = top;

	in_ptr[30] = top;
	out_ptr[30] = bottom - OUT_LEN;

	in_ptr[31] = out_ptr[30];

	*iofm_size_max = iofm_max;

	return MEM_LEN;
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

#define LNUM (32)
#define CONV_LNUM (23)

void yolov2_hls_ps(network *net, float *input)
{
	int weight_offset[32] = {864, 18432, 73728, 8192, 73728,
		294912, 32768, 294912, 1179648, 131072, 1179648, 131072,
		1179648, 4718592, 524288, 4718592, 524288, 4718592, 9437184,
		9437184, 32768, 11796480, 435200, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	int beta_offset[32] = {32, 64, 128, 64, 128, 256, 128, 256, 512, 256, 512, 256, 512, 1024,
		512, 1024, 512, 1024, 1024, 1024, 64, 1024, 425, 0, 0, 0, 0, 0, 0, 0, 0, 0};

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

	uint32_t in_ptr[32];
	uint32_t out_ptr[32];
	uint32_t iofm_size_max = 0;
	uint32_t fm_mem_size = generate_iofm_offset( in_ptr, out_ptr, net, &iofm_size_max);
	float *Memory_buf = (float*)calloc(fm_mem_size, sizeof(float));

	memcpy(Memory_buf + in_ptr[0], input, 416*416*3*sizeof(float));//416x416x3 input_pic

// #ifdef REORDER_GEN
	int16_t ap16_min = 0x8000;
	int16_t ap16_max = 0x7fff;
	float ap16_range[16*2];
	for(int i=0;i<16;i++)
	{
			//printf("Q%2d:",i);
			ap16_range[2*i]   = (float)ap16_min*pow((float)2,-i);//min
			ap16_range[2*i+1] = (float)ap16_max*pow((float)2,-i);//max
			// printf("min=%.7lf,max=%.7lf\n",ap16_range[2*i],ap16_range[2*i+1]);
	}
	int bias_Qarray[LNUM];
	int weight_Qarray[LNUM];
	int q_offset[LNUM];
	int ifm_Qarray[LNUM];
	int ofm_Qarray[LNUM];
	int inter_Qarray[LNUM];

	memset(inter_Qarray, 0, sizeof(inter_Qarray));
	FILE *fout;

	fout = fopen("bias_i16c.bin","wb");
	if(!fout) printf("fopen bias_f32c.bin error\n");

	memset(q_offset, 0, sizeof(q_offset));
	memset(bias_Qarray, 0, sizeof(bias_Qarray));	

	bias_reorder_i16c(Beta_buf, beta_offset, CONV_LNUM, ap16_range, bias_Qarray, LANE_NUM, q_offset, fout);
	fclose(fout);

	fout = fopen("bias_i16c_oadd.bin","wb");
	if(!fout) printf("fopen bias_i16c_oadd.bin error\n");
	fwrite(q_offset, sizeof(int), LNUM, fout);
	fclose(fout);

	fout = fopen("bias_i16c_Q.bin","wb");
	if(!fout) printf("fopen bias_i16c_Q.bin error\n");
	fwrite(bias_Qarray, sizeof(int), LNUM, fout);
	fclose(fout);	

	memset(q_offset, 0, sizeof(q_offset));

	fout = fopen("kernel_w_rc_i16.bin","wb");
	if(!fout) printf("fopen kernel_w_rc_i16.bin error\n");

	float *inout_fixed_buf = (float *)calloc(sizeof(float), iofm_size_max);
	printf("iofm_size_max=%d\n", iofm_size_max);	
// #endif

	int offset_index = 0;
	int woffset = 0;
	int boffset = 0;
	double sum_gop = 0.0;
	uint16_t ofm_w, ofm_h;
	float *ifm_ptr, *ofm_ptr;

    for(int i = 0; i < net->n; ++i)
	{
        layer l = net->layers[i];
		printf("Layer[%2d]: ",i);

		ifm_ptr = Memory_buf + in_ptr[i];
		ofm_ptr = Memory_buf + out_ptr[i];		

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

			// uint16_t TR,TC,TM,TN;
			// TC = MIN_diy(((IB_HxW-kernel_size)/kernel_stride+1),ofm_w);
			// uint16_t TCol = (TC-1)*kernel_stride + kernel_size;
			// TR = MIN_diy(((IB_HxW/TCol-kernel_size)/kernel_stride+1),ofm_h);//keep Kernel_stride>=1
			// TR = MIN_diy(TR, TrxTc/TC);
			// uint16_t TRow = (TR-1)*kernel_stride + kernel_size;
			// assert(((TR*TC)>0)&&((TR*TC)<=TrxTc));
			// assert(((TRow*TCol)>0)&&((TRow*TCol)<=IB_HxW));
			// // printf("TR=%d, TC=%d, TRow=%d, TCol=%d\n", TR, TC, TRow, TCol);

			uint16_t TR,TC,TM,TN;
			uint16_t TRow, TCol;
			assert((IB_HxW/l.size)>=l.size);
			TC = MIN_diy(((IB_HxW/l.size-l.size)/l.stride+1),ofm_w);
			TC = MIN_diy(TrxTc, TC);
			TCol = (TC-1)*l.stride + l.size;
			TR = MIN_diy(((IB_HxW/TCol-l.size)/l.stride+1),ofm_h);//keep Kernel_stride>=1
			TR = MIN_diy(TR, TrxTc/TC);
			TRow = (TR-1)*l.stride + l.size;				

			TM = MIN_diy(ofm_num,Tm);
			TN = MIN_diy(ifm_num,Tn);

			if(ltype != LT_CONV)
			{
				TM = MIN_diy(TM,TN);
				TN = 0;
			}
			bool IsNotConv = (ltype != LT_CONV);

			uint8_t  KK = kernel_size*kernel_size;
			uint32_t INumxKK = ifm_num*KK;//24bit
			uint32_t TRTC = (TR << 16) | TC;
			uint32_t TMTN = (TM << 16) | TN;
			uint32_t ofm_w_h = (ofm_w << 16) | ofm_h;
			uint32_t ifm_w_h = (ifm_w << 16) | ifm_h;
			uint32_t iofm_num = (ifm_num << 16) | ofm_num;
			uint32_t k_s_pad_ltype = (kernel_size << 24) | (kernel_stride << 16) | (pad_int << 8) | (ltype);
			uint32_t IHW = ifm_h*ifm_w;
			uint32_t OHW = ofm_h*ofm_w;
			uint32_t KK_INumxKK = (KK << 24) | INumxKK;

			uint32_t en_bits = 0x0;//enable_bits[2:0]={IsReLU, LoadBias, IsNotConv}
			if(IsNotConv)
				en_bits |= 0x1;
			if((ltype == LT_CONV) || (ltype == LT_DCONV))
				en_bits |= (0x1 << 1);

			int NToy = ceil(ofm_h*1.0f/TR);
			int NTox = ceil(ofm_w*1.0f/TC);
			int NTof = ceil(ofm_num*1.0f/TM);
			int NTcomb = NToy*NTox*NTof;

			int NTif;
			if(ltype == LT_CONV){
				NTif = ceil(ifm_num*1.0f/TN);
			}else{
				NTif = 1;
			}

			uint8_t lmode;
			int NTcomb_l;
			if(NTif==1){
				lmode = 0;
				NTcomb_l = NTcomb+2;
			}else{
				lmode = 1;
				NTcomb_l = NTcomb+1;
			}

			uint8_t weightQ, biasQ, ifmQ, ofmQ, avgQ, interQ;
			bool interQ_en = 0;
			avgQ = 0;
			// interQ = 12;					

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

		// #ifdef REORDER_GEN
					weight_reorder_i16c(Weight_buf + woffset,ifm_num,ofm_num,kernel_size,TM,TN, ltype, ap16_range, &weight_Qarray[offset_index], LANE_NUM, &q_offset[offset_index], fout);

					quantize_ifm_int16(ifm_ptr,inout_fixed_buf, ifm_h*ifm_w*ifm_num, ap16_range, &ifm_Qarray[i]);
					memcpy(ifm_ptr,inout_fixed_buf, ifm_h*ifm_w*ifm_num*sizeof(float));
		// #endif
					interQ_en = 0;	
					FPGA_Acc(ifm_ptr, ofm_ptr, Weight_buf+woffset, Beta_buf+boffset,
						k_s_pad_ltype, iofm_num, ifm_w_h, ofm_w_h, TRTC, TMTN, NToy, NTox, NTof, NTcomb, NTif, lmode, NTcomb_l, 0.0f, 1.0f,
						IHW, OHW, KK_INumxKK, en_bits, weightQ, biasQ, ifmQ, ofmQ, avgQ, interQ, interQ_en);

		// #ifdef REORDER_GEN
					quantize_ifm_int16(ofm_ptr,inout_fixed_buf, ofm_h*ofm_w*ofm_num, ap16_range, &ofm_Qarray[i]);
					memcpy(ofm_ptr,inout_fixed_buf, ofm_h*ofm_w*ofm_num*sizeof(float));

					interQ_en = 1;
					weightQ = weight_Qarray[offset_index]; biasQ = bias_Qarray[offset_index]; ifmQ = ifm_Qarray[i];
					ofmQ = ofm_Qarray[i];

					for(int tmpQ=0; tmpQ < INTERQ_MAX; tmpQ++){
						FPGA_Acc(ifm_ptr, inout_fixed_buf, Weight_buf + woffset, Beta_buf + boffset,
							k_s_pad_ltype, iofm_num, ifm_w_h, ofm_w_h, TRTC, TMTN, NToy, NTox, NTof, NTcomb, NTif, lmode, NTcomb_l, 0.0f, 1.0f,
							IHW, OHW, KK_INumxKK, en_bits, weightQ, biasQ, ifmQ, ofmQ, avgQ, tmpQ, interQ_en);
						
						interQ = diff_float_set(ofm_ptr, inout_fixed_buf, ofm_h*ofm_w*ofm_num, 1e-5, tmpQ, tmpQ==0);
					}
					inter_Qarray[i] = interQ;
					printf("weightQ=%d, biasQ=%d, ifmQ=%d, ofmQ=%d, avgQ=%d, interQ=%d\n", weightQ, biasQ, ifmQ, ofmQ, avgQ, interQ);
		// #endif							
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
		// #ifdef REORDER_GEN
					quantize_ifm_int16(ifm_ptr,inout_fixed_buf, ifm_h*ifm_w*ifm_num, ap16_range, &ifm_Qarray[i]);
					memcpy(ifm_ptr,inout_fixed_buf, ifm_h*ifm_w*ifm_num*sizeof(float));
		// #endif			
					interQ_en = 0;
					FPGA_Acc(ifm_ptr, ofm_ptr, NULL, NULL,
						k_s_pad_ltype, iofm_num, ifm_w_h, ofm_w_h, TRTC, TMTN, NToy, NTox, NTof, NTcomb, NTif, lmode, NTcomb_l, MIN_NEG, 1.0f,
						IHW, OHW, KK_INumxKK, en_bits, weightQ, biasQ, ifmQ, ofmQ, avgQ, interQ, interQ_en);
		// #ifdef REORDER_GEN
					quantize_ifm_int16(ofm_ptr,inout_fixed_buf, ofm_h*ofm_w*ofm_num, ap16_range, &ofm_Qarray[i]);
					memcpy(ofm_ptr,inout_fixed_buf, ofm_h*ofm_w*ofm_num*sizeof(float));
		// #endif
					avgQ = 15; interQ_en = 1;
					weightQ = 0; biasQ = 0; ifmQ = ifm_Qarray[i];
					ofmQ = ofm_Qarray[i]; 

					for(int tmpQ=0; tmpQ < INTERQ_MAX; tmpQ++){
						FPGA_Acc(ifm_ptr, inout_fixed_buf, NULL, NULL,
							k_s_pad_ltype, iofm_num, ifm_w_h, ofm_w_h, TRTC, TMTN, NToy, NTox, NTof, NTcomb, NTif, lmode, NTcomb_l, MIN_NEG, 1.0f,
							IHW, OHW, KK_INumxKK, en_bits, weightQ, biasQ, ifmQ, ofmQ, avgQ, tmpQ, interQ_en);
						
						interQ = diff_float_set(ofm_ptr, inout_fixed_buf, ofm_h*ofm_w*ofm_num, 1e-5, tmpQ, tmpQ==0);
					}
					inter_Qarray[i] = interQ;
					printf("weightQ=%d, biasQ=%d, ifmQ=%d, ofmQ=%d, avgQ=%d, interQ=%d\n", weightQ, biasQ, ifmQ, ofmQ, avgQ, interQ);

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
					reorg_cpu(ifm_ptr, ofm_w, ofm_h, 4, 2, ofm_ptr);
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
					forward_region_layer(l, ifm_ptr);
					break;
			}
		}
    	}
	printf("SUM_GOP=%.7lf\n",sum_gop);

// #ifdef REORDER_GEN
	fclose(fout);//weight reorg quanti
	free(inout_fixed_buf);

	FILE *fp_ioq;
	fp_ioq = fopen("i16c_inter_Q.bin", "wb");
	if(!fp_ioq) printf("i16c_inter_Q.bin\n");
	fwrite(inter_Qarray, sizeof(int), LNUM, fp_ioq);
	fclose(fp_ioq);//interQ 	

	fp_ioq = fopen("i16c_ifm_Q.bin", "wb");
	if(!fp_ioq) printf("i16c_ifm_Q.bin\n");
	fwrite(ifm_Qarray, sizeof(int), LNUM, fp_ioq);
	fclose(fp_ioq);//iofm quanti

	fp_ioq = fopen("i16c_ofm_Q.bin", "wb");
	if(!fp_ioq) printf("i16c_ofm_Q.bin\n");
	fwrite(ofm_Qarray, sizeof(int), LNUM, fp_ioq);
	fclose(fp_ioq);//iofm quanti		

	fout = fopen("kernel_w_i16c_oadd.bin","wb");
	if(!fout) printf("fopen kernel_w_i16c_oadd.bin error\n");
	fwrite(q_offset, sizeof(int), LNUM, fout);
	fclose(fout);

	fout = fopen("kernel_w_i16c_Q.bin","wb");
	if(!fout) printf("fopen kernel_w_i16c_Q.bin error\n");
	fwrite(weight_Qarray, sizeof(int), LNUM, fout);
	fclose(fout);	
// #endif

	free(Memory_buf);
	free(Weight_buf);
	free(Beta_buf);

}
///////////////////////////////////////////////////////////////////////20181229 anti-reorg ok end n4m32
