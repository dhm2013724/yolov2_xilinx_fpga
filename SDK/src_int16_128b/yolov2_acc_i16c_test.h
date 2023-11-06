
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>
#include <sys/stat.h>
#include <unistd.h>
#include "yolov2.h"

#include "acc_i16c.h"
//#include "reorder_f32.h"

int read_binfile_flt32_rb(float *buf, char *filename, int data_num)
{
	FILE *fp;
	if( (fp = fopen(filename, "rb")) == NULL)
		printf("cannot open bin_file %s\n", filename);
	int rd_num  = fread(buf, sizeof(unsigned char)*4, data_num, fp);
	fclose(fp);

	return rd_num;
}

unsigned long get_file_size(const char *filename)  
{  
    struct stat buf;  
    if(stat(filename, &buf)<0)  
    {  
        return 0;  
    }
	printf("%s's data size is %ld\n", filename, buf.st_size);  
    return (unsigned long)buf.st_size;  
}

uint32_t generate_iofm_offset(uint32_t in_ptr[32], uint32_t out_ptr[32], network *net)
{
	// const uint32_t MEM_LEN = ((416*416+208*208)*((32 + LANE_NUM -1)/LANE_NUM)*LANE_NUM);
	const uint32_t ROUTE16_LEN = (26*26*((512 + LANE_NUM -1)/LANE_NUM)*LANE_NUM);
	const uint32_t CONV27_LEN = (13*13*((256 + LANE_NUM -1)/LANE_NUM)*LANE_NUM);
	const uint32_t CONV24_LEN = (13*13*((1024 + LANE_NUM -1)/LANE_NUM)*LANE_NUM);
	const uint32_t OUT_LEN = (13*13*((425 + LANE_NUM -1)/LANE_NUM)*LANE_NUM);

	uint32_t top_max = 0, bot_max = 0;
	for(int x=0;x<8;x++)
	{
		layer l = net->layers[x*2];
		uint16_t ifm_num = l.c;
		uint16_t ifm_w = l.w;
		uint16_t ifm_h = l.h;
		uint16_t ofm_num = l.n;
		if(l.type == MAXPOOL)
			ofm_num = l.c;
		uint16_t ofm_w = l.out_w;
		uint16_t ofm_h = l.out_h;

		uint32_t ifm_offset = ((ifm_num + LANE_NUM -1)/LANE_NUM)*LANE_NUM*ifm_h*ifm_w;
		uint32_t ofm_offset = ((ofm_num + LANE_NUM -1)/LANE_NUM)*LANE_NUM*ofm_h*ofm_w;

		if(ifm_offset > top_max)
			top_max = ifm_offset;

		if(ofm_offset > bot_max)
			bot_max = ofm_offset;
	}

	uint32_t ifm_list[6] = {16, 19, 21, 23, 27, 30};
	for(int x=0;x<6;x++)
	{
		uint32_t lnum = ifm_list[x];
		layer l = net->layers[lnum];
		uint16_t ifm_num = l.c;
		uint16_t ifm_w = l.w;
		uint16_t ifm_h = l.h;

		uint32_t ifm_offset = ((ifm_num + LANE_NUM -1)/LANE_NUM)*LANE_NUM*ifm_h*ifm_w;

		if(ifm_offset > top_max)
			top_max = ifm_offset;
	}

	uint32_t ofm_list[4] = {18, 20, 22, 24};
	for(int x=0;x<4;x++)
	{
		uint32_t lnum = ofm_list[x];
		layer l = net->layers[lnum];
		uint16_t ifm_num = l.c;
		uint16_t ifm_w = l.w;
		uint16_t ifm_h = l.h;

		uint32_t ifm_offset = ((ifm_num + LANE_NUM -1)/LANE_NUM)*LANE_NUM*ifm_h*ifm_w + ROUTE16_LEN;

		if(ifm_offset > bot_max)
			bot_max = ifm_offset;
	}

	if((ROUTE16_LEN+CONV24_LEN+CONV27_LEN)> bot_max)
		bot_max = ROUTE16_LEN+CONV24_LEN+CONV27_LEN;

	printf("top_max = %d, bot_max = %d\n", top_max, bot_max);

	uint32_t MEM_LEN = top_max + bot_max;
	uint32_t top = 0;
	uint32_t bottom = top + MEM_LEN;
	for(int x=0;x<18;x++)
	{
		layer l = net->layers[x];
		uint16_t ofm_num = l.n;
		if(l.type == MAXPOOL)
			ofm_num = l.c;
		uint16_t ofm_w = l.out_w;
		uint16_t ofm_h = l.out_h;

		uint32_t ofm_offset = ((ofm_num + LANE_NUM -1)/LANE_NUM)*LANE_NUM*ofm_h*ofm_w;

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

		uint32_t ofm_offset = ((ofm_num + LANE_NUM -1)/LANE_NUM)*LANE_NUM*ofm_h*ofm_w;

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

	return MEM_LEN;
}

template<typename DT>
void reorg_cpu(DT *x, int w, int h, int c, int stride, DT *out)
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

	double time1,time2;
	time1 = what_time_is_it_now();
	copy_file2mem("kernel_w_rc_i16.bin", get_file_size("kernel_w_rc_i16.bin"), WEIGHT_BASE);
	printf("yolov2_w copy ok\n");
	copy_file2mem("bias_i16c.bin", get_file_size("bias_i16c.bin"), BETA_BASE);
	printf("yolov2_b copy ok\n");
	time2 = what_time_is_it_now();
	printf("Predicted in %f seconds.\n",time2 - time1);

	int w_Q[LNUM];
	int w_aoffset[LNUM];

	printf("%d\n", read_binfile_flt32_rb((float*)w_Q, "kernel_w_i16c_Q.bin", LNUM));
	printf("%d\n", read_binfile_flt32_rb((float*)w_aoffset, "kernel_w_i16c_oadd.bin", LNUM));	

	//bias_i16 bias_i16_Q bias_i16_offset_add
	int bias_Q[LNUM];
	int bias_aoffset[LNUM];	

	printf("%d\n", read_binfile_flt32_rb((float*)bias_Q, "bias_i16c_Q.bin", LNUM));
	printf("%d\n", read_binfile_flt32_rb((float*)bias_aoffset, "bias_i16c_oadd.bin", LNUM));

	int inter_Q[LNUM];
	printf("%d\n", read_binfile_flt32_rb((float*)inter_Q, "i16c_inter_Q.bin", LNUM));

	//i16_iofm_Q
	int IFM_Q[LNUM], OFM_Q[LNUM];
	printf("%d\n", read_binfile_flt32_rb((float*)IFM_Q, "i16c_ifm_Q.bin", LNUM));
	printf("%d\n", read_binfile_flt32_rb((float*)OFM_Q, "i16c_ofm_Q.bin", LNUM));

	if(OFM_Q[24] < OFM_Q[26])
		OFM_Q[26] = OFM_Q[24];
	else
		OFM_Q[24] = OFM_Q[26];
	// IFM_Q[17] = IFM_Q[26] = OFM_Q[16];
	IFM_Q[29] = OFM_Q[24];		

	uint32_t in_ptr[32];
	uint32_t out_ptr[32];
	uint32_t fm_mem_size = generate_iofm_offset( in_ptr, out_ptr, net);
//	int16_t *Memory_buf = (int16_t*)calloc(fm_mem_size, sizeof(int16_t));

	float *region_buffer = (float *)calloc(416*416*4,sizeof(float));
	if(!region_buffer) file_error("region_buffer error \n");
	float *region_buffer2 = (float *)calloc(13*13*432,sizeof(float));
	if(!region_buffer2) file_error("region_buffer error \n");

	// memcpy(Memory_buf + in_ptr[0], input, 416*416*3*sizeof(float));//416x416x3 input_pic
	int img_c = 3, img_h = 416, img_w =416;
//	int16_t *ifm_buf = Memory_buf + in_ptr[0];
	int16_t *ifm_buf = (int16_t *)region_buffer;
	int FM_Q = IFM_Q[0];
	printf("IFM_Q[0]=%d\n", FM_Q);
	float tmp_exp = pow(2.0,FM_Q);	
//	float* tmp_ptr_ft32 = input;
	int Ca3_d4 =  (img_c+LANE_NUM-1)/LANE_NUM;
	for(int c4=0; c4<Ca3_d4; c4++)
	for(int ih=0; ih<img_h; ih++)
	for(int iw=0; iw<img_w; iw++)
	for(int ic=0; ic<LANE_NUM; ic++)
	{
		int idx_out = c4*img_h*img_w*LANE_NUM + ih*img_w*LANE_NUM + iw*LANE_NUM + ic;
		int idx_in = (c4*LANE_NUM+ic)*img_h*img_w + ih*img_w + iw;

		if((c4*LANE_NUM+ic) >= img_c)
			ifm_buf[idx_out] = 0;
		else
			ifm_buf[idx_out] = input[idx_in]*tmp_exp;
	}
	copy_mem2dev((uint8_t *)ifm_buf, 416*416*Ca3_d4*LANE_NUM*2, MEM_BASE + in_ptr[0]*2);

	int16_t *tmp_i16_ptr0, *tmp_i16_ptr1;
	int offset_index = 0;
	int woffset = 0;
	int boffset = 0;
	double sum_gop = 0.0;
	uint16_t ofm_w, ofm_h;
//	int16_t *ifm_ptr, *ofm_ptr, *tmp_i16_ptr;
	unsigned int ifm_ptr, ofm_ptr;
//	float *region_buf = (float *)malloc(13*13*425*sizeof(float));
//	int16_t *reorg_buf = (int16_t *)malloc(26*26*64*sizeof(int16_t));

	time1 = what_time_is_it_now();
    for(int i = 0; i < net->n; ++i)
	{
        layer l = net->layers[i];
		printf("Layer[%2d]: ",i);

		ifm_ptr = MEM_BASE + in_ptr[i]*2;
		ofm_ptr = MEM_BASE + out_ptr[i]*2;
		uint32_t out_align;
		uint32_t out_left_num;

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

			int biasQ = bias_Q[offset_index];
			int weightQ = w_Q[offset_index];
			int ifmQ = IFM_Q[i];
			int ofmQ = OFM_Q[i];
			int avgQ = 0;
			int interQ = inter_Q[i];	

			if(!((ltype == LT_CONV) || (ltype == LT_DCONV))){
				weightQ = 0;
				biasQ = 0;
			}					

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

					FPGA_Acc(ifm_ptr, ofm_ptr, woffset, boffset,
						k_s_pad_ltype, iofm_num, ifm_w_h, ofm_w_h, TRTC, TMTN, NToy, NTox, NTof, NTcomb, NTif, lmode, NTcomb_l, 0, 1,
						IHW, OHW, KK_INumxKK, en_bits, weightQ, biasQ, ifmQ, ofmQ, avgQ, interQ);						

					woffset = woffset + weight_offset[offset_index] + w_aoffset[offset_index];
					boffset = boffset + beta_offset[offset_index] + bias_aoffset[offset_index];
					offset_index++;

					break;
				case MAXPOOL:
					printf("outputMemory:%8d;max          %d x %d / %d  %4d x%4d x%4d   ->  %4d x%4d x%4d\n",l.outputs, l.size, l.size, l.stride, l.w, l.h, l.c,
						 l.out_w, l.out_h, l.out_c);
					sum_gop += (4.0 * l.c * l.out_h*l.out_w)/1000000000.;
					//output_w = (l.w - l.size)/l.stride + 1 ;
					//output_h = (l.h - l.size)/l.stride + 1 ;
					//very strange that yolo's maxpool padding is all 0?!!
					FPGA_Acc(ifm_ptr, ofm_ptr, NULL, NULL,
						k_s_pad_ltype, iofm_num, ifm_w_h, ofm_w_h, TRTC, TMTN, NToy, NTox, NTof, NTcomb, NTif, lmode, NTcomb_l, 0x8000, 1,
						IHW, OHW, KK_INumxKK, en_bits, weightQ, biasQ, ifmQ, ofmQ, avgQ, interQ);
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

					tmp_i16_ptr0 = (int16_t *)region_buffer;
					tmp_i16_ptr1 = (int16_t *)region_buffer2;
					copy_dev2mem((uint8_t *)tmp_i16_ptr0, 26*26*64*2, ifm_ptr);
					for(int oc=0; oc<64;oc++)
					for(int oh=0; oh<26;oh++)
					for(int ow=0; ow<26;ow++){
						tmp_i16_ptr1[oc*26*26 + oh*26 + ow] = tmp_i16_ptr0[((oc/LANE_NUM)*26*26 + oh*26 + ow)*LANE_NUM + (oc & (LANE_NUM-1))];
					}

					reorg_cpu(tmp_i16_ptr1, ofm_w, ofm_h, 4, 2, tmp_i16_ptr0);

					out_align = ofm_ptr & 0xFFFFF000;
					out_left_num = ofm_ptr & 0xFFF;

					for(int oc=0; oc<256;oc++)
					for(int oh=0; oh<13;oh++)
					for(int ow=0; ow<13;ow++){
						tmp_i16_ptr1[((oc/LANE_NUM)*13*13 + oh*13 + ow)*LANE_NUM + (oc & (LANE_NUM-1)) + (out_left_num/2)] = tmp_i16_ptr0[oc*13*13 + oh*13 + ow];
					}
//					copy_mem2dev((uint8_t *)tmp_i16_ptr1, 13*13*256*2, ofm_ptr);
					copy_mem2dev((uint8_t *)tmp_i16_ptr1, out_left_num + 13*13*256*2, out_align);

//					tmp_i16_ptr = (int16_t *)region_buf;
//
//					for(int oc=0; oc<64;oc++)
//					for(int oh=0; oh<26;oh++)
//					for(int ow=0; ow<26;ow++){
//						reorg_buf[oc*26*26 + oh*26 + ow] = ifm_ptr[((oc/LANE_NUM)*26*26 + oh*26 + ow)*LANE_NUM + (oc & (LANE_NUM-1))];
//					}
//
//					reorg_cpu(reorg_buf, ofm_w, ofm_h, 4, 2, tmp_i16_ptr);
//
//					for(int oc=0; oc<256;oc++)
//					for(int oh=0; oh<13;oh++)
//					for(int ow=0; ow<13;ow++){
//						ofm_ptr[((oc/LANE_NUM)*13*13 + oh*13 + ow)*LANE_NUM + (oc & (LANE_NUM-1))] = tmp_i16_ptr[oc*13*13 + oh*13 + ow];
//					}

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
					tmp_exp = pow(2.0, -OFM_Q[30]);
//					for(int oc=0; oc<425;oc++)
//					for(int oh=0; oh<13;oh++)
//					for(int ow=0; ow<13;ow++){
//						region_buf[oc*13*13 + oh*13 + ow] = ifm_ptr[((oc/LANE_NUM)*13*13 + oh*13 + ow)*LANE_NUM + (oc & (LANE_NUM-1))]*tmp_exp;
//					}
//					forward_region_layer(l, region_buf);

					tmp_i16_ptr0 = (int16_t *)region_buffer;
//					tmp_i16_ptr1 = (int16_t *)region_buffer2;
					out_align = ifm_ptr & 0xFFFFF000;
					out_left_num = ifm_ptr & 0xFFF;

					copy_dev2mem((uint8_t *)tmp_i16_ptr0, 13*13*(425+LANE_NUM-1)/LANE_NUM*LANE_NUM*2 + out_left_num, out_align);
					for(int oc=0; oc<425;oc++)
					for(int oh=0; oh<13;oh++)
					for(int ow=0; ow<13;ow++){
						region_buffer2[oc*13*13 + oh*13 + ow] =
				tmp_i16_ptr0[((oc/LANE_NUM)*13*13 + oh*13 + ow)*LANE_NUM + (oc & (LANE_NUM-1)) + (out_left_num/2)]*tmp_exp;
					}

					forward_region_layer(l, region_buffer2);

					// forward_region_layer(l, ifm_ptr);
					break;
			}
		}
    	}
	time2 = what_time_is_it_now();
	printf("Inference in %f seconds.(+region)\n",time2 - time1);
	printf("SUM_GOP=%.7lf\n",sum_gop);

//	free(region_buf);
//	free(reorg_buf);
//	free(Memory_buf);
//	free(weight_buf);
//	free(beta_buf);

	free(region_buffer);
	free(region_buffer2);

}
///////////////////////////////////////////////////////////////////////20181229 anti-reorg ok end n4m32

// 845
// [248]:h=0.227699,w=0.232511,x=0.127223,y=0.502272,objectness=0.846463
// [266]:h=0.231706,w=0.103453,x=0.511701,y=0.636352,objectness=0.683597
// tvmonitor: 85%
// chair: 67%
// YOLOv2 TEST End
