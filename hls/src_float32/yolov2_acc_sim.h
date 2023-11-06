

#include "cnn.h"

#include <sys/stat.h>
#include <unistd.h>

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
		int out_w = net->layers[x].out_w;
		if(x%2==0)
		{
			in_ptr[x] = Memory_top;
			out_ptr[x] = Memory_bottom - net->layers[x].out_c *  net->layers[x].out_h * out_w;
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
		if(x%2==0)
		{
			in_ptr[x] = Memory_top;
			out_ptr[x] = Memory_bottom - ROUTE16_LEN - net->layers[x].out_c *  net->layers[x].out_h * out_w;
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
	uint32_t weight_offset[32] = {864, 18432, 73728, 8192, 73728,
		294912, 32768, 294912, 1179648, 131072, 1179648, 131072,
		1179648, 4718592, 524288, 4718592, 524288, 4718592, 9437184,
		9437184, 32768, 11796480, 435200, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	uint32_t beta_offset[32] = {32, 64, 128, 64, 128, 256, 128, 256, 512, 256, 512, 256, 512, 1024,
		512, 1024, 512, 1024, 1024, 1024, 64, 1024, 425, 0, 0, 0, 0, 0, 0, 0, 0, 0};

	uint32_t data_num = 0;

	data_num = get_file_size("weights_reorg.bin")/4;
	float *Weight_buf = (float *)calloc(data_num,sizeof(float));
	FILE *fp_w = fopen("weights_reorg.bin", "rb");
    if(!fp_w) file_error("weights_reorg.bin");
	fread(Weight_buf, sizeof(float), data_num, fp_w);
	fclose(fp_w);

	data_num = get_file_size("bias.bin")/4;
	float *Beta_buf   = (float *)calloc(data_num,sizeof(float));
	FILE *fp_b = fopen("bias.bin", "rb");
    if(!fp_b) file_error("bias.bin");
	fread(Beta_buf, sizeof(float), data_num, fp_b);
	fclose(fp_b);

	float *region_buf = (float *)calloc(13*13*425, sizeof(float));
	float *tmp_ptr = NULL;

//leave some memories for overflow, because the load_module will load extra pixels near boundary for padding
	float *Memory_buf = (float*)calloc(MEM_LEN+512*2,sizeof(float));
	float* in_ptr[32];
	float* out_ptr[32];
	generate_iofm_offset( in_ptr, out_ptr, Memory_buf, net);

	memcpy(in_ptr[0], input, 416*416*3*sizeof(float));//416x416x3 input_pic

	int offset_index = 0;
	int woffset = 0;
	int boffset = 0;
	int TR,TC,TM,TN;
	int TRow, TCol;
	int output_w,output_h;
	int mLoops;
	int NToy, NTox, NTof, NTcomb, NTif;
	uint8_t lmode;
	int NTcomb_l;

    for(uint16_t i = 0; i < net->n; ++i)
	{
        layer l = net->layers[i];
		printf("Layer[%2d]: ",i);
		switch(l.type)
		{
			case CONVOLUTIONAL:
				printf("outputMemory:%8d;BN=%d;Activation=%d;conv  %5d %2d x%2d /%2d  %4d x%4d x%4d   ->  %4d x%4d x%4d  %5.3f BFLOPs\n",l.outputs,l.batch_normalize,l.activation, l.n, l.size, l.size, l.stride, l.w, l.h, l.c, l.out_w, l.out_h, l.out_c, (2.0 * l.n * l.size*l.size*l.c/l.groups * l.out_h*l.out_w)/1000000000.);

				output_w = (l.w - l.size + 2*l.pad)/l.stride + 1;
				output_h = (l.h - l.size + 2*l.pad)/l.stride + 1;

//				TR = MIN(((OnChipIB_Height-l.size)/l.stride+1),Tr);//keep Kstride>=1
//				TR = MIN(output_h,TR);
//				TC = MIN(((OnChipIB_Width-l.size)/l.stride+1),Tc);
//				TC = MIN(output_w,TC);

				assert((IB_HxW/l.size)>=l.size);
				TC = MIN(((IB_HxW/l.size-l.size)/l.stride+1),output_w);
				TC = MIN(TrxTc, TC);
				TCol = (TC-1)*l.stride + l.size;
				TR = MIN(((IB_HxW/TCol-l.size)/l.stride+1),output_h);//keep Kernel_stride>=1
				TR = MIN(TR, TrxTc/TC);
				TRow = (TR-1)*l.stride + l.size;

				TM = MIN(l.n,Tm);
				TN = MIN(l.c,Tn);

				NToy = ceil(output_h*1.0f/TR);
				NTox = ceil(output_w*1.0f/TC);
				NTof = ceil(l.n*1.0f/TM);
				NTcomb = NToy*NTox*NTof;
				NTif = ceil(l.c*1.0f/TN);

				if(NTif==1){
					lmode = 0;
					NTcomb_l = NTcomb+2;
				}else{
					lmode = 1;
					NTcomb_l = NTcomb+1;
				}

				YOLO2_FPGA(in_ptr[i],out_ptr[i],Weight_buf+woffset,Beta_buf+boffset,
					l.c,l.n,l.size,
					l.stride,l.w,l.h,output_w, output_h, l.pad,l.activation==LEAKY?1:0,
					TM,TN,TR,TC, NToy, NTox, NTof, NTcomb, NTif, lmode, NTcomb_l, 0);

				woffset += weight_offset[offset_index];
				boffset += beta_offset[offset_index];
				offset_index++;

				break;
			case MAXPOOL:
				printf("outputMemory:%8d;max          %d x %d / %d  %4d x%4d x%4d   ->  %4d x%4d x%4d\n",l.outputs, l.size, l.size, l.stride, l.w, l.h, l.c, l.out_w, l.out_h, l.out_c);
				//output_w = (l.w - l.size)/l.stride + 1 ;
				//output_h = (l.h - l.size)/l.stride + 1 ;
				output_w = l.out_h;
				output_h = l.out_w;

//				TR = MIN(((OnChipIB_Height-l.size)/l.stride+1),Tr);//keep Kstride>=1
//				TC = MIN(((OnChipIB_Width-l.size)/l.stride+1),Tc);
//				TR = MIN(output_h,TR);
//				TC = MIN(output_w,TC);

				assert((IB_HxW/l.size)>=l.size);
				TC = MIN(((IB_HxW/l.size-l.size)/l.stride+1),output_w);
				TC = MIN(TrxTc, TC);
				TCol = (TC-1)*l.stride + l.size;
				TR = MIN(((IB_HxW/TCol-l.size)/l.stride+1),output_h);//keep Kernel_stride>=1
				TR = MIN(TR, TrxTc/TC);
				TRow = (TR-1)*l.stride + l.size;

				TM = MIN(Tm,Tn);
				TM = MIN(l.c,TM);

				NToy = ceil(output_h*1.0f/TR);
				NTox = ceil(output_w*1.0f/TC);
				NTof = ceil(l.c*1.0f/TM);
				NTcomb = NToy*NTox*NTof;
				NTif = 1;

				if(NTif==1){
					lmode = 0;
					NTcomb_l = NTcomb+2;
				}else{
					lmode = 1;
					NTcomb_l = NTcomb+1;
				}

				YOLO2_FPGA(in_ptr[i],out_ptr[i],NULL,NULL,l.c,l.c,
					l.size,l.stride,l.w,l.h, output_w, output_h, l.pad,0,TM,0,TR,TC,
					NToy, NTox, NTof, NTcomb, NTif, lmode, NTcomb_l, 1);

				break;
			case REORG:
				printf("outputMemory:%8d;reorg              /%2d  %4d x%4d x%4d   ->  %4d x%4d x%4d\n",l.outputs,  l.stride, l.w, l.h, l.c, l.out_w, l.out_h, l.out_c);
				output_w = 26;
				output_h = 32*13;

				reorg_cpu(in_ptr[i], output_w, output_h, 4, 2, out_ptr[i]);
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
				tmp_ptr = in_ptr[i];
				for(uint32_t k = 0; k<13*13*425; k++)
				{
					region_buf[k] = tmp_ptr[k];
				}
				forward_region_layer(l, region_buf);
				break;
		}
    }

	free(Memory_buf);
	free(Weight_buf);
	free(Beta_buf);
	free(region_buf);

}
///////////////////////////////////////////////////////////////////////20181229 anti-reorg ok end n4m32


//#define MEM_LEN (416*416*32+208*208*32)
//void generate_iofm_offset(float* in_ptr[32], float* out_ptr[32], float *Memory_buf, network *net)
//{
//#define ROUTE16_LEN (26*32*512)
//#define CONV27_LEN (13*16*256)
//#define CONV24_LEN (13*16*1024)
//
//	float *Memory_top = Memory_buf+512;
//	float *Memory_bottom = Memory_top + MEM_LEN;
//	int x;
//	for(x=0;x<18;x++)
//	{
//		int out_w = net->layers[x].out_w;
//		int out_w_align_256b = (out_w >> 3) << 3;
//		if(out_w & 0x7)
//			out_w_align_256b += 8;
//
//		if(x%2==0)
//		{
//			in_ptr[x] = Memory_top;
//			out_ptr[x] = Memory_bottom - net->layers[x].out_c *  net->layers[x].out_h * out_w_align_256b;
//		}
//		else
//		{
//			in_ptr[x] = out_ptr[x-1];
//			out_ptr[x] = Memory_top;
//		}
//	}
//
//	for(x=18;x<25;x++)
//	{
//		int out_w = net->layers[x].out_w;
//		int out_w_align_256b = (out_w >> 3) << 3;
//		if(out_w & 0x7)
//			out_w_align_256b += 8;
//
//		if(x%2==0)
//		{
//			in_ptr[x] = Memory_top;
//			out_ptr[x] = Memory_bottom - ROUTE16_LEN - net->layers[x].out_c *  net->layers[x].out_h * out_w_align_256b;
//		}else
//		{
//			in_ptr[x] = out_ptr[x-1];
//			out_ptr[x] = Memory_top;
//		}
//	}
//
//	in_ptr[26] = Memory_bottom - ROUTE16_LEN;
//	out_ptr[26] = Memory_top;
//
//	in_ptr[27] = Memory_top;
//	out_ptr[27] = Memory_bottom - ROUTE16_LEN - CONV24_LEN - CONV27_LEN;
//
//	in_ptr[29] = out_ptr[27];
//	out_ptr[29] = Memory_top;
//
//	in_ptr[30] = Memory_top;
//	out_ptr[30] = Memory_bottom - (net->layers[30].outputs + 3*13*425);
//
//	in_ptr[31] = out_ptr[30];
//}
//
//void reorg_cpu(float *x, int w, int h, int c, int stride, float *out)
//{
//    int i,j,k;
//    int out_c = c/(stride*stride);
//
//	for(k = 0; k < c; ++k){
//	    for(j = 0; j < h; ++j){
//		for(i = 0; i < w; ++i){
//		    int in_index  = i + w*(j + h*k);
//		    int c2 = k % out_c;
//		    int offset = k / out_c;
//		    int w2 = i*stride + offset % stride;
//		    int h2 = j*stride + offset / stride;
//		    int out_index = w2 + w*stride*(h2 + h*stride*c2);
//		    out[in_index] = x[out_index];
//		}
//	    }
//	}
//}
//
//void yolov2_hls_ps(network *net, float *input)
//{
//	int weight_offset[32] = {864, 18432, 73728, 8192, 73728,
//		294912, 32768, 294912, 1179648, 131072, 1179648, 131072,
//		1179648, 4718592, 524288, 4718592, 524288, 4718592, 9437184,
//		9437184, 32768, 11796480, 435200, 0, 0, 0, 0, 0, 0, 0, 0, 0};
//	int beta_offset[32] = {32, 64, 128, 64, 128, 256, 128, 256, 512, 256, 512, 256, 512, 1024,
//		512, 1024, 512, 1024, 1024, 1024, 64, 1024, 425, 0, 0, 0, 0, 0, 0, 0, 0, 0};
//
//	float *Weight_buf = (float *)calloc(203767168/4,sizeof(float));
//	float *Beta_buf   = (float *)calloc(43044/4,sizeof(float));
//	float *tmp_ptr_f0;
//
//	FILE *fp_w = fopen("weights_reorg.bin", "rb");
//    	if(!fp_w) file_error("weights_reorg.bin");
//
//	FILE *fp_b = fopen("bias.bin", "rb");
//    	if(!fp_b) file_error("bias.bin");
//
//	fread(Weight_buf, sizeof(float), 203767168/4, fp_w);
//	fread(Beta_buf, sizeof(float), 43044/4, fp_b);
//
//	fclose(fp_w);
//	fclose(fp_b);
//
////leave some memories for overflow, because the load_module will load extra pixels near boundary for padding
//	float *Memory_buf = (float*)calloc(MEM_LEN+512*2,sizeof(float));
//	float* in_ptr[32];
//	float* out_ptr[32];
//	generate_iofm_offset( in_ptr, out_ptr, Memory_buf, net);
//
//	memcpy(in_ptr[0], input, 416*416*3*sizeof(float));//416x416x3 input_pic
//
//	float *region_buf = (float *)calloc(13*16*425,sizeof(float));
//	float *region_buf2 = (float *)calloc(13*16*425,sizeof(float));
//
//    	int i;
//	int offset_index = 0;
//	int woffset = 0;
//	int boffset = 0;
//	int TR,TC,TM,TN;
//	int output_w,output_h;
//	int mLoops;
//
//    	for(i = 0; i < net->n; ++i)
//	{
//        	layer l = net->layers[i];
//		printf("Layer[%2d]: ",i);
//		switch(l.type)
//		{
//			case CONVOLUTIONAL:
//				printf("outputMemory:%8d;BN=%d;Activation=%d;conv  %5d %2d x%2d /%2d  %4d x%4d x%4d   ->  %4d x%4d x%4d  %5.3f BFLOPs\n",l.outputs,l.batch_normalize,l.activation, l.n, l.size, l.size, l.stride, l.w, l.h, l.c, l.out_w, l.out_h, l.out_c, (2.0 * l.n * l.size*l.size*l.c/l.groups * l.out_h*l.out_w)/1000000000.);
//
//				output_w = (l.w - l.size + 2*l.pad)/l.stride + 1;
//				output_h = (l.h - l.size + 2*l.pad)/l.stride + 1;
//
//				TR = MIN(((OnChipIB_Height-l.size)/l.stride+1),Tr);//keep Kstride>=1
//				TR = MIN(output_h,TR);
//				TC = MIN(((OnChipIB_Width-l.size)/l.stride+1),Tc);
//				TC = MIN(output_w,TC);
//				TM = MIN(l.n,Tm);
//				TN = MIN(l.c,Tn);
//				mLoops = (int)ceil(((float)l.n)/TM);
//
//				YOLO2_FPGA((ap_uint<256> *)in_ptr[i],(ap_uint<256> *)out_ptr[i],
//						(ap_uint<256> *)(Weight_buf+woffset),(ap_uint<256> *)(Beta_buf+boffset),
//					l.c,l.n,l.size,
//					l.stride,l.w,l.h,output_w, output_h, l.pad,l.activation==LEAKY?1:0,l.batch_normalize?1:0,
//					TM,TN,TR,TC, (mLoops + 1)*TM, mLoops*TM, (mLoops + 1)*TM, 0);
//
//				woffset += weight_offset[offset_index];
//				boffset += beta_offset[offset_index];
//				offset_index++;
//
//				break;
//			case MAXPOOL:
//				printf("outputMemory:%8d;max          %d x %d / %d  %4d x%4d x%4d   ->  %4d x%4d x%4d\n",l.outputs, l.size, l.size, l.stride, l.w, l.h, l.c, l.out_w, l.out_h, l.out_c);
//				//output_w = (l.w - l.size)/l.stride + 1 ;
//				//output_h = (l.h - l.size)/l.stride + 1 ;
//				output_w = l.out_h;
//				output_h = l.out_w;
//
//				TR = MIN(((OnChipIB_Height-l.size)/l.stride+1),Tr);//keep Kstride>=1
//				TC = MIN(((OnChipIB_Width-l.size)/l.stride+1),Tc);
//				TR = MIN(output_h,TR);
//				TC = MIN(output_w,TC);
//				TM = MIN(Tm,Tn);
//				TM = MIN(l.c,TM);
//				mLoops = (int)ceil(((float)l.c)/TM);
//
//				YOLO2_FPGA((ap_uint<256> *)in_ptr[i],(ap_uint<256> *)out_ptr[i],NULL,NULL,l.c,l.c,
//					l.size,l.stride,l.w,l.h, output_w, output_h, l.pad,0,0,TM,0,TR,TC, (mLoops + 2)*TM, mLoops*TM, (mLoops + 1)*TM, 1);
//
//				break;
//			case REORG:
//				printf("outputMemory:%8d;reorg              /%2d  %4d x%4d x%4d   ->  %4d x%4d x%4d\n",l.outputs,  l.stride, l.w, l.h, l.c, l.out_w, l.out_h, l.out_c);
//				output_w = 26;
//				output_h = 32*13;
//
//				TR = MIN(((OnChipIB_Height-l.stride)/l.stride+1),Tr);//keep Kstride>=1
//				TR = MIN(output_h,TR);
//				TC = MIN(((OnChipIB_Width-l.stride)/l.stride+1),Tc);
//				TC = MIN(output_w,TC);
//				//TM = 4;
//				//mLoops = 1;
//				//here Tm and Tn must >=4
//				TM = MIN(Tm,Tn);
//				TM = MIN(4,TM);
//				mLoops = (int)ceil(((float)4)/TM);
//
//			//	YOLO2_FPGA(in_ptr[i],out_ptr[i],NULL,NULL,1,4,
//			//				  l.stride,l.stride,52,32*26, output_w, output_h, 0,0,0,TM,0,TR,TC, (mLoops + 2)*TM, mLoops*TM, (mLoops + 1)*TM, 2);
//				tmp_ptr_f0 = in_ptr[i];
//				for(int k = 0; k<26*64; k++)
//					memcpy((float *)(region_buf + k*26), (float *)(tmp_ptr_f0 + k*32), 26*sizeof(float));
//				reorg_cpu(region_buf, output_w, output_h, 4, 2, region_buf2);
//				tmp_ptr_f0 = region_buf;
//				memset(region_buf, 0,  13*16*256*sizeof(float));
//				for(int k = 0; k<13*256; k++)
//					memcpy((float *)(tmp_ptr_f0 + k*16), (float *)(region_buf2 + k*13), 13*sizeof(float));
//				memcpy(out_ptr[i], tmp_ptr_f0, 13*16*256*sizeof(float));
//
//				break;
//			case ROUTE:
//				printf("outputMemory:%8d;route ",l.outputs);
//				int j;
//				for(j = 0; j < l.n; ++j){
//					printf(" %d", l.input_layers[j]);
//				}
//				printf("\n");
//				break;
//			case REGION:
//				printf("outputMemory:%8d;Detection\n",l.outputs);
//				tmp_ptr_f0 = in_ptr[i];
//				for(int k = 0; k<13*425; k++)
//					for(int j = 0; j < 16; j++)
//					{
//						if(j < 13)
//							region_buf[k*13 + j] = tmp_ptr_f0[k*16 + j];
//					}
//				forward_region_layer(l, region_buf);
//				break;
//		}
//    }
//
//	free(Memory_buf);
//	free(Weight_buf);
//	free(Beta_buf);
//	free(region_buf);
//	free(region_buf2);
//
//}
/////////////////////////////////////////////////////////////////////////20181229 anti-reorg ok end n4m32
