

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))

#DEFINE_HEADER#

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
	printf("%s's data size is %ld\n", filename, buf.st_size);  
    return (unsigned long)buf.st_size;  
}

const uint32_t IB_W = OnChipIB_Width;
const uint32_t IB_H = OnChipIB_Height;
const uint32_t IB_HxW = IB_H*IB_W;
// const uint32_t TnxIB_H = Tn*IB_H;	
// const uint32_t TnxIB_HxIB_W = Tn*IB_H*IB_W;
const uint32_t TrxTc = Tr*Tc;

void input_load(float *input,float input_buffer[Tn][OnChipIB_Height][OnChipIB_Width],int r,int c,int n,int Kstride,int Padding,int TRow,int TCol,int Input_w,int Input_h,int TN_MIN,int IHxIW,int LayerType)
{
	static float local_buf[OnChipIB_Width];

	const int Coffset = c*Kstride - Padding;
	const int Roffset = r*Kstride - Padding;
	const int CurrentOffset = n*IHxIW + Roffset*Input_w + Coffset;

	float pad_value = 0.0f;
	if(LayerType==1)
		pad_value = -1024*1024;

	for(uint16_t t1 = 0;t1 < Tn; t1++)
	for(uint16_t t2 = 0;t2 < TRow; t2++)
	{
		int ifm_offset = CurrentOffset + t1*IHxIW + t2*Input_w;
		memcpy( local_buf,(float *)(input + ifm_offset), TCol*sizeof(float));

		bool TN_Enable = t1 < TN_MIN;
		int yoffset = Roffset + t2;
		bool YEnable = (yoffset >= 0)&&(yoffset < Input_h);
		bool PEnable = YEnable&&TN_Enable;

		for(uint16_t t3 = 0;t3 < TCol; t3++)
		{
			int xoffset = Coffset + t3;
			bool XEnable = (xoffset >= 0)&&(xoffset < Input_w);
			if(XEnable&&PEnable)
			{
				input_buffer[t1][t2][t3] = local_buf[t3];
			}
			else
				input_buffer[t1][t2][t3] = pad_value;
		}
	}
}

void weight_load(float *Weight,float weight_buffer[Tm][Tn][K][K],bool weight_load_enable, uint16_t m, uint16_t n,int IFM_numxKxK, uint16_t KxK, uint16_t Ksize, uint16_t TM_MIN, uint16_t TN_MIN)
{
	static float weight_memcpy_buffer[Tm*Tn*K*K];

	if(!weight_load_enable)
		return;

	const int Woffset = m*IFM_numxKxK + n*KxK;
	int w_cp_idx = 0;
	for(uint16_t t1 = 0;t1 < TM_MIN; t1++)
		for(uint16_t t2 = 0;t2 < TN_MIN; t2++)
		{
			memcpy((float *)(weight_memcpy_buffer + w_cp_idx),(float *)(Weight + Woffset + t1*IFM_numxKxK + t2*KxK),KxK*sizeof(float));
			w_cp_idx += KxK;
		}

	w_cp_idx = 0;
	for(uint16_t t1 = 0;t1 < Tm; t1++)
	for(uint16_t t2 = 0;t2 < Tn; t2++)
	for(uint16_t t3 = 0;t3 <Ksize; t3++)
	for(uint16_t t4 = 0;t4 <Ksize; t4++)
	{
		bool Enable = (t1 < TM_MIN)&&(t2 < TN_MIN);
		if(Enable)
		{
			weight_buffer[t1][t2][t3][t4] =  weight_memcpy_buffer[w_cp_idx];
			w_cp_idx++;
		}
		else
			weight_buffer[t1][t2][t3][t4] = 0;
	}
}

void copy_input_weight(float *input,float *Weight,int IFM_num,int Input_w,int Input_h, int Ksize,int Kstride,int r,int c,int m,int n,
		int TM_MIN,int TN,int TRow,int TCol,int Padding,float input_buffer[Tn][OnChipIB_Height][OnChipIB_Width],float weight_buffer[Tm][Tn][K][K],
		bool weight_load_enable,const int IHxIW,const int KxK,const int IFM_numxKxK,const int LayerType)
{
	uint16_t TN_MIN = MIN(TN, IFM_num-n);

	input_load(input, input_buffer, r, c, n, Kstride, Padding, TRow, TCol, Input_w, Input_h, TN_MIN, IHxIW, LayerType);
	weight_load(Weight,weight_buffer,weight_load_enable,m,n,IFM_numxKxK,KxK,Ksize,TM_MIN,TN_MIN);

}

void compute(float input_buffer[Tn][OnChipIB_Height][OnChipIB_Width],float output_buffer[Tm][Tr][Tc],
		float weight_buffer[Tm][Tn][K][K],float beta_buffer[MAX_BETA_LENGTH], uint16_t n,
		uint16_t Ksize, uint16_t Kstride, uint16_t m,
		uint16_t TM_MIN, uint16_t TR_MIN, uint16_t TC_MIN)
{
	float partial_mul[Tm][Tn];
	float partial_add[Tm];

	for(uint16_t i =0;i < Ksize; i++)
	for(uint16_t j = 0;j < Ksize; j++)
	for(uint16_t tr = 0;tr < TR_MIN;tr++)
	for(uint16_t tc = 0;tc < TC_MIN;tc++)
	{
		for(uint16_t tm = 0;tm < Tm;tm++)
		{
			if(i==0&&j==0&&n==0)
				partial_add[tm] = beta_buffer[tm+m];
			else
				partial_add[tm] = output_buffer[tm][tr][tc];

			for(uint16_t tn = 0;tn <Tn;tn++)
			{
				partial_mul[tm][tn] = weight_buffer[tm][tn][i][j]*input_buffer[tn][Kstride*tr+i][Kstride*tc+j];
			}

			float partial_sum = 0;
			for(uint16_t tn = 0;tn <Tn;tn++)
			{
					partial_sum += partial_mul[tm][tn];
			}
			output_buffer[tm][tr][tc] = partial_add[tm] + partial_sum;
		}
	}
}

void write_back_output_reorg(float output_buffer[Tm][Tr][Tc], float *Output, uint16_t r, uint16_t c, uint16_t m, uint16_t Output_w, uint16_t Output_h,
		uint16_t TM_MIN,uint16_t TR_MIN, uint16_t TC_MIN, uint32_t OHxOW, bool IsNL)
{
	assert((TM_MIN >0)&&(TM_MIN <=Tm));
	assert((TR_MIN >0)&&(TR_MIN <=Tr));
	assert((TC_MIN >0)&&(TC_MIN <=Tc));

	uint32_t offset = m*OHxOW + r*Output_w + c;
	static float local_buf[Tc];

	for(uint16_t tm = 0; tm < TM_MIN; tm++)
	for(uint16_t tr = 0; tr < TR_MIN; tr++)
	{
		float tmp_out;
		for(uint16_t tc = 0;tc < TC_MIN;tc++)
		{
			float tmp = output_buffer[tm][tr][tc];
			if((tmp < 0.0f)&&IsNL)
				tmp_out = tmp*0.1f;
			else
				tmp_out = tmp;			
			local_buf[tc] = tmp_out;
		}

		uint32_t ofm_offset = tm*OHxOW + tr*Output_w + offset;
		memcpy((float *)(Output + ofm_offset), local_buf, TC_MIN*sizeof(float));	
	}
}

void pool_yolo2(float Input[Tn][OnChipIB_Height][OnChipIB_Width],float Output[Tm][Tr][Tc],
		  uint16_t Ksize, uint16_t Kstride, uint16_t TM_MIN, uint16_t TR_MIN, uint16_t TC_MIN)
{
	float tmp[Tn];

	for(uint16_t tr = 0;tr < TR_MIN;tr++)
	for(uint16_t tc = 0;tc < TC_MIN;tc++)
	for(uint16_t i =0;i < Ksize; i++)
	for(uint16_t j = 0;j < Ksize; j++)
	{
		for(uint16_t of = 0; of < Tn; of++)
		{
			if(i==0&&j==0)
				tmp[of] = -1024*1024;

			if(Input[of][tr*Kstride+i][tc*Kstride+j] > tmp[of])
				tmp[of] = Input[of][tr*Kstride+i][tc*Kstride+j];

			if(i==1&&j==1)
				Output[of][tr][tc] = tmp[of];
		}
	}
}

void YOLO2_FPGA(float *Input,float *Output,float *Weight,float *Beta, int IFM_num, int OFM_num,
							   int Ksize, int Kstride,
							   int Input_w, int Input_h, int Output_w, int Output_h, int Padding, bool IsNL, bool IsBN,
							   int TM, int TN, int TR, int TC,
							   int OFM_num_bound, int mLoopsxTM, int mLoops_a1xTM, int LayerType)
{
	assert((OFM_num > 0)&&(OFM_num <= 2048));
	assert((IFM_num > 0)&&(IFM_num <= 2048));
	assert((Kstride > 0)&&(Kstride <= S));
	assert((Ksize > 0)&&(Ksize <= K));
	assert((Input_w > 0)&&(Input_w <= 1024));
	assert((Input_h > 0)&&(Input_h <= 1024));
	assert((Output_w > 0)&&(Output_w <= 1024));
	assert((Output_h > 0)&&(Output_h <= 1024));
	assert((Padding >= 0)&&(Padding <= 4));//maybe
	assert((TM > 0)&&(TM <= Tm));
	assert((TN >= 0)&&(TN <= Tn));
	// assert((TR > 0)&&(TR <= Tr));
	// assert((TC > 0)&&(TC <= Tc));
	assert(TR > 0);
	assert(TC > 0);	
	assert((TR*TC)<=(Tr*Tc));	

	const int OHxOW = Output_h*Output_w;
	const int TRow = (TR-1)*Kstride+Ksize;
	const int TCol = (TC-1)*Kstride+Ksize;
	const int IHxIW   = Input_h*Input_w;
	const int KxK = Ksize*Ksize;
	const int IFM_numxKxK = IFM_num*KxK;

	static float input_buffer[Tn][OnChipIB_Height][OnChipIB_Width];
	static float weight_buffer[Tm][Tn][K][K];
	static float output_buffer[Tm][Tr][Tc];
	static float beta_buffer[MAX_BETA_LENGTH];

	if(LayerType==0)
		memcpy(beta_buffer,Beta, OFM_num*sizeof(float));

	for(uint16_t r = 0; r < Output_h; r += TR)
	{
		uint16_t TR_MIN = MIN(TR,Output_h-r);
		for(uint16_t c = 0; c < Output_w; c += TC)
		{
			uint16_t TC_MIN = MIN(TC,Output_w-c);
			for(uint16_t m = 0; m < OFM_num; m += TM)
			{
				uint16_t TM_MIN = MIN(TM, OFM_num-m);
				if(LayerType==0)
				{
					for(uint16_t n = 0;n < IFM_num; n += TN)
					{
						copy_input_weight(Input,Weight,IFM_num,Input_w,Input_h,Ksize,Kstride,r,c,m, n,
							TM_MIN,TN,TRow,TCol,Padding,input_buffer,weight_buffer,1,IHxIW,KxK,IFM_numxKxK,LayerType);
						compute(input_buffer,output_buffer,weight_buffer,beta_buffer, n, Ksize,Kstride,m,TM_MIN,TR_MIN,TC_MIN);
					}
				}
				else if(LayerType==1)
				{
					copy_input_weight(Input,Weight,IFM_num,Input_w,Input_h,Ksize,Kstride,r,c,m,m,
							TM_MIN,TM,TRow,TCol,0,input_buffer,weight_buffer,0,IHxIW,KxK,IFM_numxKxK,LayerType);
					pool_yolo2(input_buffer,output_buffer,Ksize,Kstride,TM_MIN,TR_MIN,TC_MIN);
				}

				write_back_output_reorg(output_buffer,Output, r, c, m,Output_w,Output_h, TM_MIN, TR_MIN, TC_MIN, OHxOW, IsNL);
			}
		}
	}
}

int Weight_reorgnaization_anti(float *Weight,float *Weight_reorg,float* Alpha,int IFM_NUM,int OFM_NUM,int Ksize,int TM,int TN,const bool IsBN)
{
	const int KxK = Ksize*Ksize;
	const int IFM_NUMxKxK = IFM_NUM*KxK;

	int m,n;
	int tm,tn,tk;

	float weight_buffer[Tm*Tn*K*K];
	float weight_buffer2[Tm*Tn*K*K];

	int TM_MIN,TN_MIN;
	int offset = 0;

	for( m = 0; m < OFM_NUM; m += TM)
	{
		TM_MIN = MIN(TM,OFM_NUM - m);
		for(n = 0;n < IFM_NUM; n += TN)
		{
			TN_MIN = MIN(TN,IFM_NUM - n);
			int Woffset = m*IFM_NUMxKxK + n*KxK;
			for(tm = 0;tm < TM_MIN; tm++)
			{
				memcpy((float *)(weight_buffer + tm*TN_MIN*KxK),
					(float *)(Weight + tm*IFM_NUMxKxK + Woffset),TN_MIN*KxK*sizeof(float));
			}

			// int TN_MINxTM_MIN = TN_MIN*TM_MIN;		
			// for(tk = 0;tk < KxK; tk++)
			// 	for(tm = 0;tm < TM_MIN; tm++)
			// 		for(tn = 0;tn < TN_MIN;tn++)
			// 		{
			// 			weight_buffer2[tk*TN_MINxTM_MIN + tm*TN_MIN + tn] = weight_buffer[tm*TN_MIN*KxK + tn*KxK + tk];
			// 		}

			// memcpy((float *)(Weight_reorg+offset),weight_buffer2,TM_MIN*TN_MIN*KxK*sizeof(float));
			memcpy((float *)(Weight_reorg+offset),weight_buffer,TM_MIN*TN_MIN*KxK*sizeof(float));
			offset += TM_MIN*TN_MIN*KxK;
		}							
	}

	return 0;
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

void yolov2_hls_ps(network *net, float *input)
{
	int weight_offset[32] = {864, 18432, 73728, 8192, 73728,
		294912, 32768, 294912, 1179648, 131072, 1179648, 131072,
		1179648, 4718592, 524288, 4718592, 524288, 4718592, 9437184,
		9437184, 32768, 11796480, 435200, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	int beta_offset[32] = {32, 64, 128, 64, 128, 256, 128, 256, 512, 256, 512, 256, 512, 1024,
		512, 1024, 512, 1024, 1024, 1024, 64, 1024, 425, 0, 0, 0, 0, 0, 0, 0, 0, 0};

	uint32_t data_num = 0;

	data_num = get_file_size("weights.bin")/4;
	float *Weight_buf = (float *)calloc(data_num,sizeof(float));
	FILE *fp_w = fopen("weights.bin", "rb");
    if(!fp_w) file_error("weights.bin");
	fread(Weight_buf, sizeof(float), data_num, fp_w);
	fclose(fp_w);

	float *Weight_reorg_buf = (float *)calloc(data_num,sizeof(float));
	FILE *fp_w_reorg = fopen("weights_reorg.bin", "wb");
    if(!fp_w_reorg) file_error("weights_reorg.bin");

	data_num = get_file_size("bias.bin")/4;
	float *Beta_buf   = (float *)calloc(data_num, sizeof(float));
	FILE *fp_b = fopen("bias.bin", "rb");
    if(!fp_b) file_error("bias.bin");
	fread(Beta_buf, sizeof(float), data_num, fp_b);
	fclose(fp_b);

	float *Memory_buf = (float*)calloc(MEM_LEN+512*2,sizeof(float));
	float* in_ptr[32];
	float* out_ptr[32];
	generate_iofm_offset( in_ptr, out_ptr, Memory_buf, net);

	memcpy(in_ptr[0], input, 416*416*3*sizeof(float));//416x416x3 input_pic

    int i;
	int offset_index = 0;
	int woffset = 0;
	int boffset = 0;
	int TR,TC,TM,TN;
	int TRow, TCol;
	int output_w,output_h;
	int mLoops;

    for(i = 0; i < net->n; ++i)
	{
        	layer l = net->layers[i];
		printf("Layer[%2d]: ",i);
		switch(l.type)
		{
			case CONVOLUTIONAL:
				printf("outputMemory:%8d;BN=%d;Activation=%d;conv  %5d %2d x%2d /%2d  %4d x%4d x%4d   ->  %4d x%4d x%4d  %5.3f BFLOPs\n", 
				l.outputs,l.batch_normalize,l.activation, l.n, l.size, l.size, l.stride, l.w, l.h, l.c, l.out_w, l.out_h, l.out_c, (2.0 * l.n * l.size*l.size*l.c/l.groups * l.out_h*l.out_w)/1000000000.);

				output_w = (l.w - l.size + 2*l.pad)/l.stride + 1;
				output_h = (l.h - l.size + 2*l.pad)/l.stride + 1;

				// assert((IB_HxW/l.size)>=l.size);
				// TC = MIN(((IB_HxW/l.size-l.size)/l.stride+1),output_w);
				// TC = MIN(TrxTc, TC);
				// TCol = (TC-1)*l.stride + l.size;
				// TR = MIN(((IB_HxW/TCol-l.size)/l.stride+1),output_h);//keep Kernel_stride>=1
				// TR = MIN(TR, TrxTc/TC);
				// TRow = (TR-1)*l.stride + l.size;

				// assert(((TR*TC)>0)&&((TR*TC)<=TrxTc));
				// assert(((TRow*TCol)>0)&&((TRow*TCol)<=IB_HxW));
				// printf("TR=%d, TC=%d, TRow=%d, TCol=%d\n", TR, TC, TRow, TCol);				

				TR = MIN(((OnChipIB_Height-l.size)/l.stride+1),Tr);//keep Kstride>=1
				TR = MIN(output_h,TR);
				TC = MIN(((OnChipIB_Width-l.size)/l.stride+1),Tc);
				TC = MIN(output_w,TC);

				TM = MIN(l.n,Tm);
				TN = MIN(l.c,Tn);
				mLoops = (int)ceil(((float)l.n)/TM);

				YOLO2_FPGA(in_ptr[i],out_ptr[i],Weight_buf+woffset,Beta_buf+boffset,
					l.c,l.n,l.size,
					l.stride,l.w,l.h,output_w, output_h, l.pad,l.activation==LEAKY?1:0,l.batch_normalize?1:0,
					TM,TN,TR,TC, (mLoops + 1)*TM, mLoops*TM, (mLoops + 1)*TM, 0);

				Weight_reorgnaization_anti(Weight_buf + woffset,Weight_reorg_buf + woffset,NULL,l.c,l.n,l.size,TM,TN,0);

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

				TR = MIN(((OnChipIB_Height-l.size)/l.stride+1),Tr);//keep Kstride>=1
				TC = MIN(((OnChipIB_Width-l.size)/l.stride+1),Tc);
				TR = MIN(output_h,TR);
				TC = MIN(output_w,TC);

				// assert((IB_HxW/l.size)>=l.size);
				// TC = MIN(((IB_HxW/l.size-l.size)/l.stride+1),output_w);
				// TC = MIN(TrxTc, TC);
				// TCol = (TC-1)*l.stride + l.size;
				// TR = MIN(((IB_HxW/TCol-l.size)/l.stride+1),output_h);//keep Kernel_stride>=1
				// TR = MIN(TR, TrxTc/TC);
				// TRow = (TR-1)*l.stride + l.size;

				// assert(((TR*TC)>0)&&((TR*TC)<=TrxTc));
				// assert(((TRow*TCol)>0)&&((TRow*TCol)<=IB_HxW));
				// printf("TR=%d, TC=%d, TRow=%d, TCol=%d\n", TR, TC, TRow, TCol);									

				TM = MIN(Tm,Tn);
				TM = MIN(l.c,TM);
				mLoops = (int)ceil(((float)l.c)/TM);

				YOLO2_FPGA(in_ptr[i],out_ptr[i],NULL,NULL,l.c,l.c,
					l.size,l.stride,l.w,l.h, output_w, output_h, l.pad,0,0,TM,0,TR,TC, (mLoops + 2)*TM, mLoops*TM, (mLoops + 1)*TM, 1);

				break;
			case REORG:
				printf("outputMemory:%8d;reorg              /%2d  %4d x%4d x%4d   ->  %4d x%4d x%4d\n",l.outputs,  l.stride, l.w, l.h, l.c, l.out_w, l.out_h, l.out_c);			
				output_w = 26;
				output_h = 32*13;
				reorg_cpu(in_ptr[i], output_w, output_h, 4, 2, out_ptr[i]);
				break;
			case ROUTE:
				printf("outputMemory:%8d;route ",l.outputs);
				int j;
				for(j = 0; j < l.n; ++j){
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

	data_num = get_file_size("weights.bin")/4;
	fwrite(Weight_reorg_buf, sizeof(float), data_num, fp_w_reorg);
	fclose(fp_w_reorg);
	free(Weight_reorg_buf);

	free(Memory_buf);
	free(Weight_buf);
	free(Beta_buf);

}
///////////////////////////////////////////////////////////////////////20181229 anti-reorg ok end n4m32
