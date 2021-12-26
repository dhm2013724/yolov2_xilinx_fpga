

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))

#DEFINE_HEADER#

//#define REORG_GEN
//#define REORG_TEST

#include <assert.h>
#include <math.h>

#define PRAGMA_SUB(x) _Pragma (#x)
#define DO_PRAGMA(x) PRAGMA_SUB(x)

static int ReadByte = 0;
static int WriteByte = 0;

typedef float IO_Dtype;

void ifm_mmcpy_row(IO_Dtype *input, IO_Dtype input_memcpy_buffer[OnChipIB_Width], int CurrentOffset, int IHxIW, int Input_w, int TCol, 
	uint8_t t1, uint8_t t2, uint8_t *t1_n, uint8_t *t2_n,bool enable)
{
	if(!enable)
		return;

	int ifm_offset = CurrentOffset + t1*IHxIW + t2*Input_w;
	memcpy( input_memcpy_buffer,(IO_Dtype *)(input + ifm_offset), TCol*sizeof(IO_Dtype));
	
	ReadByte += TCol*sizeof(IO_Dtype);

	*t1_n = t1;
	*t2_n = t2;
}

void ifm_copy_lbuf2ibuf(float input_buffer[Tn][OnChipIB_Height][OnChipIB_Width], IO_Dtype local_buf[OnChipIB_Width], int TCol, int Input_w, int Input_h, int TN_MIN, float pad_value,
	int Coffset, int Roffset, uint8_t t1, uint8_t t2, bool enable)
{
	if(!enable)
		return;

	bool TN_Enable = t1 < TN_MIN;
	int yoffset = Roffset + t2;
	bool YEnable = (yoffset >= 0)&&(yoffset < Input_h);
	bool PEnable = YEnable&&TN_Enable;

	uint8_t t3;
	for(t3 = 0;t3 < TCol; t3++)
	{
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=TCol_max)
#pragma HLS PIPELINE II=1
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

void input_load(IO_Dtype *input,float input_buffer[Tn][OnChipIB_Height][OnChipIB_Width],int r,int c,int n,int Kstride,int Padding,int TRow,int TCol,int Input_w,int Input_h,int TN_MIN,int IHxIW,int LayerType)
{
	uint8_t t1,t2;
	static IO_Dtype local_buf0[OnChipIB_Width];
	static IO_Dtype local_buf1[OnChipIB_Width];

	const int Coffset = c*Kstride - Padding;
	const int Roffset = r*Kstride - Padding;
	const int CurrentOffset = n*IHxIW + Roffset*Input_w + Coffset;

	uint8_t t1_n0, t1_n1, t2_n0, t2_n1;
	bool pp = true;
	

	float pad_value = 0.0f;
	if(LayerType==1)
		pad_value = -1024*1024;

//1. normal mode:just read one row of tile each time
//we know that the TRow is always smaller than TRow_max
//2. continous mode: read several rows of tile
//be careful, we should confirm that CColNum*TCol <= TCol_max, and
//this rows must be continous <=> which also limits that TCol >= Input_w(consider Padding)
//
	
	int TnxTRow = Tn*TRow;
	int t = 0;
	t1 = 0; t2 = 0;
	for(t = 0;t < TnxTRow+1; t++)
	{
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=TRow_max)
		if(pp)
		{
			ifm_mmcpy_row(input, local_buf0, CurrentOffset, IHxIW, Input_w, TCol, t1, t2, &t1_n0, &t2_n0, t!=TnxTRow);
			ifm_copy_lbuf2ibuf( input_buffer, local_buf1, TCol, Input_w, Input_h, TN_MIN, pad_value, Coffset, Roffset, t1_n1, t2_n1, t!=0);
			pp = false;
		}else
		{
			ifm_mmcpy_row(input, local_buf1, CurrentOffset, IHxIW, Input_w, TCol, t1, t2, &t1_n1, &t2_n1, t!=TnxTRow);
			ifm_copy_lbuf2ibuf( input_buffer, local_buf0, TCol, Input_w, Input_h, TN_MIN, pad_value, Coffset, Roffset, t1_n0, t2_n0, t!=0);
			pp = true;
		}
		
		t2++;
		if(t2==TRow)
		{
			t2 = 0;
			t1++;
		}
	}


/*
	for(t1 = 0;t1 < Tn; t1++)
		for(t2 = 0;t2 < TRow + 1; t2++)
		{
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=TRow_max)
			if(pp)
			{
				ifm_mmcpy_row(input, local_buf0, CurrentOffset, IHxIW, Input_w, TCol, t1, t2, &t1_n0, &t2_n0, t2!=TRow);
				ifm_copy_lbuf2ibuf( input_buffer, local_buf1, TCol, Input_w, Input_h, TN_MIN, pad_value, Coffset, Roffset, t1_n1, t2_n1, t2!=0);
				pp = false;
			}else
			{
				ifm_mmcpy_row(input, local_buf1, CurrentOffset, IHxIW, Input_w, TCol, t1, t2, &t1_n1, &t2_n1, t2!=TRow);
				ifm_copy_lbuf2ibuf( input_buffer, local_buf0, TCol, Input_w, Input_h, TN_MIN, pad_value, Coffset, Roffset, t1_n0, t2_n0, t2!=0);
				pp = true;
			}
		}
*/

}

/*
void ifm_mmcpy_row(IO_Dtype *input, IO_Dtype input_memcpy_buffer[OnChipIB_Width], int CurrentOffset, int IHxIW, int Input_w, int TCol, 
	uint8_t t1, uint8_t t2, uint8_t *t1_n, uint8_t *t2_n,bool enable)
{
	if(!enable)
		return;

	int ifm_offset = CurrentOffset + t1*IHxIW + t2*Input_w;
	memcpy( input_memcpy_buffer,(IO_Dtype *)(input + ifm_offset), TCol*sizeof(IO_Dtype));

	*t1_n = t1;
	*t2_n = t2;
}

void ifm_copy_lbuf2ibuf(float input_buffer[Tn][OnChipIB_Height][OnChipIB_Width], IO_Dtype local_buf[OnChipIB_Width], int TCol, int Input_w, int Input_h, int TN_MIN, float pad_value,
	int Coffset, int Roffset, uint8_t t1, uint8_t t2, bool enable)
{
	if(!enable)
		return;

	bool TN_Enable = t1 < TN_MIN;
	int yoffset = Roffset + t2;
	bool YEnable = (yoffset >= 0)&&(yoffset < Input_h);
	bool PEnable = YEnable&&TN_Enable;

	uint8_t t3;
	for(t3 = 0;t3 < TCol; t3++)
	{
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=TCol_max)
#pragma HLS PIPELINE II=1
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

void input_load(IO_Dtype *input,float input_buffer[Tn][OnChipIB_Height][OnChipIB_Width],int r,int c,int n,int Kstride,int Padding,int TRow,int TCol,int Input_w,int Input_h,int TN_MIN,int IHxIW,int LayerType)
{
	uint8_t t1,t2;
	static IO_Dtype local_buf0[OnChipIB_Width];
	static IO_Dtype local_buf1[OnChipIB_Width];

	const int Coffset = c*Kstride - Padding;
	const int Roffset = r*Kstride - Padding;
	const int CurrentOffset = n*IHxIW + Roffset*Input_w + Coffset;

	uint8_t t1_n0, t1_n1, t2_n0, t2_n1;
	bool pp = true;
	

	float pad_value = 0.0f;
	if(LayerType==1)
		pad_value = -1024*1024;

	for(t1 = 0;t1 < Tn; t1++)
		for(t2 = 0;t2 < TRow + 1; t2++)
		{
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=TRow_max)
			if(pp)
			{
				ifm_mmcpy_row(input, local_buf0, CurrentOffset, IHxIW, Input_w, TCol, t1, t2, &t1_n0, &t2_n0, t2!=TRow);
				ifm_copy_lbuf2ibuf( input_buffer, local_buf1, TCol, Input_w, Input_h, TN_MIN, pad_value, Coffset, Roffset, t1_n1, t2_n1, t2!=0);
				pp = false;
			}else
			{
				ifm_mmcpy_row(input, local_buf1, CurrentOffset, IHxIW, Input_w, TCol, t1, t2, &t1_n1, &t2_n1, t2!=TRow);
				ifm_copy_lbuf2ibuf( input_buffer, local_buf0, TCol, Input_w, Input_h, TN_MIN, pad_value, Coffset, Roffset, t1_n0, t2_n0, t2!=0);
				pp = true;
			}
		}

}
*/

/*
void input_load(float *input,float input_buffer[Tn][OnChipIB_Height][OnChipIB_Width],int r,int c,int n,int Kstride,int Padding,int TRow,int TCol,int Input_w,int Input_h,int TN_MIN,int IHxIW,int LayerType)
{
	uint8_t t1,t2,t3,t4;
	int xoffset;
	int yoffset;

	static float input_memcpy_buffer[Tn*OnChipIB_Height*OnChipIB_Width];

	const int Coffset = c*Kstride - Padding;
	const int Roffset = r*Kstride - Padding;
	const int CurrentOffset = n*IHxIW + Roffset*Input_w + Coffset;

	float pad_value = 0;
	if(LayerType==1)
		pad_value = -1024*1024;

	int input_mmcpy_offset = 0;
	for(t1 = 0;t1 < TN_MIN; t1++)
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tn)
		for(t2 = 0;t2 < TRow; t2++)
		{
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=TRow_max)
			memcpy((float *)(input_memcpy_buffer + input_mmcpy_offset),(float *)(input + CurrentOffset + t1*IHxIW + t2*Input_w),TCol*sizeof(float));
			input_mmcpy_offset += TCol;
		}

	input_mmcpy_offset = 0;
	for(t1 = 0;t1 < Tn; t1++)
		for(t2 = 0;t2 < TRow; t2++)
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=TRow_max)
			for(t3 = 0;t3 < TCol; t3++)
			{
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=TCol_max)
#pragma HLS PIPELINE II=1
				xoffset = Coffset + t3;
				yoffset = Roffset + t2;
				bool XEnable    = (xoffset >= 0)&&(xoffset < Input_w);
				bool YEnable    = (yoffset >= 0)&&(yoffset < Input_h);
				bool PaddingEnable = XEnable&&YEnable;
				if(PaddingEnable&&(t1 < TN_MIN))
					input_buffer[t1][t2][t3] = input_memcpy_buffer[input_mmcpy_offset];
				else
					input_buffer[t1][t2][t3] = pad_value;
				input_mmcpy_offset++;
			}
}
*/

void weight_load(float *Weight,float weight_buffer[Tm][Tn][K][K],bool weight_load_enable,int m,int n,int IFM_numxKxK,int KxK,int Ksize,int TM_MIN,int TN_MIN)
{
	uint8_t t1,t2,t3,t4;
	static float weight_memcpy_buffer[Tm*Tn*K*K];

	if(!weight_load_enable)
		return;

	const int Woffset = m*IFM_numxKxK + n*KxK;

	int weight_memcpy_offset = 0;
	for(t1 = 0;t1 < TM_MIN; t1++)
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tm)
		for(t2 = 0;t2 < TN_MIN; t2++)
		{
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tn)
			memcpy((float *)(weight_memcpy_buffer + weight_memcpy_offset),(float *)(Weight + Woffset + t1*IFM_numxKxK + t2*KxK),KxK*sizeof(float));
			weight_memcpy_offset += KxK;
		}

	weight_memcpy_offset = 0;
	for(t1 = 0;t1 < Tm; t1++)
		for(t2 = 0;t2 < Tn; t2++)
			for(t3 = 0;t3 <Ksize; t3++)
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=K)
				for(t4 = 0;t4 <Ksize; t4++)
				{
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=K)
#pragma HLS PIPELINE II=1
					bool Enable = (t1 < TM_MIN)&&(t2 < TN_MIN);
					if(Enable)
					{
						weight_buffer[t1][t2][t3][t4] =  weight_memcpy_buffer[weight_memcpy_offset];
						weight_memcpy_offset++;
					}
					else
						weight_buffer[t1][t2][t3][t4] = 0;
				}
}

void weight_load_reorg(float *Weight,float weight_buffer[Tm][Tn][K][K],bool weight_load_enable,int m,int n,int IFM_numxKxK,int KxK,int Ksize,int TM_MIN,int TN_MIN)
{
	uint8_t t1,t2,t3,t4;
	static float weight_memcpy_buffer[Tm*Tn*K*K];
	static int Woffset;

	assert((TM_MIN > 0)&&(TM_MIN <= Tm));
	assert((TN_MIN > 0)&&(TN_MIN <= Tn));
	assert((KxK > 0)&&(KxK <= K*K));

	if(!weight_load_enable)
		return;

	if(m==0&&n==0)
		Woffset = 0;

	uint16_t mm_offset = TM_MIN*TN_MIN*KxK;
	memcpy(weight_memcpy_buffer,(float *)(Weight + Woffset), mm_offset*sizeof(float));
	Woffset += mm_offset;

	ReadByte += mm_offset*sizeof(float);

	int weight_memcpy_offset = 0;
	float input_value;
	for(t3 = 0;t3 <Ksize; t3++)
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=K)
		for(t4 = 0;t4 <Ksize; t4++)
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=K)
			for(t1 = 0;t1 < Tm; t1++)
				for(t2 = 0;t2 < Tn; t2++)
				{
#pragma HLS PIPELINE II=1
					bool Enable = (t1 < TM_MIN)&&(t2 < TN_MIN);
					if(Enable)
					{
						weight_buffer[t1][t2][t3][t4] =  weight_memcpy_buffer[weight_memcpy_offset];
						weight_memcpy_offset++;
					}
					else
						weight_buffer[t1][t2][t3][t4] = 0;
				}
}


void copy_input_weight(float *input,float *Weight,int IFM_num,int Input_w,int Input_h, int Ksize,int Kstride,int r,int c,int m,int n,
		int TM_MIN,int TN,int TRow,int TCol,int Padding,float input_buffer[Tn][OnChipIB_Height][OnChipIB_Width],float weight_buffer[Tm][Tn][K][K],int n_next[1],
		bool enable,bool weight_load_enable,bool initialize,const int IHxIW,const int KxK,const int IFM_numxKxK,const int LayerType)
{
	if(!enable)
		return ;

	const int TN_MIN = MIN(TN, IFM_num-n);
	n_next[0] = n;

	input_load(input, input_buffer, r, c, n, Kstride, Padding, TRow, TCol, Input_w, Input_h, TN_MIN, IHxIW, LayerType);
#ifdef REORG_TEST
	weight_load_reorg(Weight,weight_buffer,weight_load_enable,m,n,IFM_numxKxK,KxK,Ksize,TM_MIN,TN_MIN);
#else
	weight_load(Weight,weight_buffer,weight_load_enable,m,n,IFM_numxKxK,KxK,Ksize,TM_MIN,TN_MIN);
#endif

}

void copy_local_beta(float beta_buffer[MAX_BETA_LENGTH],float local_beta_buffer[MAX_BETA_LENGTH],const int TM_MIN,int m)
{
	int offset;
	int tm;
	for(tm = 0,offset = m;tm < TM_MIN;tm++)
	{
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tm)
#pragma HLS PIPELINE II=1
		local_beta_buffer[tm] = beta_buffer[offset];
		offset++;
	}
}

void compute(float input_buffer[Tn][OnChipIB_Height][OnChipIB_Width],float output_buffer[Tm][Tr][Tc],
		float weight_buffer[Tm][Tn][K][K],float beta_buffer[MAX_BETA_LENGTH],int n_next[1],
		const int Ksize,const int Kstride,int m,
		const int TM_MIN,const int TR_MIN,const int TC_MIN,bool enable)
{
	static float local_beta_buffer[Tm];
#pragma HLS ARRAY_PARTITION variable=local_beta_buffer complete dim=1

	if(!enable)
	{
		copy_local_beta(beta_buffer,local_beta_buffer,TM_MIN, m);
		return;
	}

	uint8_t i,j,tr,tc,tm,tn;
	int n = n_next[0];
	float partial_mul[Tm][Tn];
#pragma HLS ARRAY_PARTITION variable=partial_mul complete dim=1
#pragma HLS ARRAY_PARTITION variable=partial_mul complete dim=2
	float partial_add[Tm];
#pragma HLS ARRAY_PARTITION variable=partial_add complete dim=1

	for(i =0;i < Ksize; i++)
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=K)
		for(j = 0;j < Ksize; j++)
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=K)
			for(tr = 0;tr < TR_MIN;tr++)
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tr)
				for(tc = 0;tc < TC_MIN;tc++)
				{
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tc)
#pragma HLS PIPELINE II=3
					for(tm = 0;tm < Tm;tm++)
					{
#pragma HLS DEPENDENCE variable=output_buffer inter false

						if(i==0&&j==0&&n==0)
							partial_add[tm] = local_beta_buffer[tm];
						else
							partial_add[tm] = output_buffer[tm][tr][tc];

						for(tn = 0;tn <Tn;tn++)
						{
							partial_mul[tm][tn] = weight_buffer[tm][tn][i][j]*input_buffer[tn][Kstride*tr+i][Kstride*tc+j];
						}

						float partial_sum = 0;
						for(tn = 0;tn <Tn;tn++)
						{
							 partial_sum += partial_mul[tm][tn];
						}
						output_buffer[tm][tr][tc] = partial_add[tm] + partial_sum;
					}

				}
}

void nonlinear_leaky(float Input[Tm][Tr][Tc],uint8_t TM_MIN, uint8_t TR_MIN, uint8_t TC_MIN,const bool IsNL)
{
	uint8_t tr,tc,tm;

	assert((TM_MIN>0)&&(TM_MIN<=Tm));
	assert((TR_MIN>0)&&(TR_MIN<=Tr));
	assert((TC_MIN>0)&&(TC_MIN<=Tc));

	if(!IsNL)
		return ;

	float tmp_out;
	for(tm = 0;tm < TM_MIN;tm++)
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tm)
		for(tr = 0;tr < TR_MIN;tr++)
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tr)
			for(tc = 0;tc < TC_MIN;tc++)
			{
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tc)
#pragma HLS PIPELINE II=1
//#pragma HLS DEPENDENCE variable=Input inter false
				float tmp = Input[tm][tr][tc];
				if(tmp < 0.0f)
					tmp_out = tmp*0.1;
				else
					tmp_out = tmp;
				Input[tm][tr][tc] = tmp_out;
			}

}

/*
void write_back_output_reorg(float output_buffer[Tm][Tr][Tc],float *Output,int r,int c,int m,uint16_t Output_w,uint16_t Output_h,
		uint8_t TM_MIN,uint8_t TR_MIN,uint8_t TC_MIN,const int OHxOW, bool IsNL, bool write_flag)
{
	if(!write_flag)
		return;

	assert((TM_MIN >0)&&(TM_MIN <=Tm));
	assert((TR_MIN >0)&&(TR_MIN <=Tr));
	assert((TC_MIN >0)&&(TC_MIN <=Tc));

	const int offset = m*OHxOW + r*Output_w + c;
	uint8_t tr,tm;

	nonlinear_leaky(output_buffer,TM_MIN,TR_MIN,TC_MIN,IsNL);

	//for(tm = 0;tm < TM_MIN;tm++)
	//	for(tr = 0;tr < TR_MIN;tr++)
	//		for(tc = 0;tc < TC_MIN;tc++)
	//		{
	//				Output[tm*OHxOW + tr*Output_w + tc + offset] = output_buffer[tm][tr][tc];
	//		}

	for(tm = 0;tm < TM_MIN;tm++)
		for(tr = 0;tr < TR_MIN;tr++)
		{
			memcpy((float *)(Output + tm*OHxOW + tr*Output_w + offset),output_buffer[tm][tr],TC_MIN*sizeof(float));
		}
}
*/

void nonlinear_leaky_row(IO_Dtype output_localbuf[Tc], float Input[Tm][Tr][Tc], uint8_t tm, uint8_t tr, uint8_t *tm_n, uint8_t *tr_n, uint8_t TC_MIN,const bool IsNL, bool enable)
{
	if(!enable)
		return ;

	uint8_t tc;
	assert((TC_MIN>0)&&(TC_MIN<=Tc));

	float tmp_out;
	for(tc = 0;tc < TC_MIN;tc++)
	{
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tc)
#pragma HLS PIPELINE II=1
//#pragma HLS DEPENDENCE variable=Input inter false
		float tmp = Input[tm][tr][tc];
		if((tmp < 0.0f)&&IsNL)
			tmp_out = tmp*0.1f;
		else
			tmp_out = tmp;			
		output_localbuf[tc] = tmp_out;
	}
	
	*tm_n = tm;
	*tr_n = tr;
}

void ofm_mmcpy_row(IO_Dtype *Output, IO_Dtype local_buf[Tc], int offset, int OHxOW, int Output_w, int TC_MIN, uint8_t tm, uint8_t tr,bool enable)
{
	if(!enable)
		return;

	int ofm_offset = tm*OHxOW + tr*Output_w + offset;
	memcpy((IO_Dtype *)(Output + ofm_offset), local_buf, TC_MIN*sizeof(IO_Dtype));

	WriteByte += TC_MIN*sizeof(float);
}

void write_back_output_reorg(float output_buffer[Tm][Tr][Tc], IO_Dtype *Output,int r,int c,int m,uint16_t Output_w,uint16_t Output_h,
		uint8_t TM_MIN,uint8_t TR_MIN,uint8_t TC_MIN,const int OHxOW, bool IsNL, bool write_flag)
{
	if(!write_flag)
		return;

	assert((TM_MIN >0)&&(TM_MIN <=Tm));
	assert((TR_MIN >0)&&(TR_MIN <=Tr));
	assert((TC_MIN >0)&&(TC_MIN <=Tc));

	const int offset = m*OHxOW + r*Output_w + c;
	static IO_Dtype local_buf0[Tc];
	static IO_Dtype local_buf1[Tc];
	uint8_t tm_n0, tm_n1, tr_n0, tr_n1;

	bool pp = true;
	uint8_t tr,tm;
	uint16_t TM_MINxTR_MIN = TM_MIN*TR_MIN;
	uint16_t t;
	tr = 0, tm = 0;
	for(t = 0;t < TM_MINxTR_MIN + 1;t++)
	{
		if(pp)
		{
			nonlinear_leaky_row( local_buf0, output_buffer, tm, tr, &tm_n0, &tr_n0, TC_MIN, IsNL, t!=TM_MINxTR_MIN);
			ofm_mmcpy_row( Output, local_buf1, offset, OHxOW, Output_w, TC_MIN, tm_n1, tr_n1, t!=0);
			pp = false;
		}else
		{
			nonlinear_leaky_row( local_buf1, output_buffer, tm, tr, &tm_n1, &tr_n1, TC_MIN, IsNL, t!=TM_MINxTR_MIN);
			ofm_mmcpy_row( Output, local_buf0, offset, OHxOW, Output_w, TC_MIN, tm_n0, tr_n0, t!=0);
			pp = true;
		}

		tr++;
		if(tr==TR_MIN)
		{
			tr = 0;
			tm++;
		}
	}

/*
	for(tm = 0;tm < TM_MIN;tm++)
		for(tr = 0;tr < TR_MIN + 1;tr++)
		{
			if(pp)
			{
				nonlinear_leaky_row( local_buf0, output_buffer, tm, tr, &tm_n0, &tr_n0, TC_MIN, IsNL, tr!=TR_MIN);
				ofm_mmcpy_row( Output, local_buf1, offset, OHxOW, Output_w, TC_MIN, tm_n1, tr_n1, tr!=0);
				pp = false;
			}else
			{
				nonlinear_leaky_row( local_buf1, output_buffer, tm, tr, &tm_n1, &tr_n1, TC_MIN, IsNL, tr!=TR_MIN);
				ofm_mmcpy_row( Output, local_buf0, offset, OHxOW, Output_w, TC_MIN, tm_n0, tr_n0, tr!=0);
				pp = true;
			}
		}
*/
}

void pool_yolo2(float Input[Tn][OnChipIB_Height][OnChipIB_Width],float Output[Tm][Tr][Tc],
		  const int Ksize,const int Kstride,
		  const int TM_MIN,const int TR_MIN,const int TC_MIN,bool enable)
{
	if(!enable)
		return;

	uint8_t i,j,tr,tc,of;
	float tmp[Tn];
#pragma HLS ARRAY_PARTITION variable=tmp complete dim=1

	for(tr = 0;tr < TR_MIN;tr++)
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tr)
		for(tc = 0;tc < TC_MIN;tc++)
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tc)
			for(i =0;i < Ksize; i++)
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=K)
				for(j = 0;j < Ksize; j++)
				{
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=K)
#pragma HLS PIPELINE II=1
					for( of = 0; of < Tn; of++)
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

void reorg_yolo2(float Input[Tn][OnChipIB_Height][OnChipIB_Width],float Output[Tm][Tr][Tc],
		  const int Ksize,const int Kstride,
		  const int TM_MIN,const int TR_MIN,const int TC_MIN,bool enable)
{
	int x, y, kx, ky;
	unsigned char Yoffset;
	unsigned char Xoffset;

	if(!enable)
		return;

    for( y = 0; y < TR_MIN; y++)
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tr)
    	for( x = 0; x < TC_MIN; x++)
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tc)
		for(ky= 0;ky < 2; ky++)
		for(kx = 0;kx < 2; kx++)
		{
#pragma HLS PIPELINE II=1
			Yoffset = (y << 1) + ky;
			Xoffset = (x << 1) + kx;

			int in_index  = (ky << 1) + kx;
			Output[in_index][y][x] = Input[0][Yoffset][Xoffset];
		}
}

void intra_pingpong_wrapper(float *Input,float *Weight, float output_buffer[Tm][Tr][Tc],float beta_buffer[MAX_BETA_LENGTH],
								 float input_buffer0[Tn][OnChipIB_Height][OnChipIB_Width],float input_buffer1[Tn][OnChipIB_Height][OnChipIB_Width],
								 int IFM_num,int Input_w,int Input_h,int OFM_num,int Ksize,int Kstride,
								 int TMP_R,int TMP_C,int TMP_M,int TM_MIN,int TR_MIN,int TC_MIN,int TN,int TRow,int TCol,int Padding,
								 int IHxIW,int KxK,int IFM_numxKxK,int LayerType,int TM,int TMP_X_next[1],int TX_MIN_next[1],bool pingpongx,bool input_flag,bool process_flag)
{
	static float weight_buffer0[Tm][Tn][K][K];
#pragma HLS ARRAY_PARTITION variable=weight_buffer0 complete dim=1
#pragma HLS ARRAY_PARTITION variable=weight_buffer0 complete dim=2

	static float weight_buffer1[Tm][Tn][K][K];
#pragma HLS ARRAY_PARTITION variable=weight_buffer1 complete dim=1
#pragma HLS ARRAY_PARTITION variable=weight_buffer1 complete dim=2

	static int NOP[1];
	static int tmp_x;
	static int tmp_tx_min;

	if(LayerType==0)
	{
		if(!input_flag)
			return;

		TMP_X_next[0] = TMP_M;//consider by the inner-out loop
		TX_MIN_next[0] = TM_MIN;// like above

		bool pingpong = 0;
		int n0[1];
		int n1[1];
		int n;
		for(n = 0;n < IFM_num+TN; n += TN)
		{
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=2048)
			if(pingpong == 1)
			{
				copy_input_weight(Input,Weight,IFM_num,Input_w,Input_h,Ksize,Kstride,TMP_R,TMP_C,TMP_M, n,
					TM_MIN,TN,TRow,TCol,Padding,input_buffer1,weight_buffer1, n1, n < IFM_num,1,(TMP_M==0)&&(n==0),IHxIW,KxK,IFM_numxKxK,LayerType);
				compute(input_buffer0,output_buffer,weight_buffer0,beta_buffer, n0,Ksize,Kstride,TMP_M,TM_MIN,TR_MIN,TC_MIN, n!=0);
				pingpong = 0;
			}else
			{
				copy_input_weight(Input,Weight,IFM_num,Input_w,Input_h,Ksize,Kstride,TMP_R,TMP_C,TMP_M, n,
					TM_MIN,TN,TRow,TCol,Padding,input_buffer0,weight_buffer0, n0, n < IFM_num,1,(TMP_M==0)&&(n==0),IHxIW,KxK,IFM_numxKxK,LayerType);
				compute(input_buffer1,output_buffer,weight_buffer1,beta_buffer, n1,Ksize,Kstride,TMP_M,TM_MIN,TR_MIN,TC_MIN, n!=0);
				pingpong = 1;
			}
		}
	}
	else if(LayerType==1)
	{
		if(pingpongx==0)
		{
			TMP_X_next[0] = tmp_x;
			TX_MIN_next[0] = tmp_tx_min;
			tmp_x = TMP_M;
			tmp_tx_min = TM_MIN;

			copy_input_weight(Input,Weight,IFM_num,Input_w,Input_h,Ksize,Kstride,TMP_R,TMP_C,TMP_M,TMP_M,
				TM_MIN,TM,TRow,TCol,0,input_buffer0,weight_buffer0,NOP,input_flag,0,0,IHxIW,KxK,IFM_numxKxK,LayerType);
			pool_yolo2(input_buffer1,output_buffer,Ksize,Kstride,TX_MIN_next[0],TR_MIN,TC_MIN,process_flag);
		}else
		{
			TMP_X_next[0] = tmp_x;
			TX_MIN_next[0] = tmp_tx_min;
			tmp_x = TMP_M;
			tmp_tx_min = TM_MIN;

			copy_input_weight(Input,Weight,IFM_num,Input_w,Input_h,Ksize,Kstride,TMP_R,TMP_C,TMP_M,TMP_M,
				TM_MIN,TM,TRow,TCol,0,input_buffer1,weight_buffer1,NOP,input_flag,0,0,IHxIW,KxK,IFM_numxKxK,LayerType);
			pool_yolo2(input_buffer0,output_buffer,Ksize,Kstride,TX_MIN_next[0],TR_MIN,TC_MIN,process_flag);
		}

	}
	else if(LayerType==2)
	{
		if(pingpongx==0)
		{
			TMP_X_next[0] = tmp_x;
			TX_MIN_next[0] = tmp_tx_min;
			tmp_x = TMP_M;
			tmp_tx_min = TM_MIN;

			copy_input_weight(Input,Weight,IFM_num,Input_w,Input_h,Ksize,Kstride,TMP_R,TMP_C,TMP_M,TMP_M,
				TM_MIN,TM,TRow,TCol,0,input_buffer0,weight_buffer0,NOP,input_flag,0,0,IHxIW,KxK,IFM_numxKxK,LayerType);
			reorg_yolo2(input_buffer1,output_buffer,Ksize,Kstride,TX_MIN_next[0],TR_MIN,TC_MIN,process_flag);
		}else
		{
			TMP_X_next[0] = tmp_x;
			TX_MIN_next[0] = tmp_tx_min;
			tmp_x = TMP_M;
			tmp_tx_min = TM_MIN;

			copy_input_weight(Input,Weight,IFM_num,Input_w,Input_h,Ksize,Kstride,TMP_R,TMP_C,TMP_M,TMP_M,
				TM_MIN,TM,TRow,TCol,0,input_buffer1,weight_buffer1,NOP,input_flag,0,0,IHxIW,KxK,IFM_numxKxK,LayerType);
			reorg_yolo2(input_buffer0,output_buffer,Ksize,Kstride,TX_MIN_next[0],TR_MIN,TC_MIN,process_flag);
		}

	}

}

void YOLO2_FPGA(float *Input,float *Output,float *Weight,float *Beta, int IFM_num, int OFM_num,
							   int Ksize, int Kstride,
							   int Input_w, int Input_h, int Output_w, int Output_h, int Padding, bool IsNL, bool IsBN,
							   int TM, int TN, int TR, int TC,
							   int OFM_num_bound, int mLoopsxTM, int mLoops_a1xTM, int LayerType)
{
#pragma HLS INTERFACE m_axi depth=512 port=Input    offset=slave bundle=DATA_BUS num_read_outstanding=1 num_write_outstanding=1 max_read_burst_length=64 max_write_burst_length=64
#pragma HLS INTERFACE m_axi depth=512 port=Output    offset=slave bundle=DATA_BUS num_read_outstanding=1 num_write_outstanding=1 max_read_burst_length=64 max_write_burst_length=64
#pragma HLS INTERFACE m_axi depth=512 port=Weight offset=slave bundle=DATA_BUS1 num_read_outstanding=1 max_read_burst_length=128
#pragma HLS INTERFACE m_axi depth=512 port=Beta   offset=slave bundle=DATA_BUS1 num_read_outstanding=1 max_read_burst_length=128

#pragma HLS INTERFACE s_axilite register port=return bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=IFM_num bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=OFM_num bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=Ksize bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=Kstride bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=Input_w bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=Input_h bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=Output_w bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=Output_h bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=Padding bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=IsNL bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=IsBN bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=TM bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=TN bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=TR bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=TC bundle=CTRL_BUS

#pragma HLS INTERFACE s_axilite register port=OFM_num_bound bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=mLoopsxTM bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=mLoops_a1xTM bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=LayerType bundle=CTRL_BUS

#pragma HLS INTERFACE s_axilite register port=Input bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=Output bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=Weight bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=Beta bundle=CTRL_BUS

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
	assert((TR > 0)&&(TR <= Tr));
	assert((TC > 0)&&(TC <= Tc));

	const int OHxOW = Output_h*Output_w;
	const int TRow = (TR-1)*Kstride+Ksize;
	const int TCol = (TC-1)*Kstride+Ksize;
	const int IHxIW   = Input_h*Input_w;
	const int KxK = Ksize*Ksize;
	const int IFM_numxKxK = IFM_num*KxK;

	static float input_buffer0[Tn][OnChipIB_Height][OnChipIB_Width];
#pragma HLS ARRAY_PARTITION variable=input_buffer0 complete dim=1
	static float input_buffer1[Tn][OnChipIB_Height][OnChipIB_Width];
#pragma HLS ARRAY_PARTITION variable=input_buffer1 complete dim=1
	static float output_buffer[Tm][Tr][Tc];
#pragma HLS ARRAY_PARTITION variable=output_buffer complete dim=1
	static float output_buffer1[Tm][Tr][Tc];
#pragma HLS ARRAY_PARTITION variable=output_buffer1 complete dim=1
	static float beta_buffer[MAX_BETA_LENGTH];

/////////////////////////////////param
	int r, c, m;
	int TM_MIN,TR_MIN,TC_MIN;
///////////////////////////////////////

	int m0[1], m1[1];
	int TM_MIN0[1], TM_MIN1[1];
	bool pingpongm;

	if(LayerType==0)
		memcpy(beta_buffer,Beta, OFM_num*sizeof(float));

	for(r = 0; r < Output_h; r += TR)
	{
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=1024)
		TR_MIN = MIN(TR,Output_h-r);
		for(c = 0; c < Output_w; c += TC)
		{
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=1024)
			TC_MIN = MIN(TC,Output_w-c);
			pingpongm = 0;
			for(m = 0; m < OFM_num_bound; m += TM)
			{
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=2048)
				TM_MIN = MIN(TM, OFM_num-m);
				bool Mne0 = (m!=0);
				bool Mne1 = (m!=TM);
				bool MnemLps = (m!=mLoopsxTM);
				bool MneMLps_a1 = (m!=mLoops_a1xTM);
				bool input_flag = LayerType ? MnemLps&&MneMLps_a1: MnemLps;
				bool process_flag = LayerType ? Mne0&&MneMLps_a1 : MnemLps;
				bool write_flag = LayerType ? Mne0&&Mne1 : Mne0;

				if(pingpongm==0)
				{
					intra_pingpong_wrapper(Input,Weight,output_buffer1,beta_buffer,input_buffer0,input_buffer1,
									IFM_num, Input_w, Input_h, OFM_num, Ksize, Kstride,
									r, c, m, TM_MIN, TR_MIN, TC_MIN, TN, TRow, TCol, Padding,IHxIW,KxK,IFM_numxKxK,LayerType,TM, m1,TM_MIN1, pingpongm, input_flag, process_flag);

					write_back_output_reorg(output_buffer,Output, r, c, m0[0],Output_w,Output_h, TM_MIN0[0], TR_MIN, TC_MIN, OHxOW, IsNL, write_flag);
					pingpongm = 1;
				}else
				{
					intra_pingpong_wrapper(Input,Weight,output_buffer,beta_buffer,input_buffer0,input_buffer1,
									IFM_num, Input_w, Input_h, OFM_num, Ksize, Kstride,
									r, c, m, TM_MIN, TR_MIN, TC_MIN, TN, TRow, TCol, Padding,IHxIW,KxK,IFM_numxKxK,LayerType,TM, m0,TM_MIN0, pingpongm, input_flag, process_flag);

					write_back_output_reorg(output_buffer1,Output, r, c, m1[0],Output_w,Output_h, TM_MIN1[0], TR_MIN, TC_MIN, OHxOW, IsNL, write_flag);
					pingpongm = 0;
				}

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

			int TN_MINxTM_MIN = TN_MIN*TM_MIN;
			for(tk = 0;tk < KxK; tk++)
				for(tm = 0;tm < TM_MIN; tm++)
					for(tn = 0;tn < TN_MIN;tn++)
					{
						weight_buffer2[tk*TN_MINxTM_MIN + tm*TN_MIN + tn] = weight_buffer[tm*TN_MIN*KxK + tn*KxK + tk];
					}

			memcpy((float *)(Weight_reorg+offset),weight_buffer2,TM_MIN*TN_MIN*KxK*sizeof(float));
			offset += TM_MIN*TN_MIN*KxK;
		}							
	}

	return 0;
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

	float *Weight_buf = (float *)calloc(203767168/4,sizeof(float));
	float *Beta_buf   = (float *)calloc(43044/4,sizeof(float));

#ifdef REORG_TEST
	FILE *fp_w = fopen("weights_reorg.bin", "rb");
    	if(!fp_w) file_error("weights_reorg.bin");
#else
	FILE *fp_w = fopen("weights.bin", "rb");
    	if(!fp_w) file_error("weights.bin");
#endif

#ifdef REORG_GEN
	float *Weight_reorg_buf = (float *)calloc(203767168/4,sizeof(float));
	FILE *fp_w_reorg = fopen("weights_reorg.bin", "wb");
    	if(!fp_w_reorg) file_error("weights_reorg.bin");
#endif
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

	memcpy(in_ptr[0], input, 416*416*3*sizeof(float));//416x416x3 input_pic

    	int i;
	int offset_index = 0;
	int woffset = 0;
	int boffset = 0;
	int TR,TC,TM,TN;
	int output_w,output_h;
	int mLoops;
	double sum_gop = 0.0;
	double sum_gop_max_2x2 = 0.0;
	double sum_gop_conv2d_3x3 = 0.0;
	double sum_gop_conv2d_1x1 = 0.0;
	uint32_t sum_params_conv2d_3x3 = 0;
	uint32_t sum_params_conv2d_1x1 = 0;

	double sum_conv2d_3x3_rd, sum_conv2d_3x3_wr;
	double sum_conv2d_1x1_rd, sum_conv2d_1x1_wr;
	double sum_maxpool_2x2_rd, sum_maxpool_2x2_wr;
	double reorg_rd, reorg_wr;
	sum_conv2d_3x3_rd = sum_conv2d_3x3_wr = 0;
	sum_conv2d_1x1_rd = sum_conv2d_1x1_wr = 0;
	sum_maxpool_2x2_rd = sum_maxpool_2x2_wr = 0;

	ReadByte = 0;
	WriteByte = 0;

    	for(i = 0; i < net->n; ++i)
	{
        	layer l = net->layers[i];
		printf("Layer[%2d]: ",i);
		switch(l.type)
		{
			case CONVOLUTIONAL:
				printf("outputMemory:%8d;BN=%d;Activation=%d;conv  %5d %2d x%2d /%2d  %4d x%4d x%4d   ->  %4d x%4d x%4d  %5.3f BFLOPs\n",l.outputs,l.batch_normalize,l.activation, l.n, l.size, l.size, l.stride, l.w, l.h, l.c, l.out_w, l.out_h, l.out_c, (2.0 * l.n * l.size*l.size*l.c/l.groups * l.out_h*l.out_w)/1000000000.);
				sum_gop += (2.0 * l.n * l.size*l.size*l.c/l.groups * l.out_h*l.out_w)/1000000000.;
				if(l.size==1)
				{
					sum_params_conv2d_1x1 += (l.out_c + l.n * l.size*l.size*l.c);
					sum_gop_conv2d_1x1 += (2.0 * l.n * l.size*l.size*l.c * l.out_h*l.out_w)/1000000000.;
					sum_conv2d_1x1_rd += (3.0 * l.n * l.size*l.size*l.c * l.out_h*l.out_w)/1000000000.;
					sum_conv2d_1x1_wr += (1.0 * l.n * l.size*l.size*l.c * l.out_h*l.out_w)/1000000000.;
				}else if(l.size==3)
				{
					sum_gop_conv2d_3x3 += (2.0 * l.n * l.size*l.size*l.c/l.groups * l.out_h*l.out_w)/1000000000.;
					sum_params_conv2d_3x3 += (l.out_c + l.n * l.size*l.size*l.c);
					sum_conv2d_3x3_rd += (3.0 * l.n * l.size*l.size*l.c * l.out_h*l.out_w)/1000000000.;
					sum_conv2d_3x3_wr += (1.0 * l.n * l.size*l.size*l.c * l.out_h*l.out_w)/1000000000.;
				}

				output_w = (l.w - l.size + 2*l.pad)/l.stride + 1;
				output_h = (l.h - l.size + 2*l.pad)/l.stride + 1;

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
#ifdef REORG_GEN
				Weight_reorgnaization_anti(Weight_buf + woffset,Weight_reorg_buf + woffset,NULL,l.c,l.n,l.size,TM,TN,0);
#endif
				woffset += weight_offset[offset_index];
				boffset += beta_offset[offset_index];
				offset_index++;

				break;
			case MAXPOOL:
				printf("outputMemory:%8d;max          %d x %d / %d  %4d x%4d x%4d   ->  %4d x%4d x%4d\n",l.outputs, l.size, l.size, l.stride, l.w, l.h, l.c, l.out_w, l.out_h, l.out_c);
				sum_gop += (4.0 * l.c * l.out_h*l.out_w)/1000000000.;				
				sum_gop_max_2x2 += (4.0 * l.c * l.out_h*l.out_w)/1000000000.;
				sum_maxpool_2x2_rd += (4.0 * l.c * l.out_h*l.out_w)/1000000000.;
				sum_maxpool_2x2_wr += (1.0 * l.c * l.out_h*l.out_w)/1000000000.;
				//output_w = (l.w - l.size)/l.stride + 1 ;
				//output_h = (l.h - l.size)/l.stride + 1 ;
				output_w = l.out_h;
				output_h = l.out_w;

				TR = MIN(((OnChipIB_Height-l.size)/l.stride+1),Tr);//keep Kstride>=1
				TC = MIN(((OnChipIB_Width-l.size)/l.stride+1),Tc);
				TR = MIN(output_h,TR);
				TC = MIN(output_w,TC);
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
				reorg_rd = reorg_wr = 4*26*32*13/1000000000.;

				TR = MIN(((OnChipIB_Height-l.stride)/l.stride+1),Tr);//keep Kstride>=1
				TR = MIN(output_h,TR);
				TC = MIN(((OnChipIB_Width-l.stride)/l.stride+1),Tc);
				TC = MIN(output_w,TC);
				//TM = 4;
				//mLoops = 1;
				//here Tm and Tn must >=4 
				TM = MIN(Tm,Tn);
				TM = MIN(4,TM);
				mLoops = (int)ceil(((float)4)/TM);

				YOLO2_FPGA(in_ptr[i],out_ptr[i],NULL,NULL,1,4,
							  l.stride,l.stride,52,32*26, output_w, output_h, 0,0,0,TM,0,TR,TC, (mLoops + 2)*TM, mLoops*TM, (mLoops + 1)*TM, 2);

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
	printf("SUM_GOP=%.7lf, conv2d_3x3=%.7lf, conv2d_1x1=%.7lf, max_2x2=%.7lf\n",sum_gop, sum_gop_conv2d_3x3, sum_gop_conv2d_1x1, sum_gop_max_2x2);
	printf("SUM_params = %d, conv2d_3x3= %d, conv2d_1x1 %d\n", sum_params_conv2d_3x3 + sum_params_conv2d_1x1, sum_params_conv2d_3x3, sum_params_conv2d_1x1);
	printf("MA conv2d_3x3 rd= %.7lf, wr= %.7lf\n", sum_conv2d_3x3_rd, sum_conv2d_3x3_wr);
	printf("MA conv2d_1x1 rd= %.7lf, wr= %.7lf\n", sum_conv2d_1x1_rd, sum_conv2d_1x1_wr);
	printf("MA maxpool_2x2 rd= %.7lf, wr= %.7lf\n", sum_maxpool_2x2_rd, sum_maxpool_2x2_wr);
	printf("MA reorg rd= %.7lf, wr= %.7lf\n", reorg_rd, reorg_wr);

	printf("ReadByte = %d, WriteByte = %d, Sum = %d\n", ReadByte, WriteByte, ReadByte + WriteByte);
	printf("Operational Intensity = %g\n", sum_gop/((ReadByte + WriteByte)/1000000000.0));

#ifdef REORG_GEN
	fwrite(Weight_reorg_buf, sizeof(float), 203767168/4, fp_w_reorg);
	fclose(fp_w_reorg);
	free(Weight_reorg_buf);
#endif
	free(Memory_buf);
	free(Weight_buf);
	free(Beta_buf);

}
///////////////////////////////////////////////////////////////////////20181229 anti-reorg ok end n4m32
