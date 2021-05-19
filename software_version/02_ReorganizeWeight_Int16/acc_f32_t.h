
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#define MIN_diy(x,y) ((x) < (y) ? (x) : (y))
#define MAX_diy(x,y) ((x) > (y) ? (x) : (y))

#define FALSE 0
#define TRUE  1

#define VALID 0
#define SAME  1

#define LT_CONV  0
#define LT_MAXPOOL 1

#define MIN_NEG (-1024*1024)

#DEFINE_HEADER#

/*
#define HW_S 3
#define K 7
#define Tn 2
#define Tm 32
#define Tr 28
#define Tc 32
#define MAX_BETA_LENGTH 2048
#define REORDER_TEST

#define OnChipIB_Width  ((Tc-1)*HW_S+K)
#define OnChipIB_Height ((Tr-1)*HW_S+K)
*/

#define PRAGMA_SUB(x) _Pragma (#x)
#define DO_PRAGMA(x) PRAGMA_SUB(x)

void ifm_mmcpy_row(float *input, float input_memcpy_buffer[OnChipIB_Width], int CurrentOffset, uint32_t IHxIW, uint16_t Input_w, uint16_t TCol, 
	uint8_t t1, uint8_t t2, uint8_t *t1_n, uint8_t *t2_n,bool enable)
{
	if(!enable)
		return;

	int ifm_offset = CurrentOffset + t1*IHxIW + t2*Input_w;
	memcpy( input_memcpy_buffer,(float *)(input + ifm_offset), TCol*sizeof(float));

	t1_n[0] = t1;
	t2_n[0] = t2;
}

void ifm_copy_lbuf2ibuf(float input_buffer[Tn][OnChipIB_Height][OnChipIB_Width], float local_buf[OnChipIB_Width], uint16_t TCol, uint16_t Input_w, uint16_t Input_h, uint8_t TN_MIN, float pad_val,
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
//DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=TCol_max)
#pragma HLS PIPELINE II=1
		int xoffset = Coffset + t3;
		bool XEnable = (xoffset >= 0)&&(xoffset < Input_w);
		if(XEnable&&PEnable)
		{
			input_buffer[t1][t2][t3] = local_buf[t3];
		}
		else
			input_buffer[t1][t2][t3] = pad_val;
	}
}

void input_load(float *input,float input_buffer[Tn][OnChipIB_Height][OnChipIB_Width], uint16_t r, uint16_t c, uint16_t n, uint8_t Kstride, uint8_t Padding, uint16_t TRow, uint16_t TCol,
		 uint16_t Input_w, uint16_t Input_h, uint8_t TN_MIN, uint32_t IHxIW, float pad_val)
{
	uint8_t t1,t2;
	static float local_buf0[OnChipIB_Width];
	static float local_buf1[OnChipIB_Width];

	const int Coffset = c*Kstride - Padding;
	const int Roffset = r*Kstride - Padding;
	const int CurrentOffset = n*IHxIW + Roffset*Input_w + Coffset;

	uint8_t t1_n0[1], t1_n1[1], t2_n0[1], t2_n1[1];
	bool pp = true;

	uint32_t TnxTRow;
//	if(IsNotConv)
	TnxTRow = TN_MIN*TRow;//dont need fill all ifm buf, because result is not cross-channel
//	else
//		TnxTRow = Tn*TRow;
	uint32_t TnxTRow_a1 = TnxTRow + 1;

	uint32_t t = 0;
	t1 = 0; t2 = 0;
	for(t = 0;t < TnxTRow_a1; t++)
	{
//DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=TRow_max+1)
		if(pp)
		{
			ifm_mmcpy_row(input, local_buf0, CurrentOffset, IHxIW, Input_w, TCol, t1, t2, t1_n0, t2_n0, t!=TnxTRow);
			ifm_copy_lbuf2ibuf( input_buffer, local_buf1, TCol, Input_w, Input_h, TN_MIN, pad_val, Coffset, Roffset, t1_n1[0], t2_n1[0], t!=0);
			pp = false;
		}else
		{
			ifm_mmcpy_row(input, local_buf1, CurrentOffset, IHxIW, Input_w, TCol, t1, t2, t1_n1, t2_n1, t!=TnxTRow);
			ifm_copy_lbuf2ibuf( input_buffer, local_buf0, TCol, Input_w, Input_h, TN_MIN, pad_val, Coffset, Roffset, t1_n0[0], t2_n0[0], t!=0);
			pp = true;
		}
		
		t2++;
		if(t2==TRow)
		{
			t2 = 0;
			t1++;
		}
	}

}

void weight_load(float *Weight,float weight_buffer[Tm][Tn][K][K],int m,int n, uint32_t IFM_numxKxK, uint8_t KxK, uint8_t Ksize, uint8_t TM_MIN, uint8_t TN_MIN, bool enable)
{
	if(!enable)
		return;

	uint8_t t1,t2,t3,t4;
	static float weight_memcpy_buffer[Tm*Tn*K*K];

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

void weight_load_reorg(float *Weight,float weight_buffer[Tm][Tn][K][K],int m,int n, uint32_t IFM_numxKxK, uint8_t KxK, uint8_t Ksize, uint8_t TM_MIN, uint8_t TN_MIN, bool enable)
{
	if(!enable)
		return;

	uint8_t t1,t2,t3,t4;
	static float weight_memcpy_buffer[Tm*Tn*K*K];
	static int Woffset;

	assert((TM_MIN > 0)&&(TM_MIN <= Tm));
	assert((TN_MIN > 0)&&(TN_MIN <= Tn));
	assert((KxK > 0)&&(KxK <= K*K));

	if(m==0&&n==0)
		Woffset = 0;

	uint16_t mm_offset = TM_MIN*TN_MIN*KxK;
	memcpy(weight_memcpy_buffer,(float *)(Weight + Woffset), mm_offset*sizeof(float));
	Woffset += mm_offset;

	int weight_memcpy_offset = 0;
	float input_value;
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

void nonlinear_leaky_row(float output_localbuf[Tc], float Input[Tm][Tr][Tc], float local_beta_buffer[Tm], uint8_t tm, uint8_t tr, uint8_t *tm_n, uint8_t *tr_n,
			 uint8_t TC_MIN, bool LoadBias, bool IsNL, bool enable)
{
	if(!enable)
		return ;

	uint8_t tc;
	assert((TC_MIN>0)&&(TC_MIN<=Tc));

	float tmp_out, tmp, tmp1;
	for(tc = 0;tc < TC_MIN;tc++)
	{
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tc)
#pragma HLS PIPELINE II=1
		tmp = Input[tm][tr][tc];
		if(LoadBias)
			tmp1 = tmp + local_beta_buffer[tm];
		else
			tmp1 = tmp;

		if(IsNL)
		{
			if(tmp1 < 0.0f)
				tmp_out = tmp1*0.1f;
			else
				tmp_out = tmp1;
		}
		else
			tmp_out = tmp1;
		output_localbuf[tc] = tmp_out;
	}
	
	tm_n[0] = tm;
	tr_n[0] = tr;
}

void ofm_mmcpy_row(float *Output, float local_buf[Tc], uint32_t offset, uint32_t OHxOW, uint16_t Output_w, uint8_t TC_MIN, uint8_t tm, uint8_t tr,bool enable)
{
	if(!enable)
		return;

	uint32_t ofm_offset = tm*OHxOW + tr*Output_w + offset;
	memcpy((float *)(Output + ofm_offset), local_buf, TC_MIN*sizeof(float));
}

void copy_local_beta(float beta_buffer[MAX_BETA_LENGTH],float *local_beta_buffer, uint8_t TM_MIN, uint16_t m)
{
	int offset = m;
	int tm;
	for(tm = 0;tm < TM_MIN;tm++)
	{
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tm)
#pragma HLS PIPELINE II=1
		local_beta_buffer[tm] = beta_buffer[offset];
		offset++;
	}
}

void write_back_output_reorg(float output_buffer[Tm][Tr][Tc], float bias_buffer[MAX_BETA_LENGTH], float *Output, uint16_t r,uint16_t c,uint16_t m, uint16_t Output_w,
		uint8_t TM_MIN, uint8_t TR_MIN, uint8_t TC_MIN, uint32_t OHxOW, bool IsNL, bool LoadBias, bool enable)
{
	assert(Tn <= Tm);//for local beta buffer

	static float local_beta_buffer[Tm];
	static float local_buf0[Tc];
	static float local_buf1[Tc];

	if(!enable)
		return;

	assert((TM_MIN >0)&&(TM_MIN <=Tm));
	assert((TR_MIN >0)&&(TR_MIN <=Tr));
	assert((TC_MIN >0)&&(TC_MIN <=Tc));

	if(LoadBias)
	{
		copy_local_beta( bias_buffer, local_beta_buffer, TM_MIN, m);
	}

	uint32_t offset = m*OHxOW + r*Output_w + c;
	uint16_t TM_MINxTR_MIN = TM_MIN*TR_MIN;
	uint16_t TM_MINxTR_MIN_a1 = TM_MINxTR_MIN + 1;
	
	uint8_t tm_n0[1], tm_n1[1], tr_n0[1], tr_n1[1];
	bool pp = true;
	uint8_t tr,tm;
	uint16_t t;
	tr = 0, tm = 0;
	for(t = 0;t < TM_MINxTR_MIN_a1;t++)//flatten embedded loops
	{
		if(pp)
		{
			nonlinear_leaky_row( local_buf0, output_buffer, local_beta_buffer, tm, tr, tm_n0, tr_n0, TC_MIN, LoadBias, IsNL, t!=TM_MINxTR_MIN);
			ofm_mmcpy_row( Output, local_buf1, offset, OHxOW, Output_w, TC_MIN, tm_n1[0], tr_n1[0], t!=0);
			pp = false;
		}else
		{
			nonlinear_leaky_row( local_buf1, output_buffer, local_beta_buffer, tm, tr, tm_n1, tr_n1, TC_MIN, LoadBias, IsNL, t!=TM_MINxTR_MIN);
			ofm_mmcpy_row( Output, local_buf0, offset, OHxOW, Output_w, TC_MIN, tm_n0[0], tr_n0[0], t!=0);
			pp = true;
		}

		tr++;
		if(tr==TR_MIN)
		{
			tr = 0;
			tm++;
		}
	}

}

void maxpool_tile(float Input[Tn][OnChipIB_Height][OnChipIB_Width],float Output[Tm][Tr][Tc],
		  uint8_t Ksize, uint8_t Kstride, uint8_t TM_MIN, uint8_t TR_MIN, uint8_t TC_MIN, bool enable)
{
	if(!enable)
		return;

	uint8_t i,j,tr,tc,of;
	float tmp[Tn];
#pragma HLS ARRAY_PARTITION variable=tmp complete dim=1

	for(i =0;i < Ksize; i++)
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=K)
		for(j = 0;j < Ksize; j++)
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=K)
			for(tr = 0;tr < TR_MIN;tr++)
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tr)
				for(tc = 0;tc < TC_MIN;tc++)
				{
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tc)
#pragma HLS PIPELINE II=1
					for( of = 0; of < Tn; of++)
					{
						if(i==0&&j==0)
							tmp[of] = MIN_NEG;
						else
							tmp[of] = Output[of][tr][tc];

						float tmp_in = Input[of][tr*Kstride+i][tc*Kstride+j];

						if(tmp_in > tmp[of])
							tmp[of] = tmp_in;

						Output[of][tr][tc] = tmp[of];
					}
				}

}

void conv2d_tile(float input_buffer[Tn][OnChipIB_Height][OnChipIB_Width],float output_buffer[Tm][Tr][Tc],
		float weight_buffer[Tm][Tn][K][K],int n_next,
		const int Ksize,const int Kstride,int m,
		const int TM_MIN,const int TR_MIN,const int TC_MIN,bool enable)
{
	if(!enable)
	{
		return;
	}

	uint8_t i,j,tr,tc,tm,tn;
	const bool ne0 = (n_next==0);

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

						if(i==0&&j==0&&ne0)
							partial_add[tm] = 0;
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

void ifm_weight_load_wrapper(float *ifm, float *weight, float ifm_buffer[Tn][OnChipIB_Height][OnChipIB_Width], float weight_buffer[Tm][Tn][K][K], 
	int ifm_num, int ofm_num, int ifm_w, int ifm_h, int pad_int,
	int TM, int TN, int tr, int tc, int tm, int tn, int tn_next[1], int ksize, int kstride, int ltype, int TM_MIN, int TRow, int TCol, int IHW, int INumxKK, int KK,
	bool enable, float pad_val, bool LoadBias)
{
	if(!enable)
		return;
	
	tn_next[0] = tn;

	int TN_MIN, TMN_MIN;
	int tmn;
	if(ltype == LT_CONV)
	{
		tmn = tn;
		TN_MIN = MIN_diy(TN, ifm_num-tn);
		TMN_MIN = TN_MIN;
	}
	else
	{
		tmn = tm;
		TN_MIN = 1;//DCONV
		TMN_MIN = TM_MIN;
	}

	input_load(ifm, ifm_buffer, tr, tc, tmn, kstride, pad_int, TRow, TCol, ifm_w, ifm_h, TMN_MIN, IHW, pad_val);

#ifdef REORDER_TEST
	weight_load_reorg(weight,weight_buffer,tm,tn,INumxKK,KK,ksize,TM_MIN,TN_MIN, LoadBias);
#else
	      weight_load(weight,weight_buffer,tm,tn,INumxKK,KK,ksize,TM_MIN,TN_MIN, LoadBias);
#endif

}

void load_compute_wrapper(float *ifm, float *weight, float ofm_buffer[Tm][Tr][Tc], int ksize, int kstride, int ifm_num, int ifm_w, int ifm_h,int ofm_num,
	 int pad_int, int ltype, int TRow, int TCol, int IHW, int KK, int INumxKK, int TC_MIN, int TR_MIN, int TM_MIN, int TM, int TN, int tm, int tr, int tc,
	 int TMP_X_next[1],int TX_MIN_next[1],bool pingpongx,bool input_flag,bool process_flag, float pad_val, bool LoadBias)
{
	static float ifm_buffer0[Tn][OnChipIB_Height][OnChipIB_Width];
#pragma HLS ARRAY_PARTITION variable=ifm_buffer0 complete dim=1
	static float ifm_buffer1[Tn][OnChipIB_Height][OnChipIB_Width];
#pragma HLS ARRAY_PARTITION variable=ifm_buffer1 complete dim=1
	static float weight_buffer0[Tm][Tn][K][K];
#pragma HLS ARRAY_PARTITION variable=weight_buffer0 complete dim=1
#pragma HLS ARRAY_PARTITION variable=weight_buffer0 complete dim=2
	static float weight_buffer1[Tm][Tn][K][K];
#pragma HLS ARRAY_PARTITION variable=weight_buffer1 complete dim=1
#pragma HLS ARRAY_PARTITION variable=weight_buffer1 complete dim=2

	static int n0[1];
	static int n1[1];

	static int tmp_x;
	static int tmp_tx_min;

	if(ltype == LT_CONV)
	{
		if(!input_flag)
			return;

		TMP_X_next[0] = tm;//consider by the inner-out loop
		TX_MIN_next[0] = TM_MIN;// like above

		bool pingpong = false;
		for(int tn = 0; tn < ifm_num+TN; tn += TN)
		{
			if(pingpong)
			{
				ifm_weight_load_wrapper(ifm, weight, ifm_buffer1, weight_buffer1, ifm_num, ofm_num, ifm_w, ifm_h,
					pad_int, TM, TN, tr, tc, tm, tn, n1, ksize, kstride, ltype, TM_MIN, TRow, TCol, IHW, INumxKK, KK, tn < ifm_num, pad_val, LoadBias);
				conv2d_tile(ifm_buffer0, ofm_buffer, weight_buffer0, n0[0], ksize, kstride, tm, TM_MIN, TR_MIN, TC_MIN, tn!=0);
				pingpong = false;
			}else
			{
				ifm_weight_load_wrapper(ifm, weight, ifm_buffer0, weight_buffer0, ifm_num, ofm_num, ifm_w, ifm_h,
					pad_int, TM, TN, tr, tc, tm, tn, n0, ksize, kstride, ltype, TM_MIN, TRow, TCol, IHW, INumxKK, KK, tn < ifm_num, pad_val, LoadBias);
				conv2d_tile(ifm_buffer1, ofm_buffer, weight_buffer1, n1[0], ksize, kstride, tm, TM_MIN, TR_MIN, TC_MIN, tn!=0);
				pingpong = true;
			}
		}
	}else if(ltype == LT_MAXPOOL)
	{
		if(pingpongx)
		{
			ifm_weight_load_wrapper(ifm, weight, ifm_buffer0, weight_buffer0, ifm_num, ofm_num, ifm_w, ifm_h,
				pad_int, TM, TN, tr, tc, tm, 0, n0, ksize, kstride, ltype, TM_MIN, TRow, TCol, IHW, INumxKK, KK, input_flag, pad_val, LoadBias);
			maxpool_tile(ifm_buffer1, ofm_buffer, ksize, kstride, tmp_tx_min, TR_MIN, TC_MIN, process_flag);

			TMP_X_next[0] = tmp_x;
			TX_MIN_next[0] = tmp_tx_min;
			tmp_x = tm;
			tmp_tx_min = TM_MIN;
		}else
		{
			ifm_weight_load_wrapper(ifm, weight, ifm_buffer1, weight_buffer1, ifm_num, ofm_num, ifm_w, ifm_h,
				pad_int, TM, TN, tr, tc, tm, 0, n1, ksize, kstride, ltype, TM_MIN, TRow, TCol, IHW, INumxKK, KK, input_flag, pad_val, LoadBias);
			maxpool_tile(ifm_buffer0, ofm_buffer, ksize, kstride, tmp_tx_min, TR_MIN, TC_MIN, process_flag);

			TMP_X_next[0] = tmp_x;
			TX_MIN_next[0] = tmp_tx_min;
			tmp_x = tm;
			tmp_tx_min = TM_MIN;
		}
	}
}

void FPGA_Acc(float *ifm, float *ofm, float *weight, float *bias, uint32_t k_s_pad_ltype, uint32_t iofm_num, uint32_t ifm_w_h, uint32_t ofm_w_h,
	uint32_t TRTC, uint32_t TMTN, int32_t OFM_num_bound, int32_t mLoopsxTM, int32_t mLoops_a1xTM, float pad_val, uint32_t TRowTCol,
	uint32_t IHW, uint32_t OHW, uint32_t KK_INumxKK, uint32_t en_bits)//enable_bits[2:0]={IsReLU, LoadBias, IsNotConv}
{
#pragma HLS INTERFACE m_axi depth=512 port=ifm    offset=slave bundle=DATA_BUS num_read_outstanding=1 num_write_outstanding=1 max_read_burst_length=64 max_write_burst_length=64
#pragma HLS INTERFACE m_axi depth=512 port=ofm    offset=slave bundle=DATA_BUS num_read_outstanding=1 num_write_outstanding=1 max_read_burst_length=64 max_write_burst_length=64
#pragma HLS INTERFACE m_axi depth=512 port=weight offset=slave bundle=DATA_BUS1 num_read_outstanding=1 max_read_burst_length=128
#pragma HLS INTERFACE m_axi depth=512 port=bias   offset=slave bundle=DATA_BUS1 num_read_outstanding=1 max_read_burst_length=128

#pragma HLS INTERFACE s_axilite register port=return bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=k_s_pad_ltype bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=iofm_num bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=ifm_w_h bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=ofm_w_h bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=TRTC bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=TMTN bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=OFM_num_bound bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=mLoopsxTM bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=mLoops_a1xTM bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=pad_val bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=TRowTCol bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=IHW bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=OHW bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=KK_INumxKK bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=en_bits bundle=CTRL_BUS

#pragma HLS INTERFACE s_axilite register port=ifm bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=ofm bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=weight bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=bias bundle=CTRL_BUS

	uint16_t ifm_w = (ifm_w_h >> 16) & 0xFFFF;
	uint16_t ifm_h = ifm_w_h & 0xFFFF;
	uint16_t ofm_w = (ofm_w_h >> 16) & 0xFFFF;
	uint16_t ofm_h = ofm_w_h & 0xFFFF;
	uint16_t TR = (TRTC >> 16) & 0xFFFF;
	uint16_t TC = TRTC & 0xFFFF;
	uint16_t TM = (TMTN >> 16) & 0xFFFF;
	uint16_t TN = TMTN & 0xFFFF;
	uint16_t ifm_num = (iofm_num >> 16) & 0xFFFF;
	uint16_t ofm_num = iofm_num & 0xFFFF;
	uint8_t ksize = (k_s_pad_ltype >> 24) & 0xFF;
	uint8_t kstride = (k_s_pad_ltype >> 16) & 0xFF;
	uint8_t pad_int = (k_s_pad_ltype >> 8) & 0xFF;
	uint8_t ltype = k_s_pad_ltype & 0xFFFF;
	uint16_t TRow = (TRowTCol >> 16) & 0xFFFF;
	uint16_t TCol = TRowTCol & 0xFFFF;
	uint8_t KK = (KK_INumxKK >> 24) & 0xFF;
	uint32_t INumxKK = KK_INumxKK & 0xFFFFFF;//24bit
	bool IsReLU    = (en_bits >> 2) & 0x1;
	bool LoadBias  = (en_bits >> 1) & 0x1;
	bool IsNotConv =  en_bits       & 0x1;

	assert((ofm_num > 0)&&(ofm_num <= 4096));
	assert((ifm_num > 0)&&(ifm_num < 65536));
	assert((kstride > 0)&&(kstride <= HW_S));
	assert((ksize > 0)&&(ksize < 8));//maybe include some pool ops
	assert((ifm_w > 0)&&(ifm_w <= 2048));
	assert((ifm_h > 0)&&(ifm_h <= 2048));
	assert((ofm_w > 0)&&(ofm_w <= 2048));
	assert((ofm_h > 0)&&(ofm_h <= 2048));
	assert((pad_int >= 0)&&(pad_int <= 4));//maybe
	assert((TM > 0)&&(TM <= Tm));
	assert((TN >= 0)&&(TN <= Tn));
	assert((TR > 0)&&(TR <= Tr));
	assert((TC > 0)&&(TC <= Tc));

	static float ofm_buffer0[Tm][Tr][Tc];
#pragma HLS ARRAY_PARTITION variable=ofm_buffer0 complete dim=1
	static float ofm_buffer1[Tm][Tr][Tc];
#pragma HLS ARRAY_PARTITION variable=ofm_buffer1 complete dim=1
	static float bias_buffer[MAX_BETA_LENGTH];
	
	if(LoadBias)
		memcpy(bias_buffer, bias, ofm_num*sizeof(float));

	int tr, tc, tm;
	int TR_MIN, TC_MIN, TM_MIN;

	int m0[1], m1[1];
	int TM_MIN0[1], TM_MIN1[1];
	bool pingpongm;

	for(tr = 0; tr < ofm_h; tr += TR)
	{
		TR_MIN = MIN_diy(TR,ofm_h-tr);
		for(tc = 0; tc < ofm_w; tc += TC)
		{
		    TC_MIN = MIN_diy(TC,ofm_w-tc);
		    pingpongm = false;
		    for(tm = 0; tm < OFM_num_bound; tm += TM)
		    {
			TM_MIN = MIN_diy(TM,ofm_num-tm);
			bool input_flag =  (tm < mLoopsxTM);
			bool process_flag = (tm > 0)&&(tm < mLoops_a1xTM);//CONV dont use, the same as input_flag
			bool write_flag = IsNotConv ? (tm > TM) : (tm > 0);

			if(!pingpongm)
			{
				load_compute_wrapper(ifm, weight, ofm_buffer1, ksize, kstride, ifm_num, ifm_w, ifm_h, ofm_num,
		 			pad_int, ltype, TRow, TCol, IHW, KK, INumxKK, TC_MIN, TR_MIN, TM_MIN, TM, TN, tm, tr, tc, 
					m1, TM_MIN1, pingpongm, input_flag, process_flag, pad_val, LoadBias);

				write_back_output_reorg( ofm_buffer0, bias_buffer, ofm, tr, tc, m0[0], ofm_w, TM_MIN0[0], TR_MIN, TC_MIN, OHW, IsReLU, LoadBias, write_flag);
				pingpongm = true;
			}else
			{
				load_compute_wrapper(ifm, weight, ofm_buffer0, ksize, kstride, ifm_num, ifm_w, ifm_h, ofm_num,
		 			pad_int, ltype, TRow, TCol, IHW, KK, INumxKK, TC_MIN, TR_MIN, TM_MIN, TM, TN, tm, tr, tc,
					m0, TM_MIN0, pingpongm, input_flag, process_flag, pad_val, LoadBias);

				write_back_output_reorg( ofm_buffer1, bias_buffer, ofm, tr, tc, m1[0], ofm_w, TM_MIN1[0], TR_MIN, TC_MIN, OHW, IsReLU, LoadBias, write_flag);
				pingpongm = false;
			}
		    }	
		}
	}
}

