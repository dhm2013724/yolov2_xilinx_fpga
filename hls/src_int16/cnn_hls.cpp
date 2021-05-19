
#include "cnn.h"

void ifm_mmcpy_row(int32_t *input, int32_t local_buf[OnChipIB_Width/2+3], int CurrentOffset, uint32_t IHxIW, int IW_align_256b, uint16_t TCol,
	uint8_t t1, uint8_t t2, uint8_t *t1_n, uint8_t *t2_n, uint8_t *bn_n, bool enable)
{
	if(!enable)
		return;

	int ifm_offset = CurrentOffset + t1*IHxIW + t2*IW_align_256b;
	int ifm_trans_offset = (ifm_offset >> 1);
	uint8_t begin_num = ifm_offset & 0x1;
	uint16_t TCol_a = TCol + begin_num;
	uint16_t loop_cnts = TCol_a >> 1;
	if(TCol_a & 0x1)
		loop_cnts++;

	assert((loop_cnts > 0)&&(loop_cnts <= (OnChipIB_Width/2+3)));
	memcpy(local_buf, (int32_t *)(input + ifm_trans_offset), loop_cnts*sizeof(int32_t));

	*t1_n = t1;
	*t2_n = t2;
	*bn_n = begin_num;
}

void ifm_copy_lbuf2ibuf(int16_t input_buffer[Tn][OnChipIB_Height][OnChipIB_Width], int32_t local_buf[OnChipIB_Width/2+3], uint16_t TCol, uint16_t Input_w, uint16_t Input_h, uint8_t TN_MIN,
	int16_t pad_val, int Coffset, int Roffset, uint8_t t1, uint8_t t2, uint8_t bn, bool enable)
{
	if(!enable)
		return;

	bool TN_Enable = t1 < TN_MIN;
	int yoffset = Roffset + t2;
	bool YEnable = (yoffset >= 0)&&(yoffset < Input_h);
	bool PEnable = YEnable&&TN_Enable;

	uint16_t cnt = 1;
	uint8_t bn_local = bn;
	int16_t buf_256b[2];
#pragma HLS ARRAY_PARTITION variable=buf_256b complete dim=1
	buf_256b[0] = local_buf[0] & 0xFFFF;
	buf_256b[1] = (local_buf[0] >> 16) & 0xFFFF;
	for(uint8_t t3 = 0;t3 < TCol; t3++)
	{
//DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=TCol_max)
#pragma HLS PIPELINE II=1
		int xoffset = Coffset + t3;
		bool XEnable = (xoffset >= 0)&&(xoffset < Input_w);
		if(XEnable&&PEnable)
		{
			input_buffer[t1][t2][t3] = buf_256b[bn_local];
		}
		else
			input_buffer[t1][t2][t3] = pad_val;
		bn_local++;
		if(bn_local==2)
		{
			bn_local = 0;
			buf_256b[0] = local_buf[cnt] & 0xFFFF;
			buf_256b[1] = (local_buf[cnt] >> 16) & 0xFFFF;
			cnt++;
		}
	}

	assert((cnt>=0)&&(cnt<=(OnChipIB_Width/2+3)));
}

void input_load(int32_t *input, int16_t input_buffer[Tn][OnChipIB_Height][OnChipIB_Width], uint16_t r, uint16_t c, uint16_t n, uint8_t Kstride, uint8_t Padding,
		uint16_t TRow, uint16_t TCol, uint16_t Input_w, int IW_align_256b, uint16_t Input_h, uint8_t TN_MIN, uint32_t IHxIW, int16_t pad_val)
{
	uint8_t t1,t2;
	static int32_t local_buf0[OnChipIB_Width/2+3];
	static int32_t local_buf1[OnChipIB_Width/2+3];

	const int Coffset = c*Kstride - Padding;
	const int Roffset = r*Kstride - Padding;
	const int CurrentOffset = n*IHxIW + Roffset*IW_align_256b + Coffset;

	uint8_t t1_n0, t1_n1, t2_n0, t2_n1;
	uint8_t bn_n0, bn_n1;
	bool pp = true;

	int TnxTRow = Tn*TRow;
	int t = 0;
	t1 = 0; t2 = 0;
	for(t = 0;t < TnxTRow+1; t++)
	{
//DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=TRow_max+1)
		if(pp)
		{
			ifm_mmcpy_row(input, local_buf0, CurrentOffset, IHxIW, IW_align_256b, TCol, t1, t2, &t1_n0, &t2_n0, &bn_n0, t!=TnxTRow);
			ifm_copy_lbuf2ibuf( input_buffer, local_buf1, TCol, Input_w, Input_h, TN_MIN, pad_val, Coffset, Roffset, t1_n1, t2_n1, bn_n1, t!=0);
			pp = false;
		}else
		{
			ifm_mmcpy_row(input, local_buf1, CurrentOffset, IHxIW, IW_align_256b, TCol, t1, t2, &t1_n1, &t2_n1, &bn_n1, t!=TnxTRow);
			ifm_copy_lbuf2ibuf( input_buffer, local_buf0, TCol, Input_w, Input_h, TN_MIN, pad_val, Coffset, Roffset, t1_n0, t2_n0, bn_n0, t!=0);
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

void weight_load_reorg(int32_t *Weight,int16_t weight_buffer[Tm][Tn][K][K],int m,int n, uint32_t IFM_numxKxK, uint8_t KxK, uint8_t Ksize, uint8_t TM_MIN, uint8_t TN_MIN, bool enable)
{
	if(!enable)
		return;

	uint8_t t1,t2,t3,t4;
	static int32_t local_buf[(Tm*Tn*K*K+1)/2+1];
	static int Woffset;

	assert((TM_MIN > 0)&&(TM_MIN <= Tm));
	assert((TN_MIN > 0)&&(TN_MIN <= Tn));
	assert((KxK > 0)&&(KxK <= K*K));

	if(m==0&&n==0)
		Woffset = 0;

	uint16_t mm_offset = (TM_MIN*TN_MIN*KxK+1) >> 1;
	assert((mm_offset > 0)&&(mm_offset <= (Tm*Tn*K*K+1)/2));
	memcpy(local_buf, (int32_t *)(Weight + Woffset), mm_offset*sizeof(int32_t));
	Woffset += mm_offset;

	uint16_t lb_cnt = 1;
	int16_t buf_256b[2];
	uint8_t cnt = 0;
	buf_256b[0] = local_buf[0] & 0xFFFF;
	buf_256b[1] = (local_buf[0] >> 16) & 0xFFFF;

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
						weight_buffer[t1][t2][t3][t4] =  buf_256b[cnt];
						cnt++;
						if(cnt==2)
						{
							cnt = 0;
							buf_256b[0] = local_buf[lb_cnt] & 0xFFFF;
							buf_256b[1] = (local_buf[lb_cnt] >> 16) & 0xFFFF;
							lb_cnt++;
						}
					}
					else
						weight_buffer[t1][t2][t3][t4] = 0;
				}
}

void nonlinear_leaky_row(int32_t local_buf[Tc/2], int32_t Input[Tm][Tr][Tc], uint8_t tm, uint8_t tr, uint8_t *tm_n, uint8_t *tr_n, uint8_t TC_MIN,
				bool IsNL, bool enable, int InterSubOutput, int ltype)
{
	if(!enable)
		return ;

	uint8_t tc;
//	float tmp_out, tmp, tmp1;
	assert((TC_MIN>0)&&(TC_MIN<=Tc));
	assert((InterSubOutput>0)&&(InterSubOutput<32));

	uint8_t cnt = 0;
	uint8_t bn_local = 0;//ap_uint<4> 0-15

	int32_t tmp_out, tmp1;
	int16_t tmp_int16;

	int16_t buf_256b[2];
	for(tc = 0;tc < TC_MIN;tc++)
	{
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tc)
#pragma HLS PIPELINE II=1
		tmp1 = Input[tm][tr][tc];

		if(IsNL)
		{
			if(tmp1 < 0)
				tmp_out = ((int64_t)tmp1*0xccc)>>15;
			else
				tmp_out = tmp1;
		}
		else
			tmp_out = tmp1;

		if(ltype==LT_CONV)
		{
			tmp_int16 = tmp_out >> InterSubOutput;//tmp_out*pow(2.0, OutputQ-Inter)
		}
		else
		{
			tmp_int16 = tmp_out;
		}
		buf_256b[bn_local] = tmp_int16;

		bn_local++;
		if(bn_local == 2)
		{
			bn_local = 0;
			local_buf[cnt] = ((buf_256b[1] << 16) & 0xFFFF0000)| (buf_256b[0] & 0x0000FFFF);
			cnt++;
		}
	}

	if(bn_local & 0xF)
	{
		local_buf[cnt] = ((buf_256b[1] << 16) & 0xFFFF0000)| (buf_256b[0] & 0x0000FFFF);
	}

	*tm_n = tm;
	*tr_n = tr;
}

void ofm_mmcpy_row(int32_t *Output, int32_t local_buf[Tc/2], int offset, uint32_t OHxOW, uint16_t Output_w, uint8_t TC_MIN, uint8_t tm, uint8_t tr,bool enable)
{
	if(!enable)
		return;

	int ofm_offset = tm*OHxOW + tr*Output_w + offset;
	int trans_offset = (ofm_offset >> 1);
	uint16_t loop_cnts = TC_MIN >> 1;
	if(TC_MIN & 0x1)
		loop_cnts++;

	assert((loop_cnts > 0)&&(loop_cnts <=(Tc/2)));

	memcpy((int32_t *)(Output + trans_offset), local_buf, loop_cnts*sizeof(int32_t));

}

void write_back_output_reorg(int32_t output_buffer[Tm][Tr][Tc], int32_t *Output, uint16_t r,uint16_t c,uint16_t m, uint16_t Output_w,
		uint8_t TM_MIN, uint8_t TR_MIN, uint8_t TC_MIN, uint32_t OHxOW, bool IsNL, bool enable, int InterSubOutput, int ltype)
{
	if(!enable)
		return;

	assert((TM_MIN >0)&&(TM_MIN <=Tm));
	assert((TR_MIN >0)&&(TR_MIN <=Tr));
	assert((TC_MIN >0)&&(TC_MIN <=Tc));

	const int offset = m*OHxOW + r*Output_w + c;
	static int32_t local_buf0[Tc/2];
	static int32_t local_buf1[Tc/2];
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
			nonlinear_leaky_row( local_buf0, output_buffer, tm, tr, &tm_n0, &tr_n0, TC_MIN, IsNL, t!=TM_MINxTR_MIN, InterSubOutput, ltype);
			ofm_mmcpy_row( Output, local_buf1, offset, OHxOW, Output_w, TC_MIN, tm_n1, tr_n1, t!=0);
			pp = false;
		}else
		{
			nonlinear_leaky_row( local_buf1, output_buffer, tm, tr, &tm_n1, &tr_n1, TC_MIN, IsNL, t!=TM_MINxTR_MIN, InterSubOutput, ltype);
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

}

/*
void maxpool_tile(int16_t Input[Tn][OnChipIB_Height][OnChipIB_Width],int32_t Output[Tm][Tr][Tc],
		  uint8_t Ksize, uint8_t Kstride, uint8_t TM_MIN, uint8_t TR_MIN, uint8_t TC_MIN, bool enable)
{
	if(!enable)
		return;

	uint8_t i,j,tr,tc,of;
	int16_t tmp[Tn];
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
							tmp[of] = Output[of][tr][tc] & 0xFFFF;

						int16_t tmp_in = Input[of][tr*Kstride+i][tc*Kstride+j];

						if(tmp_in > tmp[of])
							tmp[of] = tmp_in;

						Output[of][tr][tc] = tmp[of];
					}
				}

}*/

void maxpool_tile(int16_t Input[Tn][OnChipIB_Height][OnChipIB_Width],int32_t Output[Tm][Tr][Tc],
		  uint8_t Ksize, uint8_t K_1, uint8_t Kstride, uint8_t TM_MIN, uint8_t TR_MIN, uint8_t TC_MIN, bool enable)
{
	if(!enable)
		return;

	uint8_t i,j,tr,tc,of;
	int16_t tmp[Tn];
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
							tmp[of] = MIN_NEG;

						int16_t tmp_in = Input[of][tr*Kstride+i][tc*Kstride+j];

						if(tmp_in > tmp[of])
							tmp[of] = tmp_in;

						if(i==K_1&&j==K_1)
							Output[of][tr][tc] = tmp[of];
					}
				}

}

void copy_local_beta(int16_t beta_buffer[MAX_BETA_LENGTH],int32_t local_beta_buffer[Tm],const int TM_MIN,int m, int InterSubBeta)
{
	assert((InterSubBeta > 0)&&(InterSubBeta <32));

	int offset = m;
	int tm;
	for(tm = 0;tm < TM_MIN;tm++)
	{
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tm)
#pragma HLS PIPELINE II=1
		local_beta_buffer[tm] = beta_buffer[offset] << InterSubBeta;
		offset++;
	}
}

void conv2d_tile(int16_t input_buffer[Tn][OnChipIB_Height][OnChipIB_Width],int32_t output_buffer[Tm][Tr][Tc],
		int16_t weight_buffer[Tm][Tn][K][K],int16_t beta_buffer[MAX_BETA_LENGTH],int n_next,
		const int Ksize,const int Kstride,int m,
		const int TM_MIN,const int TR_MIN,const int TC_MIN,bool enable, int InterSubBeta, int WeightAddInputSubInter, int InterSubOutput)
{
	static int32_t local_beta_buffer[Tm];
#pragma HLS ARRAY_PARTITION variable=local_beta_buffer complete dim=1

	if(!enable)
	{
		copy_local_beta(beta_buffer,local_beta_buffer,TM_MIN, m, InterSubBeta);
		return;
	}

	assert((WeightAddInputSubInter >= 0)&&(WeightAddInputSubInter < 32));

	ap_uint<3> i,j;
	ap_uint<5> WeightAddInputSubInter_5b = WeightAddInputSubInter;
	ap_uint<3> Ksize_3b = Ksize;
	ap_uint<3> Kstride_3b = Kstride;
	ap_uint<8> TR_MIN_8b = TR_MIN;
	ap_uint<8> TC_MIN_8b = TC_MIN;
	ap_uint<8> tr,tc,tm,tn;
	const bool ne0 = (n_next==0);

	int32_t partial_mul[Tm][Tn];
#pragma HLS ARRAY_PARTITION variable=partial_mul complete dim=1
#pragma HLS ARRAY_PARTITION variable=partial_mul complete dim=2
	int32_t partial_add[Tm];
#pragma HLS ARRAY_PARTITION variable=partial_add complete dim=1

	for(i =0;i < Ksize_3b; i++)
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=K)
		for(j = 0;j < Ksize_3b; j++)
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=K)
			for(tr = 0;tr < TR_MIN_8b;tr++)
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tr)
				for(tc = 0;tc < TC_MIN_8b;tc++)
				{
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tc)
#pragma HLS PIPELINE II=1
					for(tm = 0;tm < Tm;tm++)
					{
#pragma HLS DEPENDENCE variable=output_buffer inter false

						if(i==0&&j==0&&ne0)
							partial_add[tm] = local_beta_buffer[tm];
						else
							partial_add[tm] = output_buffer[tm][tr][tc];

						for(tn = 0;tn <Tn;tn++)
						{
							partial_mul[tm][tn] = weight_buffer[tm][tn][i][j]*input_buffer[tn][Kstride_3b*tr+i][Kstride_3b*tc+j] >> WeightAddInputSubInter_5b;
						}

						int32_t partial_sum = 0;
						for(tn = 0;tn <Tn;tn++)
						{
							 partial_sum += partial_mul[tm][tn];
						}
						output_buffer[tm][tr][tc] = (partial_add[tm] + partial_sum);
					}
				}
}

void ifm_weight_load_wrapper(int32_t *ifm, int32_t *weight, int16_t ifm_buffer[Tn][OnChipIB_Height][OnChipIB_Width], int16_t weight_buffer[Tm][Tn][K][K],
	int ifm_num, int ofm_num, int ifm_w, int IW_align_256b, int ifm_h, int pad_int,
	int TM, int TN, int tr, int tc, int tm, int tn, int tn_next[1], int ksize, int kstride, int ltype, int TM_MIN, int TRow, int TCol, int IHW, int INumxKK, int KK,
	bool enable, int16_t pad_val, bool LoadBias)
{
	if(!enable)
		return;

	tn_next[0] = tn;

	int TN_MIN, tnm, TNM_MIN;
	if(ltype == LT_CONV)
	{
		TN_MIN = MIN_diy(TN, ifm_num-tn);
		tnm = tn;
		TNM_MIN = TN_MIN;
	}
	else
	{
		TN_MIN = 1;
		tnm = tm;
		TNM_MIN = TM_MIN;
	}

	input_load(ifm, ifm_buffer, tr, tc, tnm, kstride, pad_int, TRow, TCol, ifm_w, IW_align_256b, ifm_h, TNM_MIN, IHW, pad_val);
	weight_load_reorg(weight,weight_buffer,tm,tn,INumxKK,KK,ksize,TM_MIN,TN_MIN, LoadBias);
}

void load_compute_wrapper(int32_t *ifm, int32_t *weight, int32_t ofm_buffer[Tm][Tr][Tc], int16_t bias_buffer[MAX_BETA_LENGTH], int ksize, uint8_t K_1, int kstride, int ifm_num, int ifm_w, int IW_align_256b,
	 int ifm_h,int ofm_num, int pad_int, int ltype, int TRow, int TCol, int IHW, int KK, int INumxKK, int TC_MIN, int TR_MIN, int TM_MIN, int TM, int TN, int tm, int tr, int tc,
	 int TMP_X_next[1],int TX_MIN_next[1],bool pingpongx,bool input_flag,bool process_flag, int16_t pad_val, bool LoadBias, int InterSubBeta, int WeightAddInputSubInter, int InterSubOutput)
{
	static int16_t ifm_buffer0[Tn][OnChipIB_Height][OnChipIB_Width];
#pragma HLS ARRAY_PARTITION variable=ifm_buffer0 complete dim=1
	static int16_t ifm_buffer1[Tn][OnChipIB_Height][OnChipIB_Width];
#pragma HLS ARRAY_PARTITION variable=ifm_buffer1 complete dim=1

	static int16_t weight_buffer0[Tm][Tn][K][K];
#pragma HLS ARRAY_PARTITION variable=weight_buffer0 complete dim=1
#pragma HLS ARRAY_PARTITION variable=weight_buffer0 complete dim=2
	static int16_t weight_buffer1[Tm][Tn][K][K];
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
				ifm_weight_load_wrapper(ifm, weight, ifm_buffer1, weight_buffer1, ifm_num, ofm_num, ifm_w, IW_align_256b, ifm_h,
					pad_int, TM, TN, tr, tc, tm, tn, n1, ksize, kstride, ltype, TM_MIN, TRow, TCol, IHW, INumxKK, KK, tn < ifm_num, pad_val, LoadBias);
				conv2d_tile(ifm_buffer0, ofm_buffer, weight_buffer0, bias_buffer, n0[0], ksize, kstride, tm, TM_MIN, TR_MIN, TC_MIN, tn!=0,
 					InterSubBeta,WeightAddInputSubInter,InterSubOutput);
				pingpong = false;
			}else
			{
				ifm_weight_load_wrapper(ifm, weight, ifm_buffer0, weight_buffer0, ifm_num, ofm_num, ifm_w, IW_align_256b, ifm_h,
					pad_int, TM, TN, tr, tc, tm, tn, n0, ksize, kstride, ltype, TM_MIN, TRow, TCol, IHW, INumxKK, KK, tn < ifm_num, pad_val, LoadBias);
				conv2d_tile(ifm_buffer1, ofm_buffer, weight_buffer1, bias_buffer, n1[0], ksize, kstride, tm, TM_MIN, TR_MIN, TC_MIN, tn!=0,
					InterSubBeta,WeightAddInputSubInter,InterSubOutput);
				pingpong = true;
			}
		}
	}else if(ltype == LT_MAXPOOL)
	{
		if(pingpongx)
		{
			ifm_weight_load_wrapper(ifm, weight, ifm_buffer0, weight_buffer0, ifm_num, ofm_num, ifm_w, IW_align_256b, ifm_h,
				pad_int, TM, TN, tr, tc, tm, 0, n0, ksize, kstride, ltype, TM_MIN, TRow, TCol, IHW, INumxKK, KK, input_flag, pad_val, LoadBias);
			maxpool_tile(ifm_buffer1, ofm_buffer, ksize, K_1, kstride, tmp_tx_min, TR_MIN, TC_MIN, process_flag);

			TMP_X_next[0] = tmp_x;
			TX_MIN_next[0] = tmp_tx_min;
			tmp_x = tm;
			tmp_tx_min = TM_MIN;
		}else
		{
			ifm_weight_load_wrapper(ifm, weight, ifm_buffer1, weight_buffer1, ifm_num, ofm_num, ifm_w, IW_align_256b, ifm_h,
				pad_int, TM, TN, tr, tc, tm, 0, n1, ksize, kstride, ltype, TM_MIN, TRow, TCol, IHW, INumxKK, KK, input_flag, pad_val, LoadBias);
			maxpool_tile(ifm_buffer0, ofm_buffer, ksize, K_1, kstride, tmp_tx_min, TR_MIN, TC_MIN, process_flag);

			TMP_X_next[0] = tmp_x;
			TX_MIN_next[0] = tmp_tx_min;
			tmp_x = tm;
			tmp_tx_min = TM_MIN;
		}
	}
}

void copy_beta(int16_t beta_buffer[MAX_BETA_LENGTH], int32_t *Beta, uint16_t OFM_NUM)
{
	static int32_t local_buf[(MAX_BETA_LENGTH+1)/2+1];
	uint16_t NUM = (OFM_NUM+1) >> 1;

	assert((NUM > 0)&&(NUM <= (MAX_BETA_LENGTH+1)/2));

	memcpy(local_buf, Beta, NUM*sizeof(int32_t));

	uint16_t cnt = 1;
	uint8_t bn = 0;
	int16_t tmp_256b[2];
#pragma HLS ARRAY_PARTITION variable=tmp_256b complete dim=1
	tmp_256b[0] = local_buf[0] & 0xFFFF;
	tmp_256b[1] = (local_buf[0] >> 16) & 0xFFFF;
	for(int t = 0; t < OFM_NUM; t++)
	{
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=MAX_BETA_LENGTH)
#pragma HLS PIPELINE II=1
		int16_t tmp_i16 = tmp_256b[bn];
		beta_buffer[t] = tmp_i16;
		bn++;
		if(bn==2)
		{
			bn = 0;
			tmp_256b[0] = local_buf[cnt] & 0xFFFF;
			tmp_256b[1] = (local_buf[cnt] >> 16) & 0xFFFF;
			cnt++;
		}
	}
}

void FPGA_Acc(int32_t *ifm, int32_t *ofm, int32_t *weight, int32_t *bias, uint32_t k_s_pad_ltype, uint32_t iofm_num, uint32_t ifm_w_h, uint32_t ofm_w_h,
	uint32_t TRTC, uint32_t TMTN, int OFM_num_bound, int mLoopsxTM, int mLoops_a1xTM, int16_t pad_val, uint32_t TRowTCol,
	uint32_t IHW, uint32_t OHW, uint32_t KK_INumxKK, uint32_t en_bits, int WeightQ, int BetaQ, int InputQ, int OutputQ)//enable_bits[2:0]={IsReLU, LoadBias, IsNotConv}
{
#pragma HLS INTERFACE m_axi depth=512 port=ifm    offset=slave bundle=DATA_BUS num_read_outstanding=1 num_write_outstanding=1 max_read_burst_length=64 max_write_burst_length=64
#pragma HLS INTERFACE m_axi depth=512 port=ofm offset=slave bundle=DATA_BUS num_read_outstanding=1 num_write_outstanding=1 max_read_burst_length=64 max_write_burst_length=64
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

#pragma HLS INTERFACE s_axilite register port=WeightQ bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=BetaQ bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=InputQ bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=OutputQ bundle=CTRL_BUS

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

	assert((ofm_num > 0)&&(ofm_num <= 2048));
	assert((ifm_num > 0)&&(ifm_num <= 2048));
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

///////////////////////////////////
	const int InterSubBeta = INTERWIDTH - BetaQ;
	const int WeightAddInputSubInter = WeightQ + InputQ - INTERWIDTH;
	const int InterSubOutput = INTERWIDTH - OutputQ;

	assert((InterSubBeta >= 0)&&(InterSubBeta < 32));
	assert((WeightAddInputSubInter >= 0)&&(WeightAddInputSubInter < 32));
	assert((InterSubOutput >= 0)&&(InterSubOutput < 32));

	static int32_t ofm_buffer0[Tm][Tr][Tc];
#pragma HLS ARRAY_PARTITION variable=ofm_buffer0 complete dim=1
	static int32_t ofm_buffer1[Tm][Tr][Tc];
#pragma HLS ARRAY_PARTITION variable=ofm_buffer1 complete dim=1
	static int16_t bias_buffer[MAX_BETA_LENGTH];

	uint16_t IW_align_256b = (ifm_w >> 1) << 1;
	if(ifm_w & 0x1)
		IW_align_256b += 2;

	uint16_t OW_align_256b = (ofm_w >> 1) << 1;
	if(ofm_w & 0x1)
		OW_align_256b += 2;

	assert((OW_align_256b%2==0)&&(IW_align_256b%2==0));
	const int OHxOW = ofm_h*OW_align_256b;
	const int IHxIW = ifm_h*IW_align_256b;

	uint8_t K_1 = ksize - 1;
/////////////////////////////////param

	if(LoadBias)
		copy_beta(bias_buffer ,bias, ofm_num);
		//memcpy(bias_buffer, bias, ofm_num*sizeof(float));

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
				load_compute_wrapper(ifm, weight, ofm_buffer1, bias_buffer, ksize, K_1, kstride, ifm_num, ifm_w, IW_align_256b, ifm_h, ofm_num,
		 			pad_int, ltype, TRow, TCol, IHxIW, KK, INumxKK, TC_MIN, TR_MIN, TM_MIN, TM, TN, tm, tr, tc,
					m1, TM_MIN1, pingpongm, input_flag, process_flag, pad_val, LoadBias,InterSubBeta,WeightAddInputSubInter,InterSubOutput);

				write_back_output_reorg( ofm_buffer0, ofm, tr, tc, m0[0], OW_align_256b, TM_MIN0[0], TR_MIN, TC_MIN, OHxOW, IsReLU,
					write_flag, InterSubOutput, ltype);
				pingpongm = true;
			}else
			{
				load_compute_wrapper(ifm, weight, ofm_buffer0, bias_buffer, ksize, K_1, kstride, ifm_num, ifm_w, IW_align_256b, ifm_h, ofm_num,
		 			pad_int, ltype, TRow, TCol, IHxIW, KK, INumxKK, TC_MIN, TR_MIN, TM_MIN, TM, TN, tm, tr, tc,
					m0, TM_MIN0, pingpongm, input_flag, process_flag, pad_val, LoadBias,InterSubBeta,WeightAddInputSubInter,InterSubOutput);

				write_back_output_reorg( ofm_buffer1, ofm, tr, tc, m1[0], OW_align_256b, TM_MIN1[0], TR_MIN, TC_MIN, OHxOW, IsReLU,
					 write_flag, InterSubOutput, ltype);
				pingpongm = false;
			}
		    }
		}
	}
}
