
#include "cnn.h"

void ifm_mmcpy_row(ap_uint<256> *input, ap_uint<256> local_buf[OnChipIB_Width/8+3], int CurrentOffset, int IHxIW, int IW_align_256b, int TCol,
	uint8_t t1, uint8_t t2, uint8_t *t1_n, uint8_t *t2_n, uint8_t *bn_n, bool enable)
{
	if(!enable)
		return;

	int ifm_offset = CurrentOffset + t1*IHxIW + t2*IW_align_256b;
//	int ifm_trans_offset = ((ifm_offset >>3) << 3);
	int ifm_trans_offset = (ifm_offset >>3);
	uint8_t begin_num = ifm_offset & 0x7;
	uint16_t TCol_a = TCol + begin_num;
	uint16_t loop_cnts = TCol_a >> 3;
	if(TCol_a & 0x7)
		loop_cnts++;
	for(int t = 0; t < loop_cnts; t++)
	{
#pragma HLS PIPELINE II=1
		local_buf[t].range(255, 0) = input[ifm_trans_offset + t].range(255, 0);
//		memcpy( local_buf[t],(float *)(input + ifm_trans_offset + t*8), 8*sizeof(float));
	}

	*t1_n = t1;
	*t2_n = t2;
	*bn_n = begin_num;
}

void ifm_copy_lbuf2ibuf(float input_buffer[Tn][OnChipIB_Height][OnChipIB_Width], ap_uint<256> local_buf[OnChipIB_Width/8+3], int TCol, int Input_w, int Input_h, int TN_MIN, float pad_value,
	int Coffset, int Roffset, uint8_t t1, uint8_t t2, uint8_t bn, bool enable)
{
	if(!enable)
		return;

	bool TN_Enable = t1 < TN_MIN;
	int yoffset = Roffset + t2;
	bool YEnable = (yoffset >= 0)&&(yoffset < Input_h);
	bool PEnable = YEnable&&TN_Enable;

	uint16_t cnt = 1;
	uint8_t bn_local = bn;
	ap_uint<256> buf_256b;
	buf_256b.range(255, 0) =  local_buf[0].range(255, 0);
//	memcpy(buf_256b, local_buf[0], 8*sizeof(float));
	for(uint8_t t3 = 0;t3 < TCol; t3++)
	{
//DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=TCol_max)
#pragma HLS PIPELINE II=1
		int xoffset = Coffset + t3;
		bool XEnable = (xoffset >= 0)&&(xoffset < Input_w);
		if(XEnable&&PEnable)
		{
			union {float f;uint32_t u;} tmp_32;
			tmp_32.u = buf_256b.range(bn_local*32+31, bn_local*32);
			input_buffer[t1][t2][t3] = tmp_32.f;
		}
		else
			input_buffer[t1][t2][t3] = pad_value;
		bn_local++;
		if(bn_local==8)
		{
			bn_local = 0;
			buf_256b.range(255, 0) = local_buf[cnt].range(255, 0);
//			memcpy(buf_256b, local_buf[cnt], 8*sizeof(float));
			cnt++;
		}
	}
}

void input_load(ap_uint<256> *input,float input_buffer[Tn][OnChipIB_Height][OnChipIB_Width],int r,int c,int n,int Kstride,int Padding,int TRow,int TCol,int Input_w, int IW_align_256b,int Input_h,int TN_MIN,int IHxIW,int LayerType)
{
	uint8_t t1,t2;
	static ap_uint<256> local_buf0[OnChipIB_Width/8+3];
	static ap_uint<256> local_buf1[OnChipIB_Width/8+3];

	const int Coffset = c*Kstride - Padding;
	const int Roffset = r*Kstride - Padding;
	const int CurrentOffset = n*IHxIW + Roffset*IW_align_256b + Coffset;

	uint8_t t1_n0, t1_n1, t2_n0, t2_n1;
	uint8_t bn_n0, bn_n1;
	bool pp = true;

	float pad_value = 0.0f;
	if(LayerType==1)
		pad_value = -1024*1024;

	int TnxTRow = Tn*TRow;
	int t = 0;
	t1 = 0; t2 = 0;
	for(t = 0;t < TnxTRow+1; t++)
	{
//DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=TRow_max)
		if(pp)
		{
			ifm_mmcpy_row(input, local_buf0, CurrentOffset, IHxIW, IW_align_256b, TCol, t1, t2, &t1_n0, &t2_n0, &bn_n0, t!=TnxTRow);
			ifm_copy_lbuf2ibuf( input_buffer, local_buf1, TCol, Input_w, Input_h, TN_MIN, pad_value, Coffset, Roffset, t1_n1, t2_n1, bn_n1, t!=0);
			pp = false;
		}else
		{
			ifm_mmcpy_row(input, local_buf1, CurrentOffset, IHxIW, IW_align_256b, TCol, t1, t2, &t1_n1, &t2_n1, &bn_n1, t!=TnxTRow);
			ifm_copy_lbuf2ibuf( input_buffer, local_buf0, TCol, Input_w, Input_h, TN_MIN, pad_value, Coffset, Roffset, t1_n0, t2_n0, bn_n0, t!=0);
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

void weight_load_reorg(ap_uint<256> *Weight,float weight_buffer[Tm][Tn][K][K],bool weight_load_enable,int m,int n,int IFM_numxKxK,int KxK,int Ksize,int TM_MIN,int TN_MIN)
{
	uint8_t t1,t2,t3,t4;
	static ap_uint<256> local_buf[(Tm*Tn*K*K)/8 + 3];
	static int Woffset;

	assert((TM_MIN > 0)&&(TM_MIN <= Tm));
	assert((TN_MIN > 0)&&(TN_MIN <= Tn));
	assert((KxK > 0)&&(KxK <= K*K));

	if(!weight_load_enable)
		return;

	if(m==0&&n==0)
		Woffset = 0;

	uint16_t mm_offset = TM_MIN*TN_MIN*KxK;

//	uint32_t trans_offset = ((Woffset >>3) << 3);
	uint32_t trans_offset = (Woffset >>3);
	uint8_t begin_num = Woffset & 0x7;//ap_uint<3>
	uint16_t TCol_a = mm_offset + begin_num;
	uint16_t loop_cnts = TCol_a >> 3;
	if(TCol_a & 0x7)
		loop_cnts++;
	for(int t = 0; t < loop_cnts; t++)
	{
#pragma HLS PIPELINE II=1
		local_buf[t].range(255, 0) = Weight[trans_offset + t].range(255, 0);
//		memcpy( local_buf[t],(float *)(Weight + trans_offset + t*8), 8*sizeof(float));
	}
	Woffset += mm_offset;

	uint16_t cnt = 1;
	uint8_t bn_local = begin_num;
	ap_uint<256> buf_256b;
	buf_256b.range(255, 0) = local_buf[0].range(255, 0);
//	memcpy(buf_256b, local_buf[0], 8*sizeof(float));

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
						union {float f; uint32_t u;} tmp_32;
						tmp_32.u = buf_256b.range(bn_local*32+31,bn_local*32);
						weight_buffer[t1][t2][t3][t4] =  tmp_32.f;
						bn_local++;
						//bn_local = bn_local % 8;//ap_uint<3> 0-7
						if(bn_local==8)
						{
							bn_local = 0;
							buf_256b.range(255, 0) = local_buf[cnt].range(255, 0);
//							memcpy(buf_256b, local_buf[cnt], 8*sizeof(float));
							cnt++;
						}
					}
					else
						weight_buffer[t1][t2][t3][t4] = 0;
				}
}


void copy_input_weight(ap_uint<256> *input,ap_uint<256> *Weight,int IFM_num,int Input_w, int IW_align_256b, int Input_h, int Ksize,int Kstride,int r,int c,int m,int n,
		int TM_MIN,int TN,int TRow,int TCol,int Padding,float input_buffer[Tn][OnChipIB_Height][OnChipIB_Width],float weight_buffer[Tm][Tn][K][K],int n_next[1],
		bool enable,bool weight_load_enable,bool initialize,const int IHxIW,const int KxK,const int IFM_numxKxK,const int LayerType)
{
	if(!enable)
		return ;

	const int TN_MIN = MIN(TN, IFM_num-n);
	n_next[0] = n;

	input_load(input, input_buffer, r, c, n, Kstride, Padding, TRow, TCol, Input_w, IW_align_256b, Input_h, TN_MIN, IHxIW, LayerType);

	weight_load_reorg(Weight,weight_buffer,weight_load_enable,m,n,IFM_numxKxK,KxK,Ksize,TM_MIN,TN_MIN);
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

void nonlinear_leaky_row(ap_uint<256> local_buf[Tc/8+2], float Input[Tm][Tr][Tc], uint8_t tm, uint8_t tr, uint8_t *tm_n, uint8_t *tr_n, uint8_t TC_MIN,const bool IsNL, bool enable)
{
	if(!enable)
		return ;

	uint8_t tc;
	assert((TC_MIN>0)&&(TC_MIN<=Tc));

	uint8_t cnt = 0;
	uint8_t bn_local = 0;//ap_uint<3> 0-7

	float tmp_out;
	for(tc = 0;tc < TC_MIN;tc++)
	{
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tc)
#pragma HLS PIPELINE II=1
		float tmp = Input[tm][tr][tc];
		if((tmp < 0.0f)&&IsNL)
			tmp_out = tmp*0.1f;
		else
			tmp_out = tmp;
		union {float f; uint32_t u;} tmp_32;
		tmp_32.f = tmp_out;
		local_buf[cnt].range(bn_local*32+31, bn_local*32) = tmp_32.u;
		bn_local++;
		//bn_local = bn_local % 8;//ap_uint<3> 0-7
		if(bn_local == 8)
		{
			bn_local = 0;
			cnt++;
		}
	}

	*tm_n = tm;
	*tr_n = tr;
}

void ofm_mmcpy_row(ap_uint<256> *Output, ap_uint<256> local_buf[Tc/8+2], int offset, int OHxOW, int Output_w, int TC_MIN, uint8_t tm, uint8_t tr,bool enable)
{
	if(!enable)
		return;

	int ofm_offset = tm*OHxOW + tr*Output_w + offset;
//	int trans_offset = (ofm_offset >> 3) << 3;
	int trans_offset = (ofm_offset >> 3);
	uint16_t loop_cnts = TC_MIN >> 3;
	if(TC_MIN & 0x7)
		loop_cnts++;
	for(int t = 0;t < loop_cnts; t++)
	{
#pragma HLS PIPELINE II=1
		Output[trans_offset + t].range(255, 0) = local_buf[t].range(255, 0);
//		memcpy((float *)(Output + trans_offset + t*8), local_buf[t], 8*sizeof(float));
	}

}

void write_back_output_reorg(float output_buffer[Tm][Tr][Tc], ap_uint<256> *Output,int r,int c,int m,uint16_t Output_w,uint16_t Output_h,
		uint8_t TM_MIN,uint8_t TR_MIN,uint8_t TC_MIN,const int OHxOW, bool IsNL, bool write_flag)
{
	if(!write_flag)
		return;

	assert((TM_MIN >0)&&(TM_MIN <=Tm));
	assert((TR_MIN >0)&&(TR_MIN <=Tr));
	assert((TC_MIN >0)&&(TC_MIN <=Tc));

	const int offset = m*OHxOW + r*Output_w + c;
	static ap_uint<256> local_buf0[Tc/8+2];
	static ap_uint<256> local_buf1[Tc/8+2];
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

void intra_pingpong_wrapper(ap_uint<256> *Input,ap_uint<256> *Weight, float output_buffer[Tm][Tr][Tc],float beta_buffer[MAX_BETA_LENGTH],
								 float input_buffer0[Tn][OnChipIB_Height][OnChipIB_Width],float input_buffer1[Tn][OnChipIB_Height][OnChipIB_Width],
								 int IFM_num,int Input_w, int IW_align_256b, int Input_h,int OFM_num,int Ksize,int Kstride,
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
				copy_input_weight(Input,Weight,IFM_num,Input_w,IW_align_256b,Input_h,Ksize,Kstride,TMP_R,TMP_C,TMP_M, n,
					TM_MIN,TN,TRow,TCol,Padding,input_buffer1,weight_buffer1, n1, n < IFM_num,1,(TMP_M==0)&&(n==0),IHxIW,KxK,IFM_numxKxK,LayerType);
				compute(input_buffer0,output_buffer,weight_buffer0,beta_buffer, n0,Ksize,Kstride,TMP_M,TM_MIN,TR_MIN,TC_MIN, n!=0);
				pingpong = 0;
			}else
			{
				copy_input_weight(Input,Weight,IFM_num,Input_w,IW_align_256b,Input_h,Ksize,Kstride,TMP_R,TMP_C,TMP_M, n,
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

			copy_input_weight(Input,Weight,IFM_num,Input_w,IW_align_256b,Input_h,Ksize,Kstride,TMP_R,TMP_C,TMP_M,TMP_M,
				TM_MIN,TM,TRow,TCol,0,input_buffer0,weight_buffer0,NOP,input_flag,0,0,IHxIW,KxK,IFM_numxKxK,LayerType);
			pool_yolo2(input_buffer1,output_buffer,Ksize,Kstride,TX_MIN_next[0],TR_MIN,TC_MIN,process_flag);
		}else
		{
			TMP_X_next[0] = tmp_x;
			TX_MIN_next[0] = tmp_tx_min;
			tmp_x = TMP_M;
			tmp_tx_min = TM_MIN;

			copy_input_weight(Input,Weight,IFM_num,Input_w,IW_align_256b,Input_h,Ksize,Kstride,TMP_R,TMP_C,TMP_M,TMP_M,
				TM_MIN,TM,TRow,TCol,0,input_buffer1,weight_buffer1,NOP,input_flag,0,0,IHxIW,KxK,IFM_numxKxK,LayerType);
			pool_yolo2(input_buffer0,output_buffer,Ksize,Kstride,TX_MIN_next[0],TR_MIN,TC_MIN,process_flag);
		}

	}

}

void beta_copy(float beta_buffer[MAX_BETA_LENGTH], ap_uint<256> *Beta, int OFM_num)
{
	static ap_uint<256> local_buf[MAX_BETA_LENGTH/8+2];
	int loop_cnts = OFM_num >> 3;
	if(OFM_num & 0x7)
		loop_cnts++;

	for(int t = 0; t < loop_cnts; t++)
	{
#pragma HLS PIPELINE II=1
		local_buf[t].range(255, 0) = Beta[t].range(255, 0);
//		memcpy(local_buf[t], (float *)(Beta + t*8), 8*sizeof(float));
	}

	uint8_t cnt = 1;
	uint8_t bn = 0;
	ap_uint<256> tmp_256b;
	tmp_256b.range(255, 0) = local_buf[0].range(255, 0);
//	memcpy(tmp_256b, local_buf[0], 8*sizeof(float));
	for(int t = 0; t < OFM_num; t++)
	{
#pragma HLS PIPELINE II=1
		union {float f; uint32_t u;} tmp_32;
		tmp_32.u = tmp_256b.range(bn*32+31,bn*32);
		beta_buffer[t] = tmp_32.f;
//		beta_buffer[t] = tmp_256b[bn];
		bn++;
		if(bn==8)
		{
			bn = 0;
			tmp_256b.range(255, 0) = local_buf[cnt].range(255, 0);
//			memcpy(tmp_256b, local_buf[cnt], 8*sizeof(float));
			cnt++;
		}
	}



}

void YOLO2_FPGA(ap_uint<256> *Input,ap_uint<256> *Output,ap_uint<256> *Weight,ap_uint<256> *Beta, int IFM_num, int OFM_num,
							   int Ksize, int Kstride,
							   int Input_w, int Input_h, int Output_w, int Output_h, int Padding, bool IsNL, bool IsBN,
							   int TM, int TN, int TR, int TC,
							   int OFM_num_bound, int mLoopsxTM, int mLoops_a1xTM, int LayerType)
{
#pragma HLS INTERFACE m_axi depth=512 port=Input    offset=slave bundle=DATA_BUS num_read_outstanding=1 num_write_outstanding=1 max_read_burst_length=32 max_write_burst_length=32
#pragma HLS INTERFACE m_axi depth=512 port=Output    offset=slave bundle=DATA_BUS num_read_outstanding=1 num_write_outstanding=1 max_read_burst_length=32 max_write_burst_length=32
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

	const int TRow = (TR-1)*Kstride+Ksize;
	const int TCol = (TC-1)*Kstride+Ksize;
	const int KxK = Ksize*Ksize;
	const int IFM_numxKxK = IFM_num*KxK;
	//const int IHxIW   = Input_h*Input_w;
	//const int OHxOW = Output_h*Output_w;

	uint16_t IW_align_256b = (Input_w >> 3) << 3;
	if(Input_w & 0x7)
		IW_align_256b += 8;
	uint16_t OW_align_256b = (Output_w >> 3) << 3;
	if(Output_w & 0x7)
		OW_align_256b += 8;

	assert((OW_align_256b%8==0)&&(IW_align_256b%8==0));
	const int OHxOW = Output_h*OW_align_256b;
	const int IHxIW = Input_h*IW_align_256b;

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
		beta_copy( beta_buffer, Beta, OFM_num);
		//memcpy(beta_buffer,Beta, OFM_num*sizeof(float));

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
									IFM_num, Input_w, IW_align_256b, Input_h, OFM_num, Ksize, Kstride,
									r, c, m, TM_MIN, TR_MIN, TC_MIN, TN, TRow, TCol, Padding,IHxIW,KxK,IFM_numxKxK,LayerType,TM, m1,TM_MIN1, pingpongm, input_flag, process_flag);

					write_back_output_reorg(output_buffer,Output, r, c, m0[0],OW_align_256b,Output_h, TM_MIN0[0], TR_MIN, TC_MIN, OHxOW, IsNL, write_flag);
					pingpongm = 1;
				}else
				{
					intra_pingpong_wrapper(Input,Weight,output_buffer,beta_buffer,input_buffer0,input_buffer1,
									IFM_num, Input_w, IW_align_256b, Input_h, OFM_num, Ksize, Kstride,
									r, c, m, TM_MIN, TR_MIN, TC_MIN, TN, TRow, TCol, Padding,IHxIW,KxK,IFM_numxKxK,LayerType,TM, m0,TM_MIN0, pingpongm, input_flag, process_flag);

					write_back_output_reorg(output_buffer1,Output, r, c, m1[0],OW_align_256b,Output_h, TM_MIN1[0], TR_MIN, TC_MIN, OHxOW, IsNL, write_flag);
					pingpongm = 0;
				}

			}
		}
	}
}
