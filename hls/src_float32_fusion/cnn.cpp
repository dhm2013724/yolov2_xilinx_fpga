
#include "cnn.h"

void input_load(float *input,float input_buffer[Tn][IB_HxW], uint16_t r, uint16_t c, uint16_t n, uint8_t Kstride, uint8_t ksize,
		 uint8_t Padding, uint16_t TR_MIN, uint16_t TC_MIN,
		 uint16_t Input_w, uint16_t Input_h, uint8_t TN_MIN, uint32_t IHxIW, float pad_val,
		 bool r_init_en, bool c_init_en, uint16_t TCol_MIN_r)
{
	static uint16_t tl_x, br_x, px_l, px_r, TCol_MIN;
	static uint16_t tl_y, br_y, px_b, px_t, TRow_MIN;
	static int16_t Coffset, Roffset;
	static uint16_t px_lr, col_len, px_tb, row_len;
	static bool IsColCont, IsRowCont;

	if(c_init_en){
		Coffset = c*Kstride - Padding;
		// TCol_MIN = (TC_MIN-1)*Kstride + ksize;
		TCol_MIN = TCol_MIN_r;

		if(Coffset < 0){
			tl_x = 0;
			px_l = -Coffset;
		}
		else{
			tl_x = Coffset;
			px_l = 0;
		}

		br_x = Coffset + TCol_MIN_r -1;
		if(br_x < Input_w){
			px_r = 0;
		}else{
			px_r = br_x - Input_w + 1;
		}
		px_lr = (px_l + px_r);
		col_len = TCol_MIN_r - px_lr;
		IsColCont = (col_len == Input_w);
	}

	if(r_init_en){
		Roffset = r*Kstride - Padding;
		TRow_MIN = (TR_MIN-1)*Kstride + ksize;

		if(Roffset < 0){
			tl_y = 0;
			px_t = -Roffset;
		}
		else{
			tl_y = Roffset;
			px_t = 0;
		}
		br_y = Roffset + TRow_MIN -1;
		if(br_y < Input_h){
			px_b = 0;
		}else{
			px_b = br_y - Input_h + 1;
		}
		px_tb = (px_t + px_b);
		row_len = TRow_MIN - px_tb;
		IsRowCont = (row_len == Input_h);
	}

	bool IsAllCont = IsRowCont && IsColCont;

	uint32_t burstlen;
	uint32_t offset_base;
	uint16_t T1_MIN, T2_MIN;
	if(IsAllCont){
		burstlen = TN_MIN*row_len*col_len;
		offset_base = n*IHxIW;
		T1_MIN = 1;
		T2_MIN = 1;
	}else if(IsColCont){
		burstlen = row_len*col_len;
		offset_base = n*IHxIW + tl_y*Input_w;
		T1_MIN = TN_MIN;
		T2_MIN = 1;
	}else{
		burstlen = col_len;
		offset_base = n*IHxIW + tl_y*Input_w + tl_x;
		T1_MIN = TN_MIN;
		T2_MIN = row_len;
	}

        uint16_t t1, t2, t3;
        t1 = 0; t2 = 0; t3 = 0;
        for( uint16_t t1m = 0;t1m < T1_MIN; t1m++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tn)
        for( uint16_t t2m = 0;t2m < T2_MIN; t2m++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=IB_H)
                uint32_t offset_ex = offset_base;
                t3 = 0;
                if(IsAllCont){
                        t1 = 0; t2 = 0;
                }else if(IsColCont){
                        t1 = t1m; t2 = 0;
                        offset_ex += (t1m*IHxIW);
                }else{
                        t1 = t1m; t2 = t2m;
                        offset_ex += ( t1m*IHxIW + t2m*Input_w);
                }

                float *ifm_addr = input + offset_ex;
                for(uint32_t tbl=0; tbl<burstlen; tbl++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=TnxIB_HxIB_W)
#pragma HLS PIPELINE II=1
					uint32_t ifm_idx = (t2+px_t)*TCol_MIN + t3+px_l;
					input_buffer[t1][ifm_idx] = ifm_addr[tbl];

					t3++;
					if(t3==col_len){
						t3 = 0;
						t2++;
						if(t2==row_len){
								t2 = 0;
								t1++;
					}}
                }

        }}

	if(c_init_en&&px_lr){
		for(uint16_t t2 = 0;t2 < row_len; t2++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=IB_H)
		for(uint16_t t3 = 0;t3 < px_lr; t3++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=IB_W)
#pragma HLS PIPELINE II=1
			uint16_t t2_idx = t2+px_t;
			uint16_t t3_idx;
			if(t3 < px_l)
				t3_idx = t3;
			else
				t3_idx = t3+col_len;
			uint32_t ifm_idx = t2_idx*TCol_MIN + t3_idx;

			for(uint16_t t1 = 0;t1 < Tn; t1++){
//				if(t1<TN_MIN)
					input_buffer[t1][ifm_idx] = pad_val;
			}
		}}
	}

//	bool roc_init_en = r_init_en || c_init_en;
	if(px_tb){//just now, I dont know why could not set (roc_init_en&&px_tb), result dismatched.
		for(uint16_t t2 = 0;t2 < px_tb; t2++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=IB_H)
		for(uint16_t t3 = 0;t3 < TCol_MIN; t3++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=IB_W)
#pragma HLS PIPELINE II=1
			uint16_t t2_idx;
			if(t2 < px_t)
				t2_idx = t2;
			else
				t2_idx = t2+row_len;
			uint32_t ifm_idx = t2_idx*TCol_MIN + t3;

			for(uint16_t t1 = 0;t1 < Tn; t1++){
				if(t1<TN_MIN)
					input_buffer[t1][ifm_idx] = pad_val;
			}
		}}
	}

//	uint8_t TN_left = Tn - TN_MIN;
//	if(IsPad0ExtraFM&&TN_left){
//		for(uint16_t t2 = 0;t2 < TRow_MIN; t2++){
//DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=IB_H)
//		for(uint16_t t3 = 0;t3 < TCol_MIN; t3++){
//DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=IB_W)
//#pragma HLS PIPELINE II=1
//			uint32_t ifm_idx = t2*TCol_MIN + t3;
//			for(uint16_t t1 = 0;t1 < Tn; t1++){
//				if(t1>=TN_MIN)
//					input_buffer[t1][ifm_idx] = 0;
//			}
//		}}
//	}

}

const uint32_t TmxTnxKxK = Tm*Tn*K*K;

void weight_load_reorg(float *Weight,float weight_buffer[Tm][Tn][K][K],bool w_enable,uint16_t m,uint16_t n, uint32_t IFM_numxKxK, uint8_t KxK, uint8_t Ksize, uint8_t TM_MIN, uint8_t TN_MIN)
{
	if(!w_enable)
		return;

	static uint32_t Woffset;

	assert((TM_MIN > 0)&&(TM_MIN <= Tm));
	assert((TN_MIN > 0)&&(TN_MIN <= Tn));
	assert((KxK > 0)&&(KxK <= K*K));

	if(m==0&&n==0)
		Woffset = 0;

	float *weight_addr = Weight + Woffset;
	uint32_t burstlen = TM_MIN*TN_MIN*KxK;
	Woffset += burstlen;
	assert(burstlen <= TmxTnxKxK);

	uint8_t t1,t2,t3,t4;
	for(t3 = 0;t3 <Ksize; t3++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=K)
	for(t4 = 0;t4 <Ksize; t4++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=K)
#pragma HLS PIPELINE II=1
	for(t1 = 0;t1 < Tm; t1++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tm)
	for(t2 = 0;t2 < Tn; t2++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tn)
		weight_buffer[t1][t2][t3][t4] = 0;

	}}}}

	t1 = 0; t2 = 0; t3 = 0; t4 = 0;
	for(uint16_t t=0; t<burstlen; t++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=TmxTnxKxK)
#pragma HLS PIPELINE II=1
		float tmp_in = weight_addr[t];
		weight_buffer[t1][t2][t3][t4] = tmp_in;

		t2++;
		if(t2==TN_MIN){
			t2=0; t1++;
			if(t1==TM_MIN){
				t1=0; t4++;
				if(t4==Ksize){
					t4=0; t3++;
				}
			}
		}
	}

}

void copy_input_weight(float *input,float *Weight,int IFM_num,int Input_w, int Input_h, bool r_init_en, bool c_init_en, uint16_t TCol_MIN,
		int Ksize,int Kstride,int r,int c,int m,int n,
		int TM_MIN,int TN,int TR_MIN,int TC_MIN, int Padding,float input_buffer[Tn][IB_HxW],float weight_buffer[Tm][Tn][K][K],
		bool weight_load_enable,const int IHxIW,const int KxK,const int IFM_numxKxK,const int LayerType, bool enable)
{
	if(!enable)
		return;

	int TN_MIN = MIN(TN, IFM_num-n);
	bool IsPad0ExtraFM = (LayerType==0);
	float pad_val = 0.0f;
	if(!IsPad0ExtraFM)
		pad_val = -1024*1024;

	input_load(input, input_buffer, r, c, n, Kstride, Ksize, Padding, TR_MIN, TC_MIN, Input_w, Input_h, TN_MIN, IHxIW, pad_val, r_init_en, c_init_en, TCol_MIN);
	weight_load_reorg(Weight,weight_buffer,weight_load_enable,m,n,IFM_numxKxK,KxK,Ksize,TM_MIN,TN_MIN);
}

void conv2d_tile(float input_buffer[Tn][IB_HxW],float output_buffer[Tm][TrxTc],
		float weight_buffer[Tm][Tn][K][K], uint16_t n_next,
		uint8_t Ksize, uint8_t Kstride, uint16_t TCol_MIN,
		uint16_t TM_MIN, uint16_t TR_MIN, uint16_t TC_MIN,bool enable)
{
	if(!enable)
	{
		return;
	}

//	uint16_t TCol_MIN_l = (TC_MIN-1)*Kstride + Ksize;

	uint8_t i,j,tm,tn;
	uint16_t tr,tc;
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
					uint32_t ifm_idx = (Kstride*tr+i)*TCol_MIN + Kstride*tc+j;
					for(tm = 0;tm < Tm;tm++)
					{
#pragma HLS DEPENDENCE variable=output_buffer inter false

						if(i==0&&j==0&&ne0)
							partial_add[tm] = 0;
						else
							partial_add[tm] = output_buffer[tm][tr*TC_MIN + tc];

						for(tn = 0;tn <Tn;tn++)
						{
							partial_mul[tm][tn] = weight_buffer[tm][tn][i][j]*input_buffer[tn][ifm_idx];
						}

						float partial_sum = 0;
						for(tn = 0;tn <Tn;tn++)
						{
							 partial_sum += partial_mul[tm][tn];
						}
						output_buffer[tm][tr*TC_MIN + tc] = partial_add[tm] + partial_sum;
					}

				}
}

float postproc(float ofm_in, float bias_in, bool IsBias, bool IsNL){
	float tmp_out, tmp, tmp0;
	tmp = ofm_in;
	if(IsBias)
		tmp0 = tmp + bias_in;
	else
		tmp0 = tmp;

	if(IsNL&&(tmp0<0.0f))
	{
		tmp_out = tmp0*0.1f;
	}
	else
		tmp_out = tmp0;

	return tmp_out;
}


const uint32_t OFM_BLmax = Tm*Tr*Tc;
const uint32_t T12MINmax = Tm*Tr;

void aff_maxpool_2x2s2(float ofm_buf[Tm][TrxTc], uint16_t TCol_MIN, uint16_t TR_MIN, uint16_t TC_MIN, bool enable)
{
	if(!enable)
		return;

	assert(TR_MIN>0);
	assert(TC_MIN>0);

	float tmp[Tm];
#pragma HLS ARRAY_PARTITION variable=tmp complete dim=1

	for(uint16_t tr = 0;tr < TR_MIN;tr++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=TrxTc)
	for(uint16_t tc = 0;tc < TC_MIN;tc++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=TrxTc)
	for(uint16_t i = 0;i < 2; i++)
	for(uint16_t j = 0;j < 2; j++)
	{
#pragma HLS PIPELINE II=2
		uint32_t ifm_idx = (tr*2+i)*TCol_MIN + (tc*2+j);
		bool init_flag = (i==0&&j==0);
		bool wb_flag = (i==1&&j==1);

		for(uint16_t of = 0; of < Tm; of++)
		{//this DEPENDENCE, maybe useful
#pragma HLS DEPENDENCE variable=ofm_buf inter false

			float tmp_in = ofm_buf[of][ifm_idx];
			bool max_flag = (tmp_in > tmp[of]);
			if(init_flag||max_flag)
				tmp[of] = tmp_in;

			if(wb_flag)
				ofm_buf[of][tr*TC_MIN + tc] = tmp[of];
		}
	}}}

}

void write_back_output_reorg(float output_buffer[Tm][TrxTc], float bias_buffer[MAX_BETA_LENGTH], /*float bias_buffer[MAX_BETA_LENGTH],*/
		float *Output, uint16_t r_r,uint16_t c_r,uint16_t m,
		uint16_t ofm_num, uint16_t ofm_h_r, uint16_t ofm_w_r,
		uint8_t TM_MIN, uint16_t TR_MIN_r, uint16_t TC_MIN_r, uint32_t OHxOW_r, bool IsNL, bool IsBias, bool enable, bool IsPoolAff)
{
	if(!enable)
		return;

	assert((TM_MIN >0)&&(TM_MIN <=Tm));
	// assert((TR_MIN >0)&&(TR_MIN <=Tr));
	// assert((TC_MIN >0)&&(TC_MIN <=Tc));
	uint16_t r, c, TR_MIN, TC_MIN;
	uint32_t OHxOW;
	uint16_t ofm_h, ofm_w;

	if(IsPoolAff){
		assert(TR_MIN%2==0);
		assert(TC_MIN%2==0);
		r = r_r/2; c = c_r/2; TR_MIN = TR_MIN_r/2; TC_MIN = TC_MIN_r/2;
		ofm_h = ofm_h_r/2; ofm_w = ofm_w_r/2;
		OHxOW = OHxOW_r/4;
	}else
	{
		r = r_r; c = c_r; TR_MIN = TR_MIN_r; TC_MIN = TC_MIN_r;
		ofm_h = ofm_h_r; ofm_w = ofm_w_r;
		OHxOW = OHxOW_r;
	}

	aff_maxpool_2x2s2(output_buffer, TC_MIN_r, TR_MIN, TC_MIN, IsPoolAff);

	uint32_t offset = m*OHxOW + r*ofm_w + c;

	bool IsColCont = (ofm_w==TC_MIN);//equal FMCont
	bool IsRowCont = (ofm_h==TR_MIN);
	bool IsAllCont = IsColCont && IsRowCont;

	int burstlen;
	uint16_t T1_MIN, T2_MIN;
	if(IsAllCont){
		burstlen = TM_MIN*TR_MIN*TC_MIN;
		T1_MIN = 1;
		T2_MIN = 1;
	}else if(IsColCont){
		burstlen = TR_MIN*TC_MIN;
		T1_MIN = TM_MIN;
		T2_MIN = 1;
	}else{
		burstlen = TC_MIN;
		T1_MIN = TM_MIN;
		T2_MIN = TR_MIN;
	}

	uint32_t Tcomb_TC = T1_MIN*T2_MIN;

	uint16_t tr,tm,tc;
	uint16_t t1, t2;

	t1 = 0; t2 = 0;
	for(uint32_t tcb = 0; tcb < Tcomb_TC; tcb++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=T12MINmax)
		tc = 0;
		if(IsAllCont){
			tm = 0;
			tr = 0;
		}else if(IsColCont){
			tm = t1;
			tr = 0;
		}else{
			tm = t1;
			tr = t2;
		}

		uint32_t ofm_offset0 = offset;
		uint32_t ofm_offset;
		if(IsAllCont){
			ofm_offset = ofm_offset0;
		}else if(IsColCont){
			ofm_offset = ofm_offset0 + t1*OHxOW;
		}else{
			ofm_offset = ofm_offset0 + t1*OHxOW + t2*ofm_w;
		}

		float *OFM_base = Output + ofm_offset;
		for(int tbl_min=0; tbl_min < burstlen; tbl_min++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=OFM_BLmax)
#pragma HLS PIPELINE II=1
			float tmp_out = postproc(output_buffer[tm][tr*TC_MIN + tc], bias_buffer[m+tm], IsBias, IsNL);

			tc++;
			if(tc==TC_MIN){
				tc = 0;
				tr++;
				if(tr==TR_MIN){
					tr = 0;
					tm++;
				}
			}

			OFM_base[tbl_min] = tmp_out;
		}

		t2++;
		if(t2==T2_MIN){
			t2=0;
			t1++;
		}
	}
}


void pool_yolo2(float Input[Tn][IB_HxW],float Output[Tm][TrxTc],
		  uint8_t Ksize, uint8_t Kstride, uint16_t TCol_MIN, uint16_t TM_MIN, uint16_t TR_MIN, uint16_t TC_MIN, bool enable)
{
	if(!enable)
		return;

	float tmp[Tn];
#pragma HLS ARRAY_PARTITION variable=tmp complete dim=1

	for(uint16_t tr = 0;tr < TR_MIN;tr++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tr)
	for(uint16_t tc = 0;tc < TC_MIN;tc++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tc)
	for(uint8_t i =0;i < Ksize; i++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=K)
	for(uint8_t j = 0;j < Ksize; j++)
	{
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=K)
#pragma HLS PIPELINE II=2
		uint32_t ifm_idx = (tr*Kstride+i)*TCol_MIN + tc*Kstride+j;
		for(uint16_t of = 0; of < Tn; of++)
		{
			if(i==0&&j==0)
				tmp[of] = -1024*1024;

			float tmp_in = Input[of][ifm_idx];
			if(tmp_in > tmp[of])
				tmp[of] = tmp_in;

			if(i==1&&j==1)
				Output[of][tr*TC_MIN + tc] = tmp[of];
		}
	}}}}
}

void load_compute_wrapper(float *ifm, float *weight, float ofm_buffer[Tm][TrxTc], int ksize, int kstride, int ifm_num, int ifm_w, int ifm_h,int ofm_num,
	 int ofm_h, int ofm_w, int pad_int, int ltype, int IHW, int KK, int INumxKK, int TM, int TN, int TR, int TC, int tm_r, int tr_r, int tc_r,
	 int tx_n1[3],int TX_MIN_n1[3],bool pp,bool in_flag,bool proc_flag, bool LoadBias, int NTif, uint8_t lmode, bool enable, bool r_init_en, bool c_init_en)
{
	static float ifm_buffer0[Tn][IB_HxW];
#pragma HLS ARRAY_PARTITION variable=ifm_buffer0 complete dim=1
	static float ifm_buffer1[Tn][IB_HxW];
#pragma HLS ARRAY_PARTITION variable=ifm_buffer1 complete dim=1
	static float weight_buffer0[Tm][Tn][K][K];
#pragma HLS ARRAY_PARTITION variable=weight_buffer0 complete dim=1
#pragma HLS ARRAY_PARTITION variable=weight_buffer0 complete dim=2
	static float weight_buffer1[Tm][Tn][K][K];
#pragma HLS ARRAY_PARTITION variable=weight_buffer1 complete dim=1
#pragma HLS ARRAY_PARTITION variable=weight_buffer1 complete dim=2

	static uint16_t tx_n[2][3];//tmtrtc, between load2compute
#pragma HLS ARRAY_PARTITION variable=tx_n complete dim=1
#pragma HLS ARRAY_PARTITION variable=tx_n complete dim=2
	static uint16_t TX_MIN_n[2][3];
#pragma HLS ARRAY_PARTITION variable=TX_MIN_n complete dim=1
#pragma HLS ARRAY_PARTITION variable=TX_MIN_n complete dim=2

	if(!enable)
		return ;

	static uint16_t TR_MIN, TC_MIN, TM_MIN;
	static uint16_t TCol_MIN;
	static uint16_t TCol_MIN0[1], TCol_MIN1[1];

	if(c_init_en){
		TC_MIN = MIN(TC,ofm_w-tc_r);
		TCol_MIN = (TC_MIN-1)*kstride + ksize;
	}
	if(r_init_en){
		TR_MIN = MIN(TR,ofm_h-tr_r);
	}
	TM_MIN = MIN(TM,ofm_num-tm_r);

	if(lmode==0){
		if(pp){
			if(ltype == 0)
			{
				copy_input_weight(ifm,weight,ifm_num,ifm_w,ifm_h, r_init_en, c_init_en, TCol_MIN, ksize,kstride,tr_r,tc_r,tm_r, 0,
					TM_MIN,TN,TR_MIN,TC_MIN,pad_int,ifm_buffer0,weight_buffer0,1,IHW,KK,INumxKK,ltype,in_flag);
				conv2d_tile(ifm_buffer1, ofm_buffer, weight_buffer1, 0, ksize, kstride, TCol_MIN1[0],
						TX_MIN_n[1][0], TX_MIN_n[1][1], TX_MIN_n[1][2], proc_flag);
			}else if(ltype == 1)
			{
				copy_input_weight(ifm,weight,ifm_num,ifm_w,ifm_h, r_init_en, c_init_en, TCol_MIN, ksize,kstride,tr_r,tc_r,tm_r,tm_r,
										TM_MIN,TM,TR_MIN,TC_MIN,0,ifm_buffer0,weight_buffer0,0,IHW,KK,INumxKK,ltype,in_flag);
				pool_yolo2(ifm_buffer1, ofm_buffer, ksize, kstride, TCol_MIN1[0], TX_MIN_n[1][0], TX_MIN_n[1][1], TX_MIN_n[1][2], proc_flag);
			}

			TCol_MIN0[0] = TCol_MIN;
			tx_n[0][0] = tm_r; tx_n[0][1] = tr_r; tx_n[0][2] = tc_r;
			TX_MIN_n[0][0] = TM_MIN; TX_MIN_n[0][1] = TR_MIN; TX_MIN_n[0][2] = TC_MIN;

			tx_n1[0] = tx_n[1][0]; tx_n1[1] = tx_n[1][1]; tx_n1[2] = tx_n[1][2];
			TX_MIN_n1[0] = TX_MIN_n[1][0]; TX_MIN_n1[1] = TX_MIN_n[1][1]; TX_MIN_n1[2] = TX_MIN_n[1][2];
		}else{
			if(ltype == 0)
			{
				copy_input_weight(ifm,weight,ifm_num,ifm_w,ifm_h, r_init_en, c_init_en, TCol_MIN, ksize,kstride,tr_r,tc_r,tm_r, 0,
					TM_MIN,TN,TR_MIN,TC_MIN,pad_int,ifm_buffer1,weight_buffer1,1,IHW,KK,INumxKK,ltype,in_flag);
				conv2d_tile(ifm_buffer0, ofm_buffer, weight_buffer0, 0, ksize, kstride, TCol_MIN0[0],
						TX_MIN_n[0][0], TX_MIN_n[0][1], TX_MIN_n[0][2], proc_flag);
			}else if(ltype == 1)
			{
				copy_input_weight(ifm,weight,ifm_num,ifm_w,ifm_h, r_init_en, c_init_en, TCol_MIN, ksize,kstride,tr_r,tc_r,tm_r,tm_r,
										TM_MIN,TM,TR_MIN,TC_MIN,0,ifm_buffer1,weight_buffer1,0,IHW,KK,INumxKK,ltype,in_flag);
				pool_yolo2(ifm_buffer0, ofm_buffer, ksize, kstride, TCol_MIN0[0], TX_MIN_n[0][0], TX_MIN_n[0][1], TX_MIN_n[0][2], proc_flag);
			}

			TCol_MIN1[0] = TCol_MIN;
			tx_n[1][0] = tm_r; tx_n[1][1] = tr_r; tx_n[1][2] = tc_r;
			TX_MIN_n[1][0] = TM_MIN; TX_MIN_n[1][1] = TR_MIN; TX_MIN_n[1][2] = TC_MIN;

			tx_n1[0] = tx_n[0][0]; tx_n1[1] = tx_n[0][1]; tx_n1[2] = tx_n[0][2];
			TX_MIN_n1[0] = TX_MIN_n[0][0]; TX_MIN_n1[1] = TX_MIN_n[0][1]; TX_MIN_n1[2] = TX_MIN_n[0][2];
		}
	}else{
			if(ltype == 0)
			{
				bool pp_tn = 1;
				uint16_t tn_n[2];
				uint16_t tn_r = 0;
				for(int tn = 0; tn < NTif+1; tn++)
				{
					bool in_flag = tn < NTif;
					bool proc_flag = tn > 0;
					// uint16_t tn_r = tn*TN;
					if(pp_tn){
						copy_input_weight(ifm,weight,ifm_num,ifm_w,ifm_h, r_init_en, c_init_en, TCol_MIN, ksize,kstride,tr_r,tc_r,tm_r, tn_r,
							TM_MIN,TN,TR_MIN,TC_MIN,pad_int,ifm_buffer0,weight_buffer0,1,IHW,KK,INumxKK,ltype,in_flag);
						conv2d_tile(ifm_buffer1, ofm_buffer, weight_buffer1, tn_n[1], ksize, kstride, TCol_MIN1[0],
								TX_MIN_n[1][0], TX_MIN_n[1][1], TX_MIN_n[1][2], proc_flag);

						TCol_MIN0[0] = TCol_MIN;
						tn_n[0] = tn_r;
						tx_n[0][0] = tm_r; tx_n[0][1] = tr_r; tx_n[0][2] = tc_r;
						TX_MIN_n[0][0] = TM_MIN; TX_MIN_n[0][1] = TR_MIN; TX_MIN_n[0][2] = TC_MIN;
						tn_r += TN;

						pp_tn = 0;
					}else{
						copy_input_weight(ifm,weight,ifm_num,ifm_w,ifm_h, r_init_en, c_init_en, TCol_MIN, ksize,kstride,tr_r,tc_r,tm_r, tn_r,
							TM_MIN,TN,TR_MIN,TC_MIN,pad_int,ifm_buffer1,weight_buffer1,1,IHW,KK,INumxKK,ltype,in_flag);
						conv2d_tile(ifm_buffer0, ofm_buffer, weight_buffer0, tn_n[0], ksize, kstride, TCol_MIN0[0],
								TX_MIN_n[0][0], TX_MIN_n[0][1], TX_MIN_n[0][2], proc_flag);

						TCol_MIN1[0] = TCol_MIN;
						tn_n[1] = tn_r;
						tx_n[1][0] = tm_r; tx_n[1][1] = tr_r; tx_n[1][2] = tc_r;
						TX_MIN_n[1][0] = TM_MIN; TX_MIN_n[1][1] = TR_MIN; TX_MIN_n[1][2] = TC_MIN;
						tn_r += TN;

						pp_tn = 1;
					}
				}
			}

			if(pp){
				tx_n1[0] = tx_n[1][0]; tx_n1[1] = tx_n[1][1]; tx_n1[2] = tx_n[1][2];
				TX_MIN_n1[0] = TX_MIN_n[1][0]; TX_MIN_n1[1] = TX_MIN_n[1][1]; TX_MIN_n1[2] = TX_MIN_n[1][2];
			}else{
				tx_n1[0] = tx_n[0][0]; tx_n1[1] = tx_n[0][1]; tx_n1[2] = tx_n[0][2];
				TX_MIN_n1[0] = TX_MIN_n[0][0]; TX_MIN_n1[1] = TX_MIN_n[0][1]; TX_MIN_n1[2] = TX_MIN_n[0][2];
			}
	}
}

void YOLO2_FPGA(float *Input,float *Output,float *Weight,float *Beta, int IFM_num, int OFM_num,
							   int Ksize, int Kstride, int Input_w, int Input_h, int Output_w, int Output_h, int Padding, bool IsNL,
							   int TM, int TN, int TR, int TC,
							   int32_t NToy, int32_t NTox, int32_t NTof, int32_t NTcomb, int32_t NTif, uint8_t lmode, int32_t NTcomb_l, int LayerType,
							   bool IsPoolAff)
{
#pragma HLS INTERFACE m_axi depth=512 port=Input    offset=slave bundle=DATA_BUS num_read_outstanding=1 num_write_outstanding=1 max_read_burst_length=256 max_write_burst_length=256
#pragma HLS INTERFACE m_axi depth=512 port=Output    offset=slave bundle=DATA_BUS num_read_outstanding=1 num_write_outstanding=1 max_read_burst_length=256 max_write_burst_length=256
#pragma HLS INTERFACE m_axi depth=512 port=Weight offset=slave bundle=DATA_BUS1 num_read_outstanding=1 max_read_burst_length=256
#pragma HLS INTERFACE m_axi depth=512 port=Beta   offset=slave bundle=DATA_BUS1 num_read_outstanding=1 max_read_burst_length=256

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
#pragma HLS INTERFACE s_axilite register port=TM bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=TN bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=TR bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=TC bundle=CTRL_BUS

#pragma HLS INTERFACE s_axilite register port=NToy bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=NTox bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=NTof bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=NTcomb bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=NTif bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=lmode bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=NTcomb_l bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=LayerType bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=IsPoolAff bundle=CTRL_BUS

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
//	assert((TR > 0)&&(TR <= Tr));
//	assert((TC > 0)&&(TC <= Tc));

	const int KxK = Ksize*Ksize;
	const int IFM_numxKxK = IFM_num*KxK;

	const int IHxIW   = Input_h*Input_w;
	const int OHxOW = Output_h*Output_w;

	static float bias_buffer[MAX_BETA_LENGTH];
	static float ofm_buffer0[Tm][TrxTc];
#pragma HLS ARRAY_PARTITION variable=ofm_buffer0 complete dim=1
	static float ofm_buffer1[Tm][TrxTc];
#pragma HLS ARRAY_PARTITION variable=ofm_buffer1 complete dim=1

	int tx_n10[3], tx_n11[3];//tmtrtc, between compute2store
#pragma HLS ARRAY_PARTITION variable=tx_n10 complete dim=1
#pragma HLS ARRAY_PARTITION variable=tx_n11 complete dim=1
	int TX_MIN_n10[3], TX_MIN_n11[3];
#pragma HLS ARRAY_PARTITION variable=TX_MIN_n10 complete dim=1
#pragma HLS ARRAY_PARTITION variable=TX_MIN_n11 complete dim=1

	bool IsBias;
	if(LayerType==0){
		memcpy(bias_buffer,Beta, OFM_num*sizeof(float));
		IsBias = 1;
	}else{
		IsBias = 0;
	}

	uint16_t tr, tc, tm;
	uint16_t TR_MIN, TC_MIN, TM_MIN;
	uint16_t tm_r, tr_r, tc_r;
	bool in_flag, proc_flag, out_flag, lc_enable;

	tr = 0; tc = 0; tm = 0;
	tm_r = 0; tr_r = 0; tc_r = 0;
	bool pp = 1;
	bool r_init_en = 1, c_init_en = 1;
	for(int ntc = 0; ntc < NTcomb_l; ntc++){
		if(lmode==0){
			in_flag = ntc < NTcomb;
			proc_flag = (ntc > 0)&& (ntc < NTcomb + 1);
			out_flag = ntc > 1;
			lc_enable = 1;
		}else{
			in_flag = 1;
			proc_flag = 1;
			out_flag = ntc > 0;
			lc_enable = ntc < NTcomb;
		}

		if(pp){
			load_compute_wrapper(Input, Weight, ofm_buffer0, Ksize, Kstride, IFM_num, Input_w, Input_h, OFM_num, Output_h, Output_w,
				Padding, LayerType, IHxIW, KxK, IFM_numxKxK, TM, TN, TR, TC, tm_r, tr_r, tc_r,
				tx_n10, TX_MIN_n10, pp, in_flag, proc_flag, IsBias, NTif, lmode, lc_enable, r_init_en, c_init_en);

			write_back_output_reorg( ofm_buffer1, bias_buffer, Output, tx_n11[1], tx_n11[2], tx_n11[0], OFM_num, Output_h, Output_w,
					TX_MIN_n11[0], TX_MIN_n11[1], TX_MIN_n11[2], OHxOW, IsNL, IsBias, out_flag, IsPoolAff);
			pp = 0;
		}else{
			load_compute_wrapper(Input, Weight, ofm_buffer1, Ksize, Kstride, IFM_num, Input_w, Input_h, OFM_num, Output_h, Output_w,
				Padding, LayerType, IHxIW, KxK, IFM_numxKxK, TM, TN, TR, TC, tm_r, tr_r, tc_r,
				tx_n11, TX_MIN_n11, pp, in_flag, proc_flag, IsBias, NTif, lmode, lc_enable, r_init_en, c_init_en);

			write_back_output_reorg( ofm_buffer0, bias_buffer, Output, tx_n10[1], tx_n10[2], tx_n10[0], OFM_num, Output_h, Output_w,
					TX_MIN_n10[0], TX_MIN_n10[1], TX_MIN_n10[2], OHxOW, IsNL, IsBias, out_flag, IsPoolAff);
			pp = 1;
		}

		tm++;
		tm_r += TM;
		r_init_en = 0; c_init_en = 0;
		if(tm==NTof)
		{
			tm=0; tc++;
			tm_r = 0; tc_r += TC; c_init_en = 1;//update col
			if(tc==NTox)
			{
				tc=0; tr++;
				tc_r = 0; tr_r += TR; r_init_en = 1;//update row
		}}
	}
}
