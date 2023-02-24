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

#define LT_DCONV 0
#define LT_CONV  1
#define LT_AVGPOOL 2
#define LT_MAXPOOL 3

#define MIN_NEG (-1024*1024)
// #define MIN_NEG_INT16 (0x8000)
#define MIN_NEG_INT32 (0x80000000)

#ifndef _DEF_IN_MAKEFILE_
#define HW_S 3
#define K 7
#define Tn 1
#define Tm 32
#define Tr 23
#define Tc 31
#define MAX_BETA_LENGTH 2048
#endif

#define OnChipIB_Width  ((Tc-1)*HW_S+K)
#define OnChipIB_Height ((Tr-1)*HW_S+K)

#define PRAGMA_SUB(x) _Pragma (#x)
#define DO_PRAGMA(x) PRAGMA_SUB(x)

const uint32_t IB_W = OnChipIB_Width;
const uint32_t IB_H = OnChipIB_Height;
const uint32_t IB_HxW = IB_H*IB_W;
const uint32_t TnxIB_H = Tn*IB_H;	
const uint32_t TnxIB_HxIB_W = Tn*IB_H*IB_W;
const uint32_t TrxTc = Tr*Tc;

int32_t f2i32(float val_f32){
	int32_t val_i32;
	union {
		float f32;
		int32_t i32;
	}tmp_32b;

	tmp_32b.f32 = val_f32;
	val_i32 = tmp_32b.i32;

	return val_i32;
}

float i2f32(int32_t val_i32){
	float val_f32;
	union {
		float f32;
		int32_t i32;
	}tmp_32b;

	tmp_32b.i32 = val_i32;
	val_f32 = tmp_32b.f32;

	return val_f32;
}

void input_load(float *input,float input_buffer[Tn][IB_HxW], uint16_t r, uint16_t c, uint16_t n, uint8_t Kstride, uint8_t ksize, uint8_t Padding,
		 uint16_t TR_MIN, uint16_t TC_MIN,
		 uint16_t Input_w, uint16_t Input_h, uint8_t TN_MIN, uint32_t IHxIW, float pad_val, bool r_init_en, bool c_init_en, uint16_t TCol_MIN_r)
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
	uint16_t T2_MIN;
	uint32_t T12_MIN;
	if(IsAllCont){
		burstlen = TN_MIN*row_len*col_len;
		offset_base = n*IHxIW;
		T2_MIN = 1;
		T12_MIN = 1;
	}else if(IsColCont){
		burstlen = row_len*col_len;
		offset_base = n*IHxIW + tl_y*Input_w;
		T2_MIN = 1;
		T12_MIN = TN_MIN;
	}else{
		burstlen = col_len;
		offset_base = n*IHxIW + tl_y*Input_w + tl_x;
		T2_MIN = row_len;
		T12_MIN = TN_MIN*row_len;		
	}
	assert(burstlen <= TnxIB_HxIB_W);

	uint16_t t1, t2, t3;
	uint16_t t1m, t2m;
	t1 = 0; t2 = 0; t3 = 0;
	t1m = 0, t2m = 0;
	for( uint32_t t12m = 0;t12m < T12_MIN; t12m++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=TnxIB_H)
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
			uint16_t t2_idx = t2+px_t;
			uint16_t t3_idx = t3+px_l;
			// input_buffer[t1][t2+px_t][t3+px_l] = ifm_addr[tbl];
			input_buffer[t1][t2_idx*TCol_MIN_r + t3_idx] = ifm_addr[tbl];

			t3++;
			if(IsColCont)
			if(t3==col_len){
				t3 = 0;
				t2++;
				if(IsAllCont)
				if(t2==row_len){
					t2 = 0;
					t1++;
			}}	
		}

		t2m++;
		if(t2m==T2_MIN){
			t2m = 0; t1m++;
		}
	}

	if(px_tb){
		for( t2 = 0;t2 < px_tb; t2++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=IB_H)
		for( t3 = 0;t3 < TCol_MIN_r; t3++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=IB_W)
#pragma HLS PIPELINE II=1
			uint16_t t2_idx;
			if(t2 < px_t)
				t2_idx = t2;
			else
				t2_idx = t2+row_len;
			uint32_t t23_idx = t2_idx*TCol_MIN_r + t3;

			for( t1 = 0;t1 < Tn; t1++){
				if(t1 < TN_MIN){
					input_buffer[t1][t23_idx] = pad_val;
				}
			}
		}}		
	}

	if(px_lr){
		for( t2 = 0;t2 < row_len; t2++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=IB_H)
		for( t3 = 0;t3 < px_lr; t3++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=IB_W)
#pragma HLS PIPELINE II=1
			uint16_t t2_idx = t2+px_t;
			uint16_t t3_idx;
			if(t3 < px_l)
				t3_idx = t3;
			else
				t3_idx = t3+col_len;
			uint32_t t23_idx = t2_idx*TCol_MIN_r + t3_idx;

			for( t1 = 0;t1 < Tn; t1++){
				if(t1 < TN_MIN){
					input_buffer[t1][t23_idx] = pad_val;
				}
			}
		}}
	}		
}

void weight_load_reorg(float *Weight,float weight_buffer[Tm][Tn][K][K],uint16_t m,uint16_t n, uint32_t IFM_numxKxK, uint8_t KxK, uint8_t Ksize, 
			uint8_t TM_MIN, uint8_t TN_MIN, bool enable)
{
const uint32_t TmxTnxKxK = Tm*Tn*K*K;	

	if(!enable)
		return;

	uint8_t t1,t2,t3,t4;
	static uint32_t Woffset;

	assert((TM_MIN > 0)&&(TM_MIN <= Tm));
	assert((TN_MIN > 0)&&(TN_MIN <= Tn));
	assert((KxK > 0)&&(KxK <= K*K));
		
	for(t3 = 0;t3 <Ksize; t3++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=K)
	for(t4 = 0;t4 <Ksize; t4++)
	{
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=K)
#pragma HLS PIPELINE II=1
		for(t1 = 0;t1 < Tm; t1++){		
		for(t2 = 0;t2 < Tn; t2++){
			weight_buffer[t1][t2][t3][t4] = 0;
		}}
	}}	

	if(m==0&&n==0)
		Woffset = 0;

	uint32_t burstlen = TM_MIN*TN_MIN*KxK;
	assert(burstlen <= TmxTnxKxK);

	float *weight_addr = Weight + Woffset;
	Woffset += burstlen;
	t1 = 0; t2 = 0; t3 = 0; t4 = 0;
	for(uint32_t tbl = 0; tbl < burstlen; tbl++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=TmxTnxKxK)
#pragma HLS PIPELINE II=1

		weight_buffer[t1][t2][t3][t4] = weight_addr[tbl];

		t4++;
		if(t4==Ksize){
			t4=0;
			t3++;
			if(t3==Ksize){
				t3=0;
				t2++;
				if(t2==TN_MIN){
					t2=0;
					t1++;
		}}}
	}
}

void weight_load(float *Weight,float weight_buffer[Tm][Tn][K][K], uint16_t m, uint16_t n, uint32_t IFM_numxKxK, uint8_t KxK, 
			uint8_t Ksize, uint8_t TM_MIN, uint8_t TN_MIN, bool enable)
{
	if(!enable)
		return;

	uint8_t t1,t2,t3,t4;
	static float w_buf[Tm*Tn*K*K];
	const int Woffset = m*IFM_numxKxK + n*KxK;

	for(t3 = 0;t3 <Ksize; t3++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=K)
	for(t4 = 0;t4 <Ksize; t4++)
	{
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=K)
#pragma HLS PIPELINE II=1
		for(t1 = 0;t1 < Tm; t1++)
		for(t2 = 0;t2 < Tn; t2++){
			weight_buffer[t1][t2][t3][t4] = 0;
		}
	}}

	int w_buf_offset = 0;
	for(t1 = 0;t1 < TM_MIN; t1++)
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tm)
		for(t2 = 0;t2 < TN_MIN; t2++)
		{
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tn)
			memcpy((float *)(w_buf + w_buf_offset),(float *)(Weight + Woffset + t1*IFM_numxKxK + t2*KxK),KxK*sizeof(float));
			w_buf_offset += KxK;
		}

	w_buf_offset = 0;
	for(t1 = 0;t1 < TM_MIN; t1++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tm)
	for(t2 = 0;t2 < TN_MIN; t2++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tn)
	for(t3 = 0;t3 <Ksize; t3++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=K)
	for(t4 = 0;t4 <Ksize; t4++)
	{
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=K)
#pragma HLS PIPELINE II=1
		weight_buffer[t1][t2][t3][t4] =  w_buf[w_buf_offset];
		w_buf_offset++;
	}}}}
}

float postproc(float ofm_in, float bias_in, float div_kk, bool LoadBias, bool IsAVGPOOL, bool IsNL){
	float tmp_out, tmp, tmp0, tmp1;
	tmp = ofm_in;
	if(LoadBias)
		tmp0 = tmp + bias_in;
	else
		tmp0 = tmp;

	tmp1 = tmp0;

	if(IsNL)
	{
		if(tmp1 < 0.0f)
			tmp_out = tmp1*0.1f;
		else
			tmp_out = tmp1;
	}
	else
		tmp_out = tmp1;
	
	return tmp_out;
}

int16_t postproc(int32_t ofm_in, int16_t bias_in, int16_t div_kk, bool LoadBias, bool IsAVGPOOL, bool IsNL, uint8_t biasQ, uint8_t ofmQ, uint8_t avgQ, uint8_t interQ){

	int16_t tmp_out;
	int32_t tmp0, tmp1, tmp2;

	int32_t tmp = ofm_in;
	if(LoadBias){
		int32_t tmp_b_i32;
		if(biasQ > interQ)
			tmp_b_i32 = (bias_in >> (biasQ - interQ));
		else
			tmp_b_i32 = (bias_in << (interQ - biasQ));
		tmp0 = tmp + tmp_b_i32;
	}
	else
		tmp0 = tmp;

	tmp1 = tmp0;

	if(IsNL)
	{
		if(tmp1 < 0)
			tmp2 = (tmp1*0xccc) >> 15;
		else
			tmp2 = tmp1;
	}
	else
		tmp2 = tmp1;

	if(ofmQ > interQ)
		tmp_out = tmp2 << (ofmQ - interQ);
	else
		tmp_out = tmp2 >> (interQ - ofmQ);
	
	return tmp_out;
}

void copy_local_beta(float bias_buffer[MAX_BETA_LENGTH],float local_beta_buffer[Tm], uint16_t TM_MIN, uint16_t m)
{	
	for(uint16_t tm = 0;tm < TM_MIN;tm++)
	{
		local_beta_buffer[tm] = bias_buffer[m+tm];
	}
}

void write_back_output_reorg(float output_buffer[Tm][TrxTc], float bias_buffer[MAX_BETA_LENGTH], /*float bias_buffer[MAX_BETA_LENGTH],*/
		float *Output, uint16_t r,uint16_t c,uint16_t m, 
		uint16_t ofm_num, uint16_t ofm_h, uint16_t ofm_w,
		uint16_t TM_MIN, uint16_t TR_MIN, uint16_t TC_MIN, uint32_t OHxOW, bool IsNL, float div_kk, bool IsAVGPOOL, bool LoadBias, bool enable,
		uint32_t biasQ,  uint32_t ofmQ, uint32_t avgQ, uint32_t interQ, bool interQ_en)
{
	const uint32_t OFM_BLmax = 256;

	static float local_bias_buffer[Tm];

	if(!enable)
		return;

	// printf("ofmQ=%d\n", ofmQ);

	assert((TM_MIN >0)&&(TM_MIN <=Tm));
	assert((TR_MIN*TC_MIN > 0)&&(TR_MIN*TC_MIN <= Tr*Tc));

	copy_local_beta( bias_buffer, local_bias_buffer, TM_MIN, m);

	uint32_t offset = m*OHxOW + r*ofm_w + c;
	bool IsColCont = (ofm_w==TC_MIN);//equal FMCont
	bool IsRowCont = (ofm_h==TR_MIN);
	bool IsAllCont = IsColCont && IsRowCont;

	uint32_t burstlen;
	uint16_t T1_MIN, T2_MIN;
	if(IsAllCont){
		// burstlen = TM_MIN*TR_MIN*TC_MIN;
		burstlen = TM_MIN*TR_MIN*ofm_w;
		T1_MIN = 1;
		T2_MIN = 1;
	}else if(IsColCont){
		burstlen = TR_MIN*ofm_w;
		T1_MIN = TM_MIN;
		T2_MIN = 1;
	}else{
		burstlen = TC_MIN;
		T1_MIN = TM_MIN;
		T2_MIN = TR_MIN;
	}

	uint8_t add_1b;
	if(burstlen & 0xFF)
		add_1b = 1;
	else
		add_1b = 0;

	uint32_t blmax_tc =  (burstlen >> 8) + add_1b;
	uint32_t Tcomb_TC = T1_MIN*T2_MIN*blmax_tc;
	uint16_t tr,tm,tc;
	uint16_t t1, t2;
	uint32_t tbl;

	t1 = 0; t2 = 0; tbl = 0;
	for(uint32_t tcb = 0; tcb < Tcomb_TC; tcb++){

		uint32_t tbl_r = (tbl << 8);
		if(tbl==0){
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
		}

		uint32_t ofm_offset0 = offset + tbl_r;
		uint32_t ofm_offset;
		if(IsAllCont){
			ofm_offset = ofm_offset0;
		}else if(IsColCont){
			ofm_offset = ofm_offset0 + t1*OHxOW;
		}else{
			ofm_offset = ofm_offset0 + t1*OHxOW + t2*ofm_w;
		}

		float *OFM_base = Output + ofm_offset;
		uint16_t TBL_MIN = MIN_diy(OFM_BLmax, burstlen-tbl_r);

		for(int tbl_min=0; tbl_min < TBL_MIN; tbl_min++){

			if(interQ_en){
				int16_t bias_in_16b = local_bias_buffer[tm]*pow(2.0, biasQ);
				OFM_base[tbl_min] = postproc(f2i32(output_buffer[tm][tr*TC_MIN + tc]), bias_in_16b,
									div_kk*pow(2.0, avgQ), LoadBias, IsAVGPOOL, IsNL, biasQ, ofmQ, avgQ, interQ)*pow(0.5, ofmQ);
			}else{
				OFM_base[tbl_min] = postproc(output_buffer[tm][tr*TC_MIN + tc], local_bias_buffer[tm], div_kk, LoadBias, IsAVGPOOL, IsNL);
			}

			tc++;
			if(IsColCont)
			if(tc==TC_MIN){
				tc = 0;
				tr++;
				if(IsAllCont)
				if(tr==TR_MIN){
					tr = 0;
					tm++;
				}				
			}
		}

		tbl++;
		if(tbl==blmax_tc){
			tbl = 0;
			t2++;
			if(t2==T2_MIN){
				t2=0;
				t1++;
			}
		}

	}
}

void pool_tile(float Input[Tn][IB_HxW],float Output[Tm][TrxTc],
	uint8_t Ksize, uint8_t Kstride, uint8_t TM_MIN, uint16_t TR_MIN, uint16_t TC_MIN, uint16_t TCol_MIN, bool enable, uint8_t pool_mode,
	uint8_t weightQ, uint32_t ifmQ, uint8_t ComQ, bool C_SL_EN, bool interQ_en)//0:MAX; 1:AVG
{
	if(!enable)
		return;

	uint8_t i,j,of;
	uint16_t tr,tc;

	bool max_flag = (pool_mode==0);

	if(interQ_en){

		float tmp_exp_ifm = pow(2.0, ifmQ);	
		for(i =0;i < Ksize; i++)
		for(j = 0;j < Ksize; j++)
		for(tr = 0;tr < TR_MIN;tr++)
		for(tc = 0;tc < TC_MIN;tc++)
		{
			bool init_flag = ((i==0)&&(j==0));
			uint32_t t23_idx = (tr*Kstride+i)*TCol_MIN + (tc*Kstride+j);
			uint32_t ofm_rc_idx = tr*TC_MIN + tc;
			for( of = 0; of < Tn; of++)
			{
				int32_t tmp_out;
				int16_t tmp_in = Input[of][t23_idx]*tmp_exp_ifm;
				int32_t tmp_in_i32;
				if(C_SL_EN)
					tmp_in_i32 = tmp_in << ComQ;
				else
					tmp_in_i32 = tmp_in >> ComQ;

				int32_t tmp0_i32;
				if(max_flag){
					if(init_flag)
						tmp0_i32 = MIN_NEG_INT32;
					else
						tmp0_i32 = f2i32(Output[of][ofm_rc_idx]);					
					
					int32_t tmp1_i32;
					if(tmp_in_i32 > tmp0_i32)
						tmp1_i32 = tmp_in_i32;
					else
						tmp1_i32 = tmp0_i32;

					tmp_out = tmp1_i32;
				}else{
					if(init_flag)
						tmp0_i32 = 0;
					else
						tmp0_i32 = f2i32(Output[of][ofm_rc_idx]);

					tmp_out = tmp0_i32 + tmp_in_i32;
				}

				Output[of][ofm_rc_idx] = i2f32(tmp_out);
			}
		}
	}else{
		for(i =0;i < Ksize; i++)
		for(j = 0;j < Ksize; j++)
		for(tr = 0;tr < TR_MIN;tr++)
		for(tc = 0;tc < TC_MIN;tc++)
		{
			bool init_flag = ((i==0)&&(j==0));
			uint32_t t23_idx = (tr*Kstride+i)*TCol_MIN + (tc*Kstride+j);
			uint32_t ofm_rc_idx = tr*TC_MIN + tc;
			for( of = 0; of < Tn; of++)
			{
				float tmp_out, tmp0, tmp1;
				float tmp_in = Input[of][t23_idx];
				if(init_flag){
					if(max_flag)
						tmp0 = MIN_NEG;
					else
						tmp0 = 0;
				}
				else
					tmp0 = Output[of][ofm_rc_idx];
				
				if(max_flag){
					if(tmp_in > tmp0)
						tmp1 = tmp_in;
					else
						tmp1 = tmp0;

					tmp_out = tmp1;
				}else{
					tmp_out = tmp0 + tmp_in;
				}

				Output[of][ofm_rc_idx] = tmp_out;
			}
		}
	}

}

void conv2d_tile(float input_buffer[Tn][IB_HxW],float output_buffer[Tm][TrxTc],
		float weight_buffer[Tm][Tn][K][K], uint16_t n_next, uint8_t Ksize, uint8_t Kstride,uint16_t m,
		uint16_t TM_MIN, uint16_t TR_MIN, uint16_t TC_MIN, uint16_t TCol_MIN, bool enable,
		uint8_t weightQ, uint32_t ifmQ, uint8_t ComQ, bool C_SL_EN, bool interQ_en)
{
	if(!enable)
	{
		return;
	}

	uint8_t i,j, tm,tn;
	uint16_t tr,tc;
	const bool ne0 = (n_next==0);

	if(interQ_en){

		int32_t partial_mul[Tm][Tn];
		int32_t partial_add[Tm];

		float tmp_exp_w = pow(2.0, weightQ);
		float tmp_exp_ifm = pow(2.0, ifmQ);		
		for(i =0;i < Ksize; i++)
		for(j = 0;j < Ksize; j++)
		for(tr = 0;tr < TR_MIN;tr++)
		for(tc = 0;tc < TC_MIN;tc++)
		{
			uint32_t t23_idx = (Kstride*tr+i)*TCol_MIN + (Kstride*tc+j);
			uint32_t ofm_rc_idx = tr*TC_MIN + tc;
			for(tm = 0;tm < Tm;tm++)
			{
				if(i==0&&j==0&&ne0)
					partial_add[tm] = 0;
				else
					partial_add[tm] = f2i32(output_buffer[tm][ofm_rc_idx]);

				for(tn = 0;tn <Tn;tn++)
				{
					int16_t w_16i = weight_buffer[tm][tn][i][j]*tmp_exp_w;
					int16_t ifm_16i = input_buffer[tn][t23_idx]*tmp_exp_ifm;
					int32_t tmp_mul_i32 = w_16i*ifm_16i;
					int32_t tmp_sl_i32;
					if(C_SL_EN)
						tmp_sl_i32 = tmp_mul_i32 << ComQ;
					else
						tmp_sl_i32 = tmp_mul_i32 >> ComQ;

					partial_mul[tm][tn] = tmp_sl_i32;		
				}

				int32_t partial_sum = 0;
				for(tn = 0;tn <Tn;tn++)
				{
					partial_sum = partial_sum + partial_mul[tm][tn];
				}
				output_buffer[tm][ofm_rc_idx] = i2f32(partial_add[tm] + partial_sum);
			}
		}	
	}else{

		float partial_mul[Tm][Tn];
		float partial_add[Tm];

		for(i =0;i < Ksize; i++)
		for(j = 0;j < Ksize; j++)
		for(tr = 0;tr < TR_MIN;tr++)
		for(tc = 0;tc < TC_MIN;tc++)
		{
			uint32_t t23_idx = (Kstride*tr+i)*TCol_MIN + (Kstride*tc+j);
			uint32_t ofm_rc_idx = tr*TC_MIN + tc;
			for(tm = 0;tm < Tm;tm++)
			{
				if(i==0&&j==0&&ne0)
					partial_add[tm] = 0;
				else
					partial_add[tm] = output_buffer[tm][ofm_rc_idx];

				for(tn = 0;tn <Tn;tn++)
				{
					partial_mul[tm][tn] = weight_buffer[tm][tn][i][j]*input_buffer[tn][t23_idx];
				}

				float partial_sum = 0;
				for(tn = 0;tn <Tn;tn++)
				{
					partial_sum += partial_mul[tm][tn];
				}
				output_buffer[tm][ofm_rc_idx] = partial_add[tm] + partial_sum;
			}
		}		
	}
}

void ifm_weight_load_wrapper(float *ifm, float *weight, float ifm_buffer[Tn][IB_HxW], float weight_buffer[Tm][Tn][K][K], 
	int ifm_num, int ofm_num, int ifm_w, int ifm_h, bool r_init_en, bool c_init_en, uint16_t TCol_MIN_r, int pad_int,
	int TM, int TN, int tr, int tc, int tm, int tn, int ksize, int kstride, int ltype, int TM_MIN, int TR_MIN, int TC_MIN, int IHW, int INumxKK, int KK,
	bool enable, float pad_val, bool LoadBias)
{
	if(!enable)
		return;

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

	input_load(ifm, ifm_buffer, tr, tc, tmn, kstride, ksize, pad_int, TR_MIN, TC_MIN, ifm_w, ifm_h, TMN_MIN, IHW, pad_val, r_init_en, c_init_en, TCol_MIN_r);

#ifdef REORDER_GEN
	weight_load(weight,weight_buffer,tm,tn,INumxKK,KK,ksize,TM_MIN,TN_MIN, LoadBias);
#else
	weight_load_reorg(weight,weight_buffer,tm,tn,INumxKK,KK,ksize,TM_MIN,TN_MIN, LoadBias);
#endif

}

void load_compute_wrapper(float *ifm, float *weight, float ofm_buffer[Tm][TrxTc], int ksize, int kstride, int ifm_num, int ifm_w, int ifm_h,int ofm_num,
	 int ofm_h, int ofm_w, int pad_int, int ltype, int IHW, int KK, int INumxKK, int TM, int TN, int TR, int TC, int tm_r, int tr_r, int tc_r,
	 uint16_t tx_n1[3],uint16_t TX_MIN_n1[3],bool pp,bool in_flag,bool proc_flag, float pad_val, bool LoadBias, int NTif, uint8_t lmode, bool enable, 
	 bool r_init_en, bool c_init_en, uint8_t weightQ, uint32_t ifmQ, uint8_t ComQ, bool C_SL_EN, bool interQ_en)
{
	static float ifm_buffer0[Tn][IB_HxW];
	static float ifm_buffer1[Tn][IB_HxW];
	static float weight_buffer0[Tm][Tn][K][K];
	static float weight_buffer1[Tm][Tn][K][K];

	static uint16_t tx_n[2][3];//tmtrtc, between load2compute
	static uint16_t TX_MIN_n[2][3];

	if(!enable)
		return ;

	static uint16_t TR_MIN, TC_MIN, TM_MIN;
	static uint16_t TCol_MIN;
	static uint16_t TCol_MIN0[1], TCol_MIN1[1];

	if(c_init_en){
		TC_MIN = MIN_diy(TC,ofm_w-tc_r);
		TCol_MIN = (TC_MIN-1)*kstride + ksize;
	}
	if(r_init_en){
		TR_MIN = MIN_diy(TR,ofm_h-tr_r);
	}
	TM_MIN = MIN_diy(TM,ofm_num-tm_r);		

	if(lmode==0){
		if(pp){
			if(ltype == LT_CONV)
			{
				ifm_weight_load_wrapper(ifm, weight, ifm_buffer0, weight_buffer0, ifm_num, ofm_num, ifm_w, ifm_h, r_init_en, c_init_en, TCol_MIN,
						pad_int, TM, TN, tr_r, tc_r, tm_r, 0, ksize, kstride, ltype, TM_MIN, TR_MIN, TC_MIN, IHW, INumxKK, KK, in_flag, pad_val, LoadBias);
				conv2d_tile(ifm_buffer1, ofm_buffer, weight_buffer1, 0, ksize, kstride, tx_n[1][0], TX_MIN_n[1][0], TX_MIN_n[1][1], TX_MIN_n[1][2], TCol_MIN1[0], proc_flag,
						weightQ, ifmQ, ComQ, C_SL_EN, interQ_en);
			}else if(ltype == LT_MAXPOOL)
			{
				ifm_weight_load_wrapper(ifm, weight, ifm_buffer0, weight_buffer0, ifm_num, ofm_num, ifm_w, ifm_h, r_init_en, c_init_en, TCol_MIN,
						pad_int, TM, TN, tr_r, tc_r, tm_r, 0, ksize, kstride, ltype, TM_MIN, TR_MIN, TC_MIN, IHW, INumxKK, KK, in_flag, pad_val, LoadBias);							
				pool_tile(ifm_buffer1, ofm_buffer, ksize, kstride, TX_MIN_n[1][0], TX_MIN_n[1][1], TX_MIN_n[1][2], TCol_MIN1[0], proc_flag, 0,
						weightQ, ifmQ, ComQ, C_SL_EN, interQ_en);
			}								

			TCol_MIN0[0] = TCol_MIN;
			tx_n[0][0] = tm_r; tx_n[0][1] = tr_r; tx_n[0][2] = tc_r; 
			TX_MIN_n[0][0] = TM_MIN; TX_MIN_n[0][1] = TR_MIN; TX_MIN_n[0][2] = TC_MIN;

			tx_n1[0] = tx_n[1][0]; tx_n1[1] = tx_n[1][1]; tx_n1[2] = tx_n[1][2]; 
			TX_MIN_n1[0] = TX_MIN_n[1][0]; TX_MIN_n1[1] = TX_MIN_n[1][1]; TX_MIN_n1[2] = TX_MIN_n[1][2];
		}else{
			if(ltype == LT_CONV)
			{
				ifm_weight_load_wrapper(ifm, weight, ifm_buffer1, weight_buffer1, ifm_num, ofm_num, ifm_w, ifm_h, r_init_en, c_init_en, TCol_MIN,
						pad_int, TM, TN, tr_r, tc_r, tm_r, 0, ksize, kstride, ltype, TM_MIN, TR_MIN, TC_MIN, IHW, INumxKK, KK, in_flag, pad_val, LoadBias);						
				conv2d_tile(ifm_buffer0, ofm_buffer, weight_buffer0, 0, ksize, kstride, tx_n[0][0], TX_MIN_n[0][0], TX_MIN_n[0][1], TX_MIN_n[0][2], TCol_MIN0[0], proc_flag,
						weightQ, ifmQ, ComQ, C_SL_EN, interQ_en);
			}else if(ltype == LT_MAXPOOL)
			{
				ifm_weight_load_wrapper(ifm, weight, ifm_buffer1, weight_buffer1, ifm_num, ofm_num, ifm_w, ifm_h, r_init_en, c_init_en, TCol_MIN,
						pad_int, TM, TN, tr_r, tc_r, tm_r, 0, ksize, kstride, ltype, TM_MIN, TR_MIN, TC_MIN, IHW, INumxKK, KK, in_flag, pad_val, LoadBias);							
				pool_tile(ifm_buffer0, ofm_buffer, ksize, kstride, TX_MIN_n[0][0], TX_MIN_n[0][1], TX_MIN_n[0][2], TCol_MIN0[0], proc_flag, 0,
						weightQ, ifmQ, ComQ, C_SL_EN, interQ_en);
			}

			TCol_MIN1[0] = TCol_MIN;					
			tx_n[1][0] = tm_r; tx_n[1][1] = tr_r; tx_n[1][2] = tc_r; 
			TX_MIN_n[1][0] = TM_MIN; TX_MIN_n[1][1] = TR_MIN; TX_MIN_n[1][2] = TC_MIN;	

			tx_n1[0] = tx_n[0][0]; tx_n1[1] = tx_n[0][1]; tx_n1[2] = tx_n[0][2]; 
			TX_MIN_n1[0] = TX_MIN_n[0][0]; TX_MIN_n1[1] = TX_MIN_n[0][1]; TX_MIN_n1[2] = TX_MIN_n[0][2];
		}
	}else{
			if(ltype == LT_CONV)
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
						ifm_weight_load_wrapper(ifm, weight, ifm_buffer0, weight_buffer0, ifm_num, ofm_num, ifm_w, ifm_h, r_init_en, c_init_en, TCol_MIN,
								pad_int, TM, TN, tr_r, tc_r, tm_r, tn_r, ksize, kstride, ltype, TM_MIN, TR_MIN, TC_MIN, IHW, INumxKK, KK, in_flag, pad_val, LoadBias);
						conv2d_tile(ifm_buffer1, ofm_buffer, weight_buffer1, tn_n[1], ksize, kstride, tx_n[1][0], TX_MIN_n[1][0], TX_MIN_n[1][1], TX_MIN_n[1][2], TCol_MIN1[0], proc_flag,
								weightQ, ifmQ, ComQ, C_SL_EN, interQ_en);

						TCol_MIN0[0] = TCol_MIN;
						tn_n[0] = tn_r;
						tx_n[0][0] = tm_r; tx_n[0][1] = tr_r; tx_n[0][2] = tc_r;
						TX_MIN_n[0][0] = TM_MIN; TX_MIN_n[0][1] = TR_MIN; TX_MIN_n[0][2] = TC_MIN;
						tn_r += TN;

						pp_tn = 0;
					}else{
						ifm_weight_load_wrapper(ifm, weight, ifm_buffer1, weight_buffer1, ifm_num, ofm_num, ifm_w, ifm_h, r_init_en, c_init_en, TCol_MIN,
								pad_int, TM, TN, tr_r, tc_r, tm_r, tn_r, ksize, kstride, ltype, TM_MIN, TR_MIN, TC_MIN, IHW, INumxKK, KK, in_flag, pad_val, LoadBias);
						conv2d_tile(ifm_buffer0, ofm_buffer, weight_buffer0, tn_n[0], ksize, kstride, tx_n[0][0], TX_MIN_n[0][0], TX_MIN_n[0][1], TX_MIN_n[0][2], TCol_MIN0[0], proc_flag,
								weightQ, ifmQ, ComQ, C_SL_EN, interQ_en);

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

void FPGA_Acc(float *ifm, float *ofm, float *weight, float* bias, uint32_t k_s_pad_ltype, uint32_t iofm_num, uint32_t ifm_w_h, uint32_t ofm_w_h,
	uint32_t TRTC, uint32_t TMTN, int32_t NToy, int32_t NTox, int32_t NTof, int32_t NTcomb, int32_t NTif, uint8_t lmode, int32_t NTcomb_l, float pad_val, float div_kk,
	uint32_t IHW, uint32_t OHW, uint32_t KK_INumxKK, uint32_t en_bits, 
	uint32_t weightQ, uint32_t biasQ, uint32_t ifmQ, uint32_t ofmQ, uint32_t avgQ, uint32_t interQ, bool interQ_en)//enable_bits[2:0]={IsReLU, LoadBias, IsNotConv}
{
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
	uint8_t KK = (KK_INumxKK >> 24) & 0xFF;
	uint32_t INumxKK = KK_INumxKK & 0xFFFFFF;//24bit
	bool IsReLU    = (en_bits >> 2) & 0x1;
	bool LoadBias  = (en_bits >> 1) & 0x1;
	bool IsNotConv =  en_bits       & 0x1;
	bool IsAVGPOOL = (ltype == LT_AVGPOOL);

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
	assert((TR*TC > 0)&&(TR*TC <= Tr*Tc));

	static float ofm_buffer0[Tm][TrxTc];
	static float ofm_buffer1[Tm][TrxTc];
	static float bias_buffer[MAX_BETA_LENGTH];

	uint16_t tx_n10[3], tx_n11[3];//tmtrtc, between compute2store
	uint16_t TX_MIN_n10[3], TX_MIN_n11[3];

	if(LoadBias)
		memcpy(bias_buffer, bias, ofm_num*sizeof(float));

	uint16_t tr, tc, tm;
	uint16_t TR_MIN, TC_MIN, TM_MIN;
	uint16_t tm_r, tr_r, tc_r;
	bool in_flag, proc_flag, out_flag, lc_enable;

	uint8_t ComQ = 0;
	bool C_SL_EN = 1;
	if(interQ_en){
		if((ltype == LT_CONV) || (ltype == LT_DCONV))
		{
			int8_t ComQ_l = interQ - (ifmQ + weightQ);
			if(ComQ_l < 0){
				C_SL_EN = 0;
				ComQ = -ComQ_l;
			}else{
				C_SL_EN = 1;
				ComQ = ComQ_l;
			}
		}else if((ltype == LT_AVGPOOL) || (ltype == LT_MAXPOOL))
		{
			int8_t ComQ_l = interQ - ifmQ;
			if(ComQ_l < 0){
				C_SL_EN = 0;
				ComQ = -ComQ_l;
			}else{
				C_SL_EN = 1;
				ComQ = ComQ_l;
			}
		}
		assert((ComQ >=0)&&(ComQ <32));	
		// printf("%d, %3.7lf, %3.7lf\n", ofmQ, pow(2.0, -ofmQ), pow(1/2.0,ofmQ));//maybe g++ compiler or lib's bug
	}

	
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
			load_compute_wrapper(ifm, weight, ofm_buffer0, ksize, kstride, ifm_num, ifm_w, ifm_h, ofm_num, ofm_h, ofm_w,
				pad_int, ltype, IHW, KK, INumxKK, TM, TN, TR, TC, tm_r, tr_r, tc_r,
				tx_n10, TX_MIN_n10, pp, in_flag, proc_flag, pad_val, LoadBias, NTif, lmode, lc_enable, r_init_en, c_init_en, weightQ, ifmQ, ComQ, C_SL_EN, interQ_en);

			write_back_output_reorg( ofm_buffer1, bias_buffer, ofm, tx_n11[1], tx_n11[2], tx_n11[0], ofm_num, ofm_h, ofm_w, TX_MIN_n11[0], TX_MIN_n11[1], TX_MIN_n11[2], OHW, IsReLU, div_kk,
				IsAVGPOOL, LoadBias, out_flag, biasQ, ofmQ, avgQ, interQ, interQ_en);

			pp = 0;
		}else{
			load_compute_wrapper(ifm, weight, ofm_buffer1, ksize, kstride, ifm_num, ifm_w, ifm_h, ofm_num, ofm_h, ofm_w,
				pad_int, ltype, IHW, KK, INumxKK, TM, TN, TR, TC, tm_r, tr_r, tc_r,
				tx_n11, TX_MIN_n11, pp, in_flag, proc_flag, pad_val, LoadBias, NTif, lmode, lc_enable, r_init_en, c_init_en, weightQ, ifmQ, ComQ, C_SL_EN, interQ_en);			

			write_back_output_reorg( ofm_buffer0, bias_buffer, ofm, tx_n10[1], tx_n10[2], tx_n10[0], ofm_num, ofm_h, ofm_w, TX_MIN_n10[0], TX_MIN_n10[1], TX_MIN_n10[2], OHW, IsReLU, div_kk,
				IsAVGPOOL, LoadBias, out_flag, biasQ, ofmQ, avgQ, interQ, interQ_en);

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
