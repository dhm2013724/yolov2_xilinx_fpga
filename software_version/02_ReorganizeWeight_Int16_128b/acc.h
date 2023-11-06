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
#define LANE_NUM 2
// #define INTER_SF 12
#define EXTRA_BIT 10
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
const uint32_t Tmax_dx = (Tm+LANE_NUM-1)/LANE_NUM;
const uint32_t Tnax_dx = (Tn+LANE_NUM-1)/LANE_NUM;

void input_load(int16_t *input,int16_t input_buffer[Tnax_dx][IB_HxW][LANE_NUM], uint16_t r, uint16_t c, uint16_t n, uint8_t Kstride, uint8_t ksize, uint8_t Padding,
		 uint16_t TR_MIN, uint16_t TC_MIN,
		 uint16_t Input_w, uint16_t Input_h, uint8_t TN_MIN, uint32_t IHxIW, int16_t pad_val, bool r_init_en, bool c_init_en, uint16_t TCol_MIN_r)
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
	uint16_t TN_MINax_dx = (TN_MIN + LANE_NUM -1)/LANE_NUM;
	uint16_t nax_dx = n/LANE_NUM;

	uint32_t burstlen;
	uint32_t offset_base;
	uint16_t T2_MIN;
	uint32_t T12_MIN;
	if(IsAllCont){
		burstlen = TN_MINax_dx*row_len*col_len;
		offset_base = nax_dx*IHxIW;
		T2_MIN = 1;
		T12_MIN = 1;
	}else if(IsColCont){
		burstlen = row_len*col_len;
		offset_base = nax_dx*IHxIW + tl_y*Input_w;
		T2_MIN = 1;
		T12_MIN = TN_MINax_dx;
	}else{
		burstlen = col_len;
		offset_base = nax_dx*IHxIW + tl_y*Input_w + tl_x;
		T2_MIN = row_len;
		T12_MIN = TN_MINax_dx*row_len;		
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

		int16_t *ifm_addr = input + offset_ex*LANE_NUM;
		for(uint32_t tbl=0; tbl<burstlen; tbl++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=TnxIB_HxIB_W)
#pragma HLS PIPELINE II=1
			uint16_t t2_idx = t2+px_t;
			uint16_t t3_idx = t3+px_l;
			// input_buffer[t1][t2+px_t][t3+px_l] = ifm_addr[tbl];
			// input_buffer[t1][t2_idx*TCol_MIN_r + t3_idx] = ifm_addr[tbl];

			int16_t tmp_i16c[LANE_NUM];
			memcpy(tmp_i16c, ifm_addr + (tbl*LANE_NUM), sizeof(int16_t)*LANE_NUM);
			for(uint32_t x=0; x<LANE_NUM; x++){
				input_buffer[t1][t2_idx*TCol_MIN_r + t3_idx][x] = tmp_i16c[x];
			}	

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
					input_buffer[t1/LANE_NUM][t23_idx][t1 & (LANE_NUM-1)] = pad_val;
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
					input_buffer[t1/LANE_NUM][t23_idx][t1 & (LANE_NUM - 1)] = pad_val;
				}
			}
		}}
	}		
	
// // if weight has set 0 within idle channel of N dim, we dont need to clear values in idel channels 	
// 	uint16_t TN_left = Tn - TN_MIN;
// 	if(TN_left){
// 		for(uint16_t t1 = 0;t1 < TN_left; t1++){
// DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tn)				
// 		for(uint16_t t2 = 0;t2 < TRow_MIN; t2++){
// DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=IB_H)
// 		for(uint16_t t3 = 0;t3 < TCol_MIN_r; t3++){
// DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=IB_W)
// #pragma HLS PIPELINE II=1
// 				input_buffer[(t1+TN_MIN)/LANE_NUM][t2*TCol_MIN_r+t3][(t1+TN_MIN) & (LANE_NUM-1)] = 0;
// 		}}}			
// 	}

}

void weight_load_reorg(int16_t *Weight,int16_t weight_buffer[Tmax_dx][Tn][K][K][LANE_NUM],uint16_t m,uint16_t n, uint32_t IFM_numxKxK, uint8_t KxK, uint8_t Ksize, 
			uint8_t TM_MIN, uint8_t TN_MIN, bool enable)
{
const uint32_t TmxTnxKxK = ((Tm+LANE_NUM-1)/LANE_NUM)*Tn*K*K*LANE_NUM;	

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
		for(t1 = 0;t1 < Tmax_dx; t1++){		
		for(t2 = 0;t2 < Tn; t2++){
			for(uint32_t x =0; x<LANE_NUM; x++){
				weight_buffer[t1][t2][t3][t4][x] = 0;
			}
		}}
	}}	

	if(m==0&&n==0)
		Woffset = 0;

	uint32_t burstlen = ((TM_MIN+LANE_NUM-1)/LANE_NUM)*TN_MIN*KxK;
	assert(burstlen <= TmxTnxKxK);

	int16_t *weight_addr = Weight + Woffset*LANE_NUM;
	Woffset += burstlen;
	t1 = 0; t2 = 0; t3 = 0; t4 = 0;
	for(uint32_t tbl = 0; tbl < burstlen; tbl++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=TmxTnxKxK)
#pragma HLS PIPELINE II=1

		// weight_buffer[t1][t2][t3][t4] = weight_addr[tbl];
		int16_t tmp_i16c[LANE_NUM];
		memcpy(tmp_i16c, weight_addr + (tbl*LANE_NUM), sizeof(int16_t)*LANE_NUM);
		for(uint32_t x=0; x<LANE_NUM; x++){
			weight_buffer[t1][t2][t3][t4][x] = tmp_i16c[x];
		}

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

	// uint8_t avgQ = 15;
	// int16_t div_kk_i16 = 0x29c;//nearly === 1.0f/KK, K=7
	if(IsAVGPOOL)
		tmp1 = (tmp0*div_kk) >> avgQ;
	else
		tmp1 = tmp0;

	const int32_t const_ISF_6 = (6 << interQ);
	if(IsNL)
	{
#ifdef RELU6
		tmp2 = MIN_diy(MAX_diy(0, tmp1), const_ISF_6);
		// tmp_out = MIN_diy(MAX_diy(0.0f,tmp1), 6.0f);
#elif LEAKY_RELU
		// if(tmp1 < 0.0f)
		// 	tmp_out = tmp1*0.1f;
		// else
		// 	tmp_out = tmp1;
		if(tmp1 < 0)
			tmp2 = (tmp1*0xccc) >> 15;
		else
			tmp2 = tmp1;
#else
		tmp2 = MAX_diy(0,tmp1);
#endif
	}
	else
		tmp2 = tmp1;

	if(ofmQ > interQ)
		tmp_out = tmp2 << (ofmQ - interQ);
	else
		tmp_out = tmp2 >> (interQ - ofmQ);
	
	return tmp_out;
}

void copy_local_beta(int16_t *bias_buffer,int16_t local_beta_buffer[Tmax_dx][LANE_NUM], uint16_t TM_MIN, uint16_t m)
{	
	for(int tm = 0;tm < TM_MIN;tm++)
	{
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tm)
#pragma HLS PIPELINE II=1
		local_beta_buffer[tm/LANE_NUM][tm & (LANE_NUM-1)] = bias_buffer[m+tm];
	}	
}

void write_back_output_reorg(int32_t output_buffer[Tmax_dx][TrxTc][LANE_NUM], int16_t *bias_buffer, /*float bias_buffer[MAX_BETA_LENGTH],*/
		int16_t *Output, uint16_t r,uint16_t c,uint16_t m, 
		uint16_t ofm_num, uint16_t ofm_h, uint16_t ofm_w,
		uint16_t TM_MIN, uint16_t TR_MIN, uint16_t TC_MIN, uint32_t OHxOW, bool IsNL, int16_t div_kk, bool IsAVGPOOL, bool LoadBias, bool enable, 
		uint8_t biasQ, uint8_t ofmQ, uint8_t avgQ, uint8_t interQ)
{
	const uint32_t OFM_BLmax = 256;

	static int16_t local_bias_buffer[Tmax_dx][LANE_NUM];

	if(!enable)
		return;

	assert((TM_MIN >0)&&(TM_MIN <=Tm));
	// assert((TR_MIN >0)&&(TR_MIN <=Tr));
	// assert((TC_MIN >0)&&(TC_MIN <=Tc));
	assert((TR_MIN*TC_MIN > 0)&&(TR_MIN*TC_MIN <= Tr*Tc));

	copy_local_beta( bias_buffer, local_bias_buffer, TM_MIN, m);

	uint32_t offset = (m/LANE_NUM)*OHxOW + r*ofm_w + c;
	bool IsColCont = (ofm_w==TC_MIN);//equal FMCont
	bool IsRowCont = (ofm_h==TR_MIN);
	bool IsAllCont = IsColCont && IsRowCont;
	uint16_t TM_MINax_dx = (TM_MIN + LANE_NUM -1)/LANE_NUM;

	uint32_t burstlen;
	uint16_t T1_MIN, T2_MIN;
	if(IsAllCont){
		// burstlen = TM_MIN*TR_MIN*TC_MIN;
		burstlen = TM_MINax_dx*TR_MIN*ofm_w;
		T1_MIN = 1;
		T2_MIN = 1;
	}else if(IsColCont){
		burstlen = TR_MIN*ofm_w;
		T1_MIN = TM_MINax_dx;
		T2_MIN = 1;
	}else{
		burstlen = TC_MIN;
		T1_MIN = TM_MINax_dx;
		T2_MIN = TR_MIN;
	}

	// assert(burstlen > 0);
	// uint32_t burstlen_a_blmax = burstlen + OFM_BLmax;
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

		int16_t *OFM_base = Output + ofm_offset*LANE_NUM;
		uint16_t TBL_MIN = MIN_diy(OFM_BLmax, burstlen-tbl_r);

		for(int tbl_min=0; tbl_min < TBL_MIN; tbl_min++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=OFM_BLmax)
#pragma HLS PIPELINE II=1
			int16_t tmp_i16c[LANE_NUM];
			for(uint32_t x=0; x<LANE_NUM; x++){
				tmp_i16c[x] = postproc(((output_buffer[tm][tr*TC_MIN + tc][x] << EXTRA_BIT) >> EXTRA_BIT), local_bias_buffer[tm][x],
									div_kk, LoadBias, IsAVGPOOL, IsNL, biasQ, ofmQ, avgQ, interQ);
			}
			memcpy(OFM_base + (tbl_min*LANE_NUM), tmp_i16c, sizeof(int16_t)*LANE_NUM);
			// OFM_base[tbl_min] = postproc(output_buffer[tm][tr*TC_MIN + tc], local_bias_buffer[tm], div_kk, LoadBias, IsAVGPOOL, IsNL);

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

// void avgpool_tile(float Input[Tn][OnChipIB_Height][OnChipIB_Width],float Output[Tm][Tr][Tc],
// 		  uint8_t Ksize, uint8_t Kstride, uint8_t TM_MIN, uint8_t TR_MIN, uint8_t TC_MIN, bool enable)
// {
// 	if(!enable)
// 		return;

// 	uint8_t i,j,tr,tc,of;
// 	float tmp[Tn];
// #pragma HLS ARRAY_PARTITION variable=tmp complete dim=1

// 	for(i =0;i < Ksize; i++)
// DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=K)
// 	for(j = 0;j < Ksize; j++)
// DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=K)
// 	for(tr = 0;tr < TR_MIN;tr++)
// DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tr)
// 		for(tc = 0;tc < TC_MIN;tc++)
// 		{
// DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tc)
// #pragma HLS PIPELINE II=1
// 			for( of = 0; of < Tn; of++)
// 			{
// 				float tmp_in = Input[of][tr*Kstride+i][tc*Kstride+j];
// 				if(i==0&&j==0)
// 					tmp[of] = 0;
// 				else
// 					tmp[of] = Output[of][tr][tc];

// 				Output[of][tr][tc] = tmp[of] + tmp_in;
// 			}
// 		}
// }

void pool_tile(int16_t Input[Tnax_dx][IB_HxW][LANE_NUM],int32_t Output[Tmax_dx][TrxTc][LANE_NUM],
	uint8_t Ksize, uint8_t Kstride, uint8_t TM_MIN, uint16_t TR_MIN, uint16_t TC_MIN, uint16_t TCol_MIN, bool enable, uint8_t pool_mode, uint8_t ComQ, bool C_SL_EN)//0:MAX; 1:AVG
{
	if(!enable)
		return;

	uint8_t i,j,of;
	uint16_t tr,tc;

	bool max_flag = (pool_mode==0);

	// float tmp_exp = pow(2.0, -ifmQ);
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
			bool init_flag = ((i==0)&&(j==0));
			uint32_t t23_idx = (tr*Kstride+i)*TCol_MIN + (tc*Kstride+j);
			uint32_t ofm_rc_idx = tr*TC_MIN + tc;
			for( of = 0; of < Tn; of++)
			{
				int32_t tmp_out;
				int16_t tmp_in = Input[of/LANE_NUM][t23_idx][of & (LANE_NUM-1)];
				// int32_t tmp_in_i32 = (tmp_in >> (ifmQ - INTER_SF));
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
						tmp0_i32 = Output[of/LANE_NUM][ofm_rc_idx][of & (LANE_NUM-1)];					
					
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
						tmp0_i32 = Output[of/LANE_NUM][ofm_rc_idx][of & (LANE_NUM-1)];

					tmp_out = tmp0_i32 + tmp_in_i32;
				}

				Output[of/LANE_NUM][ofm_rc_idx][of & (LANE_NUM-1)] = ((tmp_out << EXTRA_BIT) >> EXTRA_BIT);
			}
		}
}

// void pool_tile(int16_t Input[Tnax_dx][IB_HxW][LANE_NUM],int32_t Output[Tmax_dx][TrxTc][LANE_NUM],
// 	uint8_t Ksize, uint8_t Kstride, uint8_t TM_MIN, uint16_t TR_MIN, uint16_t TC_MIN, uint16_t TCol_MIN, bool enable, uint8_t pool_mode, uint8_t ifmQ)//0:MAX; 1:AVG
// {
// 	if(!enable)
// 		return;

// 	uint8_t i,j,of;
// 	uint16_t tr,tc;

// 	bool max_flag = (pool_mode==0);

// 	// float tmp_exp = pow(2.0, -ifmQ);
// 	for(i =0;i < Ksize; i++)
// DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=K)
// 	for(j = 0;j < Ksize; j++)
// DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=K)
// 	for(tr = 0;tr < TR_MIN;tr++)
// DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tr)
// 		for(tc = 0;tc < TC_MIN;tc++)
// 		{
// DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tc)
// #pragma HLS PIPELINE II=1
// 			bool init_flag = ((i==0)&&(j==0));
// 			uint32_t t23_idx = (tr*Kstride+i)*TCol_MIN + (tc*Kstride+j);
// 			uint32_t ofm_rc_idx = tr*TC_MIN + tc;
// 			for( of = 0; of < Tn; of++)
// 			{
// 				int32_t tmp_out;
// 				int16_t tmp_in = Input[of/LANE_NUM][t23_idx][of & (LANE_NUM-1)];
// 				if(max_flag){
// 					int32_t tmp0_i32;
// 					if(init_flag)
// 						tmp0_i32 = MIN_NEG_INT32;
// 					else
// 						tmp0_i32 = Output[of/LANE_NUM][ofm_rc_idx][of & (LANE_NUM-1)];					
					
// 					int32_t tmp1_i32;
// 					int32_t tmp2_i32 = (tmp_in >> (ifmQ - INTER_SF));
// 					if(tmp2_i32 > tmp0_i32)
// 						tmp1_i32 = tmp2_i32;
// 					else
// 						tmp1_i32 = tmp0_i32;

// 					tmp_out = tmp1_i32;
// 				}else{
// 					int32_t tmp0_i32;
// 					if(init_flag)
// 						tmp0_i32 = 0;
// 					else
// 						tmp0_i32 = Output[of/LANE_NUM][ofm_rc_idx][of & (LANE_NUM-1)];

// 					tmp_out = tmp0_i32 + (tmp_in >> (ifmQ - INTER_SF));
// 				}

// 				Output[of/LANE_NUM][ofm_rc_idx][of & (LANE_NUM-1)] = ((tmp_out << EXTRA_BIT) >> EXTRA_BIT);
// 			}
// 		}
// }

//weight[MxNxKxK] bias[M] ifm[NxRxC] 
void dconv2d_tile(int16_t input_buffer[Tnax_dx][IB_HxW][LANE_NUM], int32_t output_buffer[Tmax_dx][TrxTc][LANE_NUM],
		  int16_t weight_buffer[Tmax_dx][Tn][K][K][LANE_NUM], uint8_t Ksize, uint8_t Kstride, uint16_t m, uint16_t TM_MIN, uint16_t TR_MIN, 
		  uint16_t TC_MIN, uint16_t TCol_MIN, bool enable, uint8_t ComQ, bool C_SL_EN)
{

	if(!enable)
		return;

	assert(TM_MIN <= Tn);

	uint8_t i,j,tm;
	uint16_t tr,tc;
	int32_t partial_mul[Tn];
#pragma HLS ARRAY_PARTITION variable=partial_mul complete dim=1
	int32_t partial_add[Tn];
#pragma HLS ARRAY_PARTITION variable=partial_add complete dim=1

	// float tmp_exp = pow(2.0, -weightQ);
	// float tmp_exp0 = pow(2.0, -ifmQ);
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
					uint32_t t23_idx = (Kstride*tr+i)*TCol_MIN + (Kstride*tc+j);
					uint32_t ofm_rc_idx = tr*TC_MIN + tc;
					for(tm = 0;tm < Tn;tm++)
					{
#pragma HLS DEPENDENCE variable=output_buffer inter false
						if(i==0&&j==0)
							partial_add[tm] = 0;
						else
							partial_add[tm] = output_buffer[tm/LANE_NUM][ofm_rc_idx][tm & (LANE_NUM-1)];

						// partial_mul[tm] = (weight_buffer[tm/LANE_NUM][0][i][j][tm & (LANE_NUM-1)]*tmp_exp)*(input_buffer[tm/LANE_NUM][t23_idx][tm & (LANE_NUM-1)]*tmp_exp0);
							partial_mul[tm] = weight_buffer[tm/LANE_NUM][0][i][j][tm & (LANE_NUM-1)]*input_buffer[tm/LANE_NUM][t23_idx][tm & (LANE_NUM-1)];
						
							int32_t tmp_sl_i32;
							if(C_SL_EN)
								tmp_sl_i32 = partial_mul[tm] << ComQ;
							else
								tmp_sl_i32 = partial_mul[tm] >> ComQ;
							
							int32_t tmp_out = partial_add[tm] + tmp_sl_i32;

						// int32_t tmp_out = partial_add[tm] + (partial_mul[tm] >> (ifmQ + weightQ - INTER_SF));
						output_buffer[tm/LANE_NUM][ofm_rc_idx][tm & (LANE_NUM-1)] = ((tmp_out << EXTRA_BIT) >> EXTRA_BIT);
					}

				}
}

void conv2d_tile(int16_t input_buffer[Tnax_dx][IB_HxW][LANE_NUM],int32_t output_buffer[Tmax_dx][TrxTc][LANE_NUM],
		int16_t weight_buffer[Tmax_dx][Tn][K][K][LANE_NUM], uint16_t n_next, uint8_t Ksize, uint8_t Kstride,uint16_t m,
		uint16_t TM_MIN, uint16_t TR_MIN, uint16_t TC_MIN, uint16_t TCol_MIN, bool enable, uint8_t ComQ, bool C_SL_EN)
{
	if(!enable)
	{
		return;
	}

	uint8_t i,j, tm,tn;
	uint16_t tr,tc;
	const bool ne0 = (n_next==0);

	int32_t partial_mul[Tm][Tn];
#pragma HLS ARRAY_PARTITION variable=partial_mul complete dim=1
#pragma HLS ARRAY_PARTITION variable=partial_mul complete dim=2
	int32_t partial_add[Tm];
#pragma HLS ARRAY_PARTITION variable=partial_add complete dim=1

	// float tmp_exp = pow(2.0, -weightQ);
	// float tmp_exp0 = pow(2.0, -ifmQ);
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
					uint32_t t23_idx = (Kstride*tr+i)*TCol_MIN + (Kstride*tc+j);
					uint32_t ofm_rc_idx = tr*TC_MIN + tc;
					for(tm = 0;tm < Tm;tm++)
					{
#pragma HLS DEPENDENCE variable=output_buffer inter false

						if(i==0&&j==0&&ne0)
							partial_add[tm] = 0;
						else
							partial_add[tm] = output_buffer[tm/LANE_NUM][ofm_rc_idx][tm & (LANE_NUM -1)];

						for(tn = 0;tn <Tn;tn++)
						{
							// int32_t tmp_mul_i32 = (weight_buffer[tm/LANE_NUM][tn][i][j][tm &(LANE_NUM-1)]*tmp_exp)*(input_buffer[tn/LANE_NUM][t23_idx][tn & (LANE_NUM-1)]*tmp_exp0);
							int32_t tmp_mul_i32 = weight_buffer[tm/LANE_NUM][tn][i][j][tm &(LANE_NUM-1)]*input_buffer[tn/LANE_NUM][t23_idx][tn & (LANE_NUM-1)];
							// partial_mul[tm][tn] = tmp_mul_i32 >> (ifmQ + weightQ - INTER_SF);

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
							 partial_sum += partial_mul[tm][tn];
						}
						int32_t tmp_out = partial_add[tm] + partial_sum;
						output_buffer[tm/LANE_NUM][ofm_rc_idx][tm & (LANE_NUM-1)] = ((tmp_out << EXTRA_BIT) >> EXTRA_BIT);
					}

				}
}

void ifm_weight_load_wrapper(int16_t *ifm, int16_t *weight, int16_t ifm_buffer[Tnax_dx][IB_HxW][LANE_NUM], int16_t weight_buffer[Tmax_dx][Tn][K][K][LANE_NUM], 
	int ifm_num, int ofm_num, int ifm_w, int ifm_h, bool r_init_en, bool c_init_en, uint16_t TCol_MIN_r, int pad_int,
	int TM, int TN, int tr, int tc, int tm, int tn, int ksize, int kstride, int ltype, int TM_MIN, int TR_MIN, int TC_MIN, int IHW, int INumxKK, int KK,
	bool enable, int16_t pad_val, bool LoadBias)
{
	if(!enable)
		return;

	int TN_MIN, TMN_MIN;
	int tmn;
	// bool IsPad0ExtraFM;
	if(ltype == LT_CONV)
	{
		tmn = tn;
		TN_MIN = MIN_diy(TN, ifm_num-tn);
		TMN_MIN = TN_MIN;
		// IsPad0ExtraFM = 1;
	}
	else
	{
		tmn = tm;
		TN_MIN = 1;//DCONV
		TMN_MIN = TM_MIN;
		// IsPad0ExtraFM = 0;
	}

	input_load(ifm, ifm_buffer, tr, tc, tmn, kstride, ksize, pad_int, TR_MIN, TC_MIN, ifm_w, ifm_h, TMN_MIN, IHW, pad_val, r_init_en, c_init_en, TCol_MIN_r);

	weight_load_reorg(weight,weight_buffer,tm,tn,INumxKK,KK,ksize,TM_MIN,TN_MIN, LoadBias);
}

void load_compute_wrapper(int16_t *ifm, int16_t *weight, int32_t ofm_buffer[Tmax_dx][TrxTc][LANE_NUM], int ksize, int kstride, int ifm_num, int ifm_w, int ifm_h,int ofm_num,
	 int ofm_h, int ofm_w, int pad_int, int ltype, int IHW, int KK, int INumxKK, int TM, int TN, int TR, int TC, int tm_r, int tr_r, int tc_r,
	 uint16_t tx_n1[3],uint16_t TX_MIN_n1[3],bool pp,bool in_flag,bool proc_flag, int16_t pad_val, bool LoadBias, int NTif, uint8_t lmode, bool enable, 
	 bool r_init_en, bool c_init_en, uint8_t ComQ, bool C_SL_EN)
{
	static int16_t ifm_buffer0[Tnax_dx][IB_HxW][LANE_NUM];
// #pragma HLS ARRAY_PARTITION variable=ifm_buffer0 complete dim=1
	static int16_t ifm_buffer1[Tnax_dx][IB_HxW][LANE_NUM];
// #pragma HLS ARRAY_PARTITION variable=ifm_buffer1 complete dim=1
	static int16_t weight_buffer0[Tmax_dx][Tn][K][K][LANE_NUM];
// #pragma HLS ARRAY_PARTITION variable=weight_buffer0 complete dim=1
// #pragma HLS ARRAY_PARTITION variable=weight_buffer0 complete dim=2
	static int16_t weight_buffer1[Tmax_dx][Tn][K][K][LANE_NUM];
// #pragma HLS ARRAY_PARTITION variable=weight_buffer1 complete dim=1
// #pragma HLS ARRAY_PARTITION variable=weight_buffer1 complete dim=2

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
				conv2d_tile(ifm_buffer1, ofm_buffer, weight_buffer1, 0, ksize, kstride, tx_n[1][0], TX_MIN_n[1][0], TX_MIN_n[1][1], TX_MIN_n[1][2], TCol_MIN1[0], proc_flag, ComQ, C_SL_EN);
			}else if(ltype == LT_DCONV)
			{
				ifm_weight_load_wrapper(ifm, weight, ifm_buffer0, weight_buffer0, ifm_num, ofm_num, ifm_w, ifm_h, r_init_en, c_init_en, TCol_MIN,
						pad_int, TM, TN, tr_r, tc_r, tm_r, 0, ksize, kstride, ltype, TM_MIN, TR_MIN, TC_MIN, IHW, KK, KK, in_flag, pad_val, LoadBias);						
				dconv2d_tile(ifm_buffer1, ofm_buffer, weight_buffer1, ksize, kstride, tx_n[1][0], TX_MIN_n[1][0], TX_MIN_n[1][1], TX_MIN_n[1][2], TCol_MIN1[0], proc_flag, ComQ, C_SL_EN);
			}else if(ltype == LT_AVGPOOL)
			{
				ifm_weight_load_wrapper(ifm, weight, ifm_buffer0, weight_buffer0, ifm_num, ofm_num, ifm_w, ifm_h, r_init_en, c_init_en, TCol_MIN,
						pad_int, TM, TN, tr_r, tc_r, tm_r, 0, ksize, kstride, ltype, TM_MIN, TR_MIN, TC_MIN, IHW, INumxKK, KK, in_flag, pad_val, LoadBias);							
				pool_tile(ifm_buffer1, ofm_buffer, ksize, kstride, TX_MIN_n[1][0], TX_MIN_n[1][1], TX_MIN_n[1][2], TCol_MIN1[0], proc_flag, 1, ComQ, C_SL_EN);
			}else if(ltype == LT_MAXPOOL)
			{
				ifm_weight_load_wrapper(ifm, weight, ifm_buffer0, weight_buffer0, ifm_num, ofm_num, ifm_w, ifm_h, r_init_en, c_init_en, TCol_MIN,
						pad_int, TM, TN, tr_r, tc_r, tm_r, 0, ksize, kstride, ltype, TM_MIN, TR_MIN, TC_MIN, IHW, INumxKK, KK, in_flag, pad_val, LoadBias);							
				pool_tile(ifm_buffer1, ofm_buffer, ksize, kstride, TX_MIN_n[1][0], TX_MIN_n[1][1], TX_MIN_n[1][2], TCol_MIN1[0], proc_flag, 0, ComQ, C_SL_EN);
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
				conv2d_tile(ifm_buffer0, ofm_buffer, weight_buffer0, 0, ksize, kstride, tx_n[0][0], TX_MIN_n[0][0], TX_MIN_n[0][1], TX_MIN_n[0][2], TCol_MIN0[0], proc_flag, ComQ, C_SL_EN);
			}else if(ltype == LT_DCONV)
			{
				ifm_weight_load_wrapper(ifm, weight, ifm_buffer1, weight_buffer1, ifm_num, ofm_num, ifm_w, ifm_h, r_init_en, c_init_en, TCol_MIN,
						pad_int, TM, TN, tr_r, tc_r, tm_r, 0, ksize, kstride, ltype, TM_MIN, TR_MIN, TC_MIN, IHW, KK, KK, in_flag, pad_val, LoadBias);						
				dconv2d_tile(ifm_buffer0, ofm_buffer, weight_buffer0, ksize, kstride, tx_n[0][0], TX_MIN_n[0][0], TX_MIN_n[0][1], TX_MIN_n[0][2], TCol_MIN0[0], proc_flag, ComQ, C_SL_EN);
			}else if(ltype == LT_AVGPOOL)
			{
				ifm_weight_load_wrapper(ifm, weight, ifm_buffer1, weight_buffer1, ifm_num, ofm_num, ifm_w, ifm_h, r_init_en, c_init_en, TCol_MIN,
						pad_int, TM, TN, tr_r, tc_r, tm_r, 0, ksize, kstride, ltype, TM_MIN, TR_MIN, TC_MIN, IHW, INumxKK, KK, in_flag, pad_val, LoadBias);							
				pool_tile(ifm_buffer0, ofm_buffer, ksize, kstride, TX_MIN_n[0][0], TX_MIN_n[0][1], TX_MIN_n[0][2], TCol_MIN0[0], proc_flag, 1, ComQ, C_SL_EN);
			}else if(ltype == LT_MAXPOOL)
			{
				ifm_weight_load_wrapper(ifm, weight, ifm_buffer1, weight_buffer1, ifm_num, ofm_num, ifm_w, ifm_h, r_init_en, c_init_en, TCol_MIN,
						pad_int, TM, TN, tr_r, tc_r, tm_r, 0, ksize, kstride, ltype, TM_MIN, TR_MIN, TC_MIN, IHW, INumxKK, KK, in_flag, pad_val, LoadBias);							
				pool_tile(ifm_buffer0, ofm_buffer, ksize, kstride, TX_MIN_n[0][0], TX_MIN_n[0][1], TX_MIN_n[0][2], TCol_MIN0[0], proc_flag, 0, ComQ, C_SL_EN);
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
						conv2d_tile(ifm_buffer1, ofm_buffer, weight_buffer1, tn_n[1], ksize, kstride, tx_n[1][0], TX_MIN_n[1][0], TX_MIN_n[1][1], TX_MIN_n[1][2], TCol_MIN1[0], proc_flag, ComQ, C_SL_EN);

						TCol_MIN0[0] = TCol_MIN;
						tn_n[0] = tn_r;
						tx_n[0][0] = tm_r; tx_n[0][1] = tr_r; tx_n[0][2] = tc_r;
						TX_MIN_n[0][0] = TM_MIN; TX_MIN_n[0][1] = TR_MIN; TX_MIN_n[0][2] = TC_MIN;
						tn_r += TN;

						pp_tn = 0;
					}else{
						ifm_weight_load_wrapper(ifm, weight, ifm_buffer1, weight_buffer1, ifm_num, ofm_num, ifm_w, ifm_h, r_init_en, c_init_en, TCol_MIN,
								pad_int, TM, TN, tr_r, tc_r, tm_r, tn_r, ksize, kstride, ltype, TM_MIN, TR_MIN, TC_MIN, IHW, INumxKK, KK, in_flag, pad_val, LoadBias);
						conv2d_tile(ifm_buffer0, ofm_buffer, weight_buffer0, tn_n[0], ksize, kstride, tx_n[0][0], TX_MIN_n[0][0], TX_MIN_n[0][1], TX_MIN_n[0][2], TCol_MIN0[0], proc_flag, ComQ, C_SL_EN);

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

void FPGA_Acc(int16_t *ifm, int16_t *ofm, int16_t *weight, int16_t* bias, uint32_t k_s_pad_ltype, uint32_t iofm_num, uint32_t ifm_w_h, uint32_t ofm_w_h,
	uint32_t TRTC, uint32_t TMTN, int32_t NToy, int32_t NTox, int32_t NTof, int32_t NTcomb, int32_t NTif, uint8_t lmode, int32_t NTcomb_l, int16_t pad_val, int16_t div_kk,
	uint32_t IHW, uint32_t OHW, uint32_t KK_INumxKK, uint32_t en_bits, uint32_t weightQ, uint32_t biasQ, uint32_t ifmQ, uint32_t ofmQ, uint32_t avgQ, uint32_t interQ)//enable_bits[2:0]={IsReLU, LoadBias, IsNotConv}
{
#pragma HLS INTERFACE m_axi depth=512 port=ifm    offset=slave bundle=DATA_BUS num_read_outstanding=1 num_write_outstanding=1 max_read_burst_length=256 max_write_burst_length=256
#pragma HLS INTERFACE m_axi depth=512 port=ofm    offset=slave bundle=DATA_BUS num_read_outstanding=1 num_write_outstanding=1 max_read_burst_length=256 max_write_burst_length=256
#pragma HLS INTERFACE m_axi depth=512 port=weight offset=slave bundle=DATA_BUS1 num_read_outstanding=1 max_read_burst_length=256
#pragma HLS INTERFACE m_axi depth=512 port=bias   offset=slave bundle=DATA_BUS1 num_read_outstanding=1 max_read_burst_length=256

#pragma HLS INTERFACE s_axilite register port=return bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=k_s_pad_ltype bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=iofm_num bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=ifm_w_h bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=ofm_w_h bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=TRTC bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=TMTN bundle=CTRL_BUS

#pragma HLS INTERFACE s_axilite register port=NToy bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=NTox bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=NTof bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=NTcomb bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=NTif bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=lmode bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=NTcomb_l bundle=CTRL_BUS

#pragma HLS INTERFACE s_axilite register port=pad_val bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=div_kk bundle=CTRL_BUS
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
	// assert((TR > 0)&&(TR <= Tr));
	// assert((TC > 0)&&(TC <= Tc));
	assert((TR*TC > 0)&&(TR*TC <= Tr*Tc));
	assert(Tm%LANE_NUM==0);
	assert(Tn%LANE_NUM==0);

	static int32_t ofm_buffer0[Tmax_dx][TrxTc][LANE_NUM];
#pragma HLS ARRAY_PARTITION variable=ofm_buffer0 complete dim=1
	static int32_t ofm_buffer1[Tmax_dx][TrxTc][LANE_NUM];
#pragma HLS ARRAY_PARTITION variable=ofm_buffer1 complete dim=1
	static int16_t bias_buffer[((MAX_BETA_LENGTH+LANE_NUM-1)/LANE_NUM)*LANE_NUM];

	uint16_t tx_n10[3], tx_n11[3];//tmtrtc, between compute2store
#pragma HLS ARRAY_PARTITION variable=tx_n10 complete dim=1
#pragma HLS ARRAY_PARTITION variable=tx_n11 complete dim=1
	uint16_t TX_MIN_n10[3], TX_MIN_n11[3];
#pragma HLS ARRAY_PARTITION variable=TX_MIN_n10 complete dim=1
#pragma HLS ARRAY_PARTITION variable=TX_MIN_n11 complete dim=1

	if(LoadBias)
		memcpy(bias_buffer, bias, (((ofm_num+LANE_NUM-1)/LANE_NUM)*LANE_NUM)*sizeof(int16_t));

	uint16_t tr, tc, tm;
	uint16_t TR_MIN, TC_MIN, TM_MIN;
	uint16_t tm_r, tr_r, tc_r;
	bool in_flag, proc_flag, out_flag, lc_enable;

	// uint8_t interQ = INTER_SF;
	uint8_t ComQ = 0;
	bool C_SL_EN = 1;
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
	if(!((ComQ >=0)&&(ComQ <32))){
		printf("ComQ=%d, interQ=%d, ifmQ=%d, weightQ=%d\n", ComQ, interQ, ifmQ, weightQ);
	}
	assert((ComQ >=0)&&(ComQ <32));
	
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
				tx_n10, TX_MIN_n10, pp, in_flag, proc_flag, pad_val, LoadBias, NTif, lmode, lc_enable, r_init_en, c_init_en, ComQ, C_SL_EN);

			write_back_output_reorg( ofm_buffer1, bias_buffer, ofm, tx_n11[1], tx_n11[2], tx_n11[0], ofm_num, ofm_h, ofm_w, TX_MIN_n11[0], TX_MIN_n11[1], TX_MIN_n11[2], OHW, IsReLU, div_kk,
				IsAVGPOOL, LoadBias, out_flag, biasQ, ofmQ, avgQ, interQ);

			pp = 0;
		}else{
			load_compute_wrapper(ifm, weight, ofm_buffer1, ksize, kstride, ifm_num, ifm_w, ifm_h, ofm_num, ofm_h, ofm_w,
				pad_int, ltype, IHW, KK, INumxKK, TM, TN, TR, TC, tm_r, tr_r, tc_r,
				tx_n11, TX_MIN_n11, pp, in_flag, proc_flag, pad_val, LoadBias, NTif, lmode, lc_enable, r_init_en, c_init_en, ComQ, C_SL_EN);			

			write_back_output_reorg( ofm_buffer0, bias_buffer, ofm, tx_n10[1], tx_n10[2], tx_n10[0], ofm_num, ofm_h, ofm_w, TX_MIN_n10[0], TX_MIN_n10[1], TX_MIN_n10[2], OHW, IsReLU, div_kk,
				IsAVGPOOL, LoadBias, out_flag, biasQ, ofmQ, avgQ, interQ);

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
