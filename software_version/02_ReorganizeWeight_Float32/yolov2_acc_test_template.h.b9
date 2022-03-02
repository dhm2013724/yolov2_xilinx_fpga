

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

#define PRAGMA_SUB(x) _Pragma (#x)
#define DO_PRAGMA(x) PRAGMA_SUB(x)

const uint32_t IB_W = OnChipIB_Width;
const uint32_t IB_H = OnChipIB_Height;
const uint32_t IB_HxW = IB_H*IB_W;
const uint32_t TnxIB_H = Tn*IB_H;	
const uint32_t TnxIB_HxIB_W = Tn*IB_H*IB_W;
const uint32_t TrxTc = Tr*Tc;

void input_load(float *input,float input_buffer[Tn][IB_HxW], uint16_t r, uint16_t c, uint16_t n, uint8_t Kstride, uint8_t ksize, 
		 uint8_t Padding, uint16_t TR_MIN, uint16_t TC_MIN, uint16_t Input_w, uint16_t Input_h, uint8_t TN_MIN, uint32_t IHxIW, float pad_val, 
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
				// if(IsColCont)
				if(t3==col_len){
					t3 = 0;
					t2++;
					// if(IsAllCont)
					if(t2==row_len){
							t2 = 0;
							t1++;
				}}
			}

	}}

	if(px_lr){
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
				if(t1<TN_MIN)
					input_buffer[t1][ifm_idx] = pad_val;
			}
		}}		
	}

	if(px_tb){
		for(uint16_t t2 = 0;t2 < px_tb; t2++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=IB_H)
		for(uint16_t t3 = 0;t3 < TCol_MIN; t3++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=IB_W)
#pragma HLS PIPELINE II=1
			uint16_t t2_idx;
			uint16_t t3_idx = t3;
			if(t2 < px_t)
				t2_idx = t2;
			else
				t2_idx = t2+row_len;
			uint32_t ifm_idx = t2_idx*TCol_MIN + t3_idx;

			for(uint16_t t1 = 0;t1 < Tn; t1++){
				if(t1<TN_MIN)
					input_buffer[t1][ifm_idx] = pad_val;
			}
		}}		
	}	

//because we have set un-used kernel weight to be 0, we dont need to set 0 here.	
// 	uint8_t TN_left = Tn - TN_MIN;
// 	if(IsPad0ExtraFM&&TN_left){			
// 		for(uint16_t t2 = 0;t2 < TRow_MIN; t2++){
// DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=IB_H)
// 		for(uint16_t t3 = 0;t3 < TCol_MIN; t3++){
// DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=IB_W)
// #pragma HLS PIPELINE II=1
// 			uint32_t ifm_idx = t2*TCol_MIN + t3;
// 			for(uint16_t t1 = 0;t1 < Tn; t1++){
// 				if(t1>=TN_MIN)
// 					input_buffer[t1][ifm_idx] = 0;
// 			}
// 		}}		
// 	}

}

const uint32_t TmxTnxKxK = Tm*Tn*K*K;

void weight_load_reorg(float *Weight,float weight_buffer[Tm][Tn][K][K],bool w_enable,int m,int n, uint32_t IFM_numxKxK, uint8_t KxK, uint8_t Ksize, uint8_t TM_MIN, uint8_t TN_MIN)
{
	if(!w_enable)
		return;

	static int Woffset;

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

void copy_input_weight(float *input,float *Weight,int IFM_num,int Input_w, int Input_h, bool r_init_en, bool c_init_en, uint16_t TCol_MIN_r, 
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

	// input_load(input, input_buffer, r, c, n, Kstride, Padding, TRow, TCol, Input_w, Input_h, TN_MIN, IHxIW, LayerType);
	input_load(input, input_buffer, r, c, n, Kstride, Ksize, Padding, TR_MIN, TC_MIN, Input_w, Input_h, TN_MIN, IHxIW, pad_val, r_init_en, c_init_en, TCol_MIN_r);
	weight_load_reorg(Weight,weight_buffer,weight_load_enable,m,n,IFM_numxKxK,KxK,Ksize,TM_MIN,TN_MIN);
}

void conv2d_tile(float input_buffer[Tn][IB_HxW],float output_buffer[Tm][TrxTc],
		float weight_buffer[Tm][Tn][K][K],int n_next,
		const int Ksize,const int Kstride, uint16_t TCol_MIN,
		const int TM_MIN,const int TR_MIN,const int TC_MIN,bool enable)
{
	if(!enable)
	{
		return;
	}

	// uint16_t TCol_MIN_l = (TC_MIN-1)*Kstride + Ksize;

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

void write_back_output_reorg(float output_buffer[Tm][TrxTc], float bias_buffer[MAX_BETA_LENGTH], /*float bias_buffer[MAX_BETA_LENGTH],*/
		float *Output, uint16_t r,uint16_t c,uint16_t m, 
		uint16_t ofm_num, uint16_t ofm_h, uint16_t ofm_w,
		uint8_t TM_MIN, uint16_t TR_MIN, uint16_t TC_MIN, uint32_t OHxOW, bool IsNL, bool IsBias, bool enable)
{
	if(!enable)
		return;

	assert((TM_MIN >0)&&(TM_MIN <=Tm));
	// assert((TR_MIN >0)&&(TR_MIN <=Tr));
	// assert((TC_MIN >0)&&(TC_MIN <=Tc));

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

// void write_back_output_reorg(float output_buffer[Tm][Tr][Tc], float bias_buffer[MAX_BETA_LENGTH], /*float bias_buffer[MAX_BETA_LENGTH],*/
// 		float *Output, uint16_t r,uint16_t c,uint16_t m, 
// 		uint16_t ofm_num, uint16_t ofm_h, uint16_t ofm_w,
// 		uint8_t TM_MIN, uint8_t TR_MIN, uint8_t TC_MIN, uint32_t OHxOW, bool IsNL, bool IsBias, bool enable)
// {
// 	if(!enable)
// 		return;

// 	assert((TM_MIN >0)&&(TM_MIN <=Tm));
// 	assert((TR_MIN >0)&&(TR_MIN <=Tr));
// 	assert((TC_MIN >0)&&(TC_MIN <=Tc));
// 	assert((OFM_BLmax>=1)&&(OFM_BLmax<=256));

// 	uint32_t offset = m*OHxOW + r*ofm_w + c;

// 	bool IsColCont = (ofm_w==TC_MIN);//equal FMCont
// 	bool IsRowCont = (ofm_h==TR_MIN);
// 	bool IsAllCont = IsColCont && IsRowCont;

// 	int burstlen;
// 	uint16_t T1_MIN, T2_MIN;
// 	if(IsAllCont){
// 		burstlen = TM_MIN*TR_MIN*TC_MIN;
// 		T1_MIN = 1;
// 		T2_MIN = 1;
// 	}else if(IsColCont){
// 		burstlen = TR_MIN*TC_MIN;
// 		T1_MIN = TM_MIN;
// 		T2_MIN = 1;
// 	}else{
// 		burstlen = TC_MIN;
// 		T1_MIN = TM_MIN;
// 		T2_MIN = TR_MIN;
// 	}

// 	// uint32_t burstlen_a_blmax = burstlen + OFM_BLmax;
// 	uint8_t add_1b;
// 	if(burstlen & 0xFF)
// 		add_1b = 1;
// 	else
// 		add_1b = 0;

// 	uint32_t blmax_tc =  (burstlen >> 8) + add_1b;
// 	uint32_t Tcomb_TC = T1_MIN*T2_MIN*blmax_tc;

// 	uint16_t tr,tm,tc;
// 	uint16_t t1, t2;
// 	uint32_t tbl;
	
// 	t1 = 0; t2 = 0; tbl = 0;
// 	for(uint32_t tcb = 0; tcb < Tcomb_TC; tcb++){	

// 		uint32_t tbl_r = tbl << 8;
// 		if(tbl==0){
// 			tc = 0;
// 			if(IsAllCont){
// 				tm = 0;
// 				tr = 0;					
// 			}else if(IsColCont){
// 				tm = t1;
// 				tr = 0;
// 			}else{
// 				tm = t1;
// 				tr = t2;
// 			}
// 		}

// 		uint32_t ofm_offset0 = offset + tbl_r;
// 		uint32_t ofm_offset;
// 		if(IsAllCont){
// 			ofm_offset = ofm_offset0;
// 		}else if(IsColCont){
// 			ofm_offset = ofm_offset0 + t1*OHxOW;
// 		}else{
// 			ofm_offset = ofm_offset0 + t1*OHxOW + t2*ofm_w;
// 		}

// 		float *OFM_base = Output + ofm_offset;
// 		uint16_t TBL_MIN = MIN(OFM_BLmax, burstlen-tbl_r);

// 		for(int tbl_min=0; tbl_min < TBL_MIN; tbl_min++){
// DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=OFM_BLmax)
// #pragma HLS PIPELINE II=1		
// 			float tmp_out = postproc(output_buffer[tm][tr][tc], bias_buffer[m+tm], IsBias, IsNL);

// 			tc++;
// 			if(tc==TC_MIN){
// 				tc = 0;
// 				tr++;
// 				if(tr==TR_MIN){
// 					tr = 0;
// 					tm++;
// 				}				
// 			}

// 			OFM_base[tbl_min] = tmp_out;				
// 		}

// 		tbl++;
// 		if(tbl==blmax_tc){
// 			tbl = 0;
// 			t2++;
// 			if(t2==T2_MIN){
// 				t2=0;
// 				t1++;
// 			}
// 		}
// 	}
// }

void pool_yolo2(float Input[Tn][IB_HxW],float Output[Tm][TrxTc],
		  const int Ksize,const int Kstride, uint16_t TCol_MIN, const int TM_MIN,const int TR_MIN,const int TC_MIN, bool enable)
{
	if(!enable)
		return;

	// uint16_t TCol_MIN_l = (TC_MIN-1)*Kstride + Ksize;

	float tmp[Tn];
	for(uint16_t tr = 0;tr < TR_MIN;tr++)
	for(uint16_t tc = 0;tc < TC_MIN;tc++)
	for(uint16_t i =0;i < Ksize; i++)
	for(uint16_t j = 0;j < Ksize; j++)
	{
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
	}
}

void load_compute_wrapper(float *ifm, float *weight, float ofm_buffer[Tm][TrxTc], int ksize, int kstride, int ifm_num, int ifm_w, int ifm_h,int ofm_num,
	 int ofm_h, int ofm_w, int pad_int, int ltype, int IHW, int KK, int INumxKK, int TM, int TN, int TR, int TC, int tm_r, int tr_r, int tc_r,
	 int tx_n1[3],int TX_MIN_n1[3],bool pp,bool in_flag,bool proc_flag, bool LoadBias, int NTif, uint8_t lmode, bool enable,
	 bool r_init_en, bool c_init_en)
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

	static int tx_n[2][3];//tmtrtc, between load2compute
#pragma HLS ARRAY_PARTITION variable=tx_n complete dim=1
#pragma HLS ARRAY_PARTITION variable=tx_n complete dim=2	
	static int TX_MIN_n[2][3];
#pragma HLS ARRAY_PARTITION variable=TX_MIN_n complete dim=1
#pragma HLS ARRAY_PARTITION variable=TX_MIN_n complete dim=2	

	if(!enable)
		return ;

	// int	TR_MIN = MIN(TR,ofm_h-tr_r);
	// int	TC_MIN = MIN(TC,ofm_w-tc_r);
	// int	TM_MIN = MIN(TM,ofm_num-tm_r);
	// uint16_t TCol_MIN = (TC_MIN-1)*kstride + ksize;
	// static uint16_t TCol_MIN0[1], TCol_MIN1[1];

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
					TM_MIN,TN,TR_MIN,TC_MIN, pad_int,ifm_buffer0,weight_buffer0,1,IHW,KK,INumxKK,ltype,in_flag);
				conv2d_tile(ifm_buffer1, ofm_buffer, weight_buffer1, 0, ksize, kstride, TCol_MIN1[0], TX_MIN_n[1][0], TX_MIN_n[1][1], TX_MIN_n[1][2], proc_flag);
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
					TM_MIN,TN,TR_MIN,TC_MIN, pad_int,ifm_buffer1,weight_buffer1,1,IHW,KK,INumxKK,ltype,in_flag);					
				conv2d_tile(ifm_buffer0, ofm_buffer, weight_buffer0, 0, ksize, kstride, TCol_MIN0[0], TX_MIN_n[0][0], TX_MIN_n[0][1], TX_MIN_n[0][2], proc_flag);
			}else if(ltype == 1)
			{
				copy_input_weight(ifm,weight,ifm_num,ifm_w,ifm_h, r_init_en, c_init_en, TCol_MIN, ksize,kstride,tr_r,tc_r,tm_r,tm_r,
					TM_MIN,TM,TR_MIN,TC_MIN, 0,ifm_buffer1,weight_buffer1,0,IHW,KK,INumxKK,ltype,in_flag);						
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
				int tn_n[2];
				uint16_t tn_r = 0;
				for(int tn = 0; tn < NTif+1; tn++)
				{
					bool in_flag = tn < NTif;
					bool proc_flag = tn > 0;						
					// uint16_t tn_r = tn*TN;
					if(pp_tn){
						copy_input_weight(ifm,weight,ifm_num,ifm_w,ifm_h, r_init_en, c_init_en, TCol_MIN, ksize,kstride,tr_r,tc_r,tm_r, tn_r,
							TM_MIN,TN,TR_MIN,TC_MIN, pad_int,ifm_buffer0,weight_buffer0,1,IHW,KK,INumxKK,ltype,in_flag);
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
							TM_MIN,TN,TR_MIN,TC_MIN, pad_int,ifm_buffer1,weight_buffer1,1,IHW,KK,INumxKK,ltype,in_flag);						
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
							   int32_t NToy, int32_t NTox, int32_t NTof, int32_t NTcomb, int32_t NTif, uint8_t lmode, int32_t NTcomb_l, int LayerType)
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

			write_back_output_reorg( ofm_buffer1, bias_buffer, Output, tx_n11[1], tx_n11[2], tx_n11[0], OFM_num, Output_h, Output_w, TX_MIN_n11[0], TX_MIN_n11[1], TX_MIN_n11[2], OHxOW, IsNL, IsBias, out_flag);
			pp = 0;
		}else{
			load_compute_wrapper(Input, Weight, ofm_buffer1, Ksize, Kstride, IFM_num, Input_w, Input_h, OFM_num, Output_h, Output_w,
				Padding, LayerType, IHxIW, KxK, IFM_numxKxK, TM, TN, TR, TC, tm_r, tr_r, tc_r,
				tx_n11, TX_MIN_n11, pp, in_flag, proc_flag, IsBias, NTif, lmode, lc_enable, r_init_en, c_init_en);	

			write_back_output_reorg( ofm_buffer0, bias_buffer, Output, tx_n10[1], tx_n10[2], tx_n10[0], OFM_num, Output_h, Output_w, TX_MIN_n10[0], TX_MIN_n10[1], TX_MIN_n10[2], OHxOW, IsNL, IsBias, out_flag);
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

#define MEM_LEN (416*416*32+208*208*32)
void generate_iofm_offset(float* in_ptr[32], float* out_ptr[32], float *Memory_buf, network *net)
{
#define ROUTE16_LEN (26*26*512)
#define CONV27_LEN (13*13*256)
#define CONV24_LEN (13*13*1024)

	float *Memory_top = Memory_buf;
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
	float *Memory_buf = (float*)calloc(MEM_LEN,sizeof(float));
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

				// TR = MIN(((OnChipIB_Height-l.size)/l.stride+1),Tr);//keep Kstride>=1
				// TR = MIN(output_h,TR);
				// TC = MIN(((OnChipIB_Width-l.size)/l.stride+1),Tc);
				// TC = MIN(output_w,TC);

				TC = MIN(((IB_HxW-l.size)/l.stride+1),output_w);
				TCol = (TC-1)*l.stride + l.size;
				TR = MIN(((IB_HxW/TCol-l.size)/l.stride+1),output_h);//keep Kernel_stride>=1
				TR = MIN(TR, TrxTc/TC);
				TRow = (TR-1)*l.stride + l.size;

				// assert(((TR*TC)>0)&&((TR*TC)<=TrxTc));
				// assert(((TRow*TCol)>0)&&((TRow*TCol)<=IB_HxW));
				// printf("TR=%d, TC=%d, TRow=%d, TCol=%d\n", TR, TC, TRow, TCol);

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

				// TR = MIN(((OnChipIB_Height-l.size)/l.stride+1),Tr);//keep Kstride>=1
				// TC = MIN(((OnChipIB_Width-l.size)/l.stride+1),Tc);
				// TR = MIN(output_h,TR);
				// TC = MIN(output_w,TC);

				TC = MIN(((IB_HxW-l.size)/l.stride+1),output_w);
				TCol = (TC-1)*l.stride + l.size;
				TR = MIN(((IB_HxW/TCol-l.size)/l.stride+1),output_h);//keep Kernel_stride>=1
				TR = MIN(TR, TrxTc/TC);
				TRow = (TR-1)*l.stride + l.size;
				
				// assert(((TR*TC)>0)&&((TR*TC)<=TrxTc));
				// assert(((TRow*TCol)>0)&&((TRow*TCol)<=IB_HxW));
				// printf("TR=%d, TC=%d, TRow=%d, TCol=%d\n", TR, TC, TRow, TCol);

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
