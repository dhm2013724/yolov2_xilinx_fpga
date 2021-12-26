
#include "cnn.h"

const uint32_t IB_W = OnChipIB_Width;
const uint32_t IB_H = OnChipIB_Height;
const uint32_t TnxIB_H = Tn*IB_H;
const uint32_t TnxIB_HxIB_W = Tn*IB_H*IB_W;

void input_load(float *input,float input_buffer[Tn][OnChipIB_Height][OnChipIB_Width], uint16_t r, uint16_t c, uint16_t n, uint8_t Kstride, uint8_t ksize, uint8_t Padding,
		 uint16_t TR_MIN, uint16_t TC_MIN,
		 uint16_t Input_w, uint16_t Input_h, uint8_t TN_MIN, uint32_t IHxIW, float pad_val, bool IsPad0ExtraFM)
{
	const int16_t Coffset = c*Kstride - Padding;
	const int16_t Roffset = r*Kstride - Padding;

	uint16_t TRow_MIN = (TR_MIN-1)*Kstride + ksize;
	uint16_t TCol_MIN = (TC_MIN-1)*Kstride + ksize;

	assert(TN_MIN < 65535);
	assert(TRow_MIN < 65535);
	assert(TCol_MIN < 65535);

	uint16_t tl_y, tl_x, br_y, br_x;
	uint16_t px_l, px_r, px_b, px_t;
	if(Coffset < 0){
		tl_x = 0;
		px_l = -Coffset;
	}
	else{
		tl_x = Coffset;
		px_l = 0;
	}

	br_x = Coffset + TCol_MIN -1;
	if(br_x < Input_w){
		px_r = 0;
	}else{
		px_r = br_x - Input_w + 1;
	}

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

	uint16_t px_tb = (px_t + px_b);
	uint16_t px_lr = (px_l + px_r);
	uint16_t col_len = TCol_MIN - px_lr;
	uint16_t row_len = TRow_MIN - px_tb;
	uint16_t col_max = col_len + px_l;
	uint16_t row_max = row_len + px_t;

	bool IsColCont = (col_len == Input_w);
	bool IsRowCont = (row_len == Input_h);
	bool IsAllCont = IsRowCont && IsColCont;

	uint32_t burstlen;
	uint32_t offset_base;
	uint16_t T1_MIN, T2_MIN, T12_MIN;
	uint8_t ifm_lmode;
	if(IsAllCont){
		burstlen = TN_MIN*row_len*col_len;
		offset_base = n*IHxIW;
		T1_MIN = 1;
		T2_MIN = 1;
		ifm_lmode = 0;
		T12_MIN = 1;
	}else if(IsColCont){
		burstlen = row_len*col_len;
		offset_base = n*IHxIW + tl_y*Input_w;
		T1_MIN = TN_MIN;
		T2_MIN = 1;
		ifm_lmode = 1;
		T12_MIN = TN_MIN;
	}else{
		burstlen = col_len;
		offset_base = n*IHxIW + tl_y*Input_w + tl_x;
		T1_MIN = TN_MIN;
		T2_MIN = row_len;
		ifm_lmode = 2;
		T12_MIN = TN_MIN*row_len;
	}
	assert((ifm_lmode >=0)&&(ifm_lmode <=2));

	uint16_t t1, t2, t3;
	for(uint16_t t12=0; t12<T12_MIN; t12++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=TnxIB_H)

		uint16_t t1_in, t2_in;
		uint32_t offset_ex = offset_base;

		static uint16_t t1m, t2m;
		if(t12==0){
			t1m = 0;
			t2m = 0;
		}else{
			t2m++;
			if(t2m==T2_MIN){
				t2m = 0;
				t1m++;
			}
		}

		if(ifm_lmode==0){
			t1 = 0; t2 = 0; t3 = 0;
		}else if(ifm_lmode==1)
		{
			t1_in = t1m; t2 = 0; t3 = 0;
			offset_ex += (t1m*IHxIW);
		}else if(ifm_lmode==2)
		{
			t1_in = t1m; t2_in = t2m+px_t; t3 = 0;
			offset_ex += ( t1m*IHxIW + t2m*Input_w);
		}

		float ibuf_local[TnxIB_HxIB_W];
		float *ifm_addr = input + offset_ex;
		assert(burstlen <= TnxIB_HxIB_W);

		uint32_t ifm_offset = 0;
		if(ifm_lmode==0)
		{
			for(uint16_t t=0; t<burstlen; t++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=TnxIB_HxIB_W)
#pragma HLS PIPELINE II=1
				float tmp_in = ifm_addr[t];
				input_buffer[t1][t2+px_t][t3+px_l] = tmp_in;

				t3++;
				if(t3==col_len){
					t3=0; t2++;
					if(t2==row_len){
						t2=0; t1++;
					}
				}
			}
		}else if(ifm_lmode==1)
		{
			for(uint16_t t=0; t<burstlen; t++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=TnxIB_HxIB_W)
#pragma HLS PIPELINE II=1
				float tmp_in = ifm_addr[t];
				input_buffer[t1_in][t2+px_t][t3+px_l] = tmp_in;

				t3++;
				if(t3==col_len){
					t3=0; t2++;
				}
			}
		}else
		{
			for(uint16_t t=0; t<burstlen; t++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=TnxIB_HxIB_W)
#pragma HLS PIPELINE II=1
				float tmp_in = ifm_addr[t];
				input_buffer[t1_in][t2_in][t3+px_l] = tmp_in;

				t3++;
			}
		}
	}

	if(px_lr){
		for(uint16_t t1 = 0;t1 < TN_MIN; t1++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tn)
		for(uint16_t t2 = 0;t2 < row_len; t2++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=IB_H)
		for(uint16_t t3 = 0;t3 < px_lr; t3++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=IB_W)
#pragma HLS PIPELINE II=1
			if(t3 < px_l)
				input_buffer[t1][t2+px_t][t3] = pad_val;
			else
				input_buffer[t1][t2+px_t][t3+col_max] = pad_val;
		}}}
	}

	if(px_tb){
		for(uint16_t t1 = 0;t1 < TN_MIN; t1++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tn)
		for(uint16_t t2 = 0;t2 < px_tb; t2++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=IB_H)
		for(uint16_t t3 = 0;t3 < TCol_MIN; t3++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=IB_W)
#pragma HLS PIPELINE II=1
			if(t2 < px_t)
				input_buffer[t1][t2][t3] = pad_val;
			else
				input_buffer[t1][t2+row_len][t3] = pad_val;
		}}}
	}

	uint16_t TN_left = Tn - TN_MIN;
	if(IsPad0ExtraFM&&TN_left){
		for(uint16_t t1 = 0;t1 < TN_left; t1++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tn)
		for(uint16_t t2 = 0;t2 < TRow_MIN; t2++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=IB_H)
		for(uint16_t t3 = 0;t3 < TCol_MIN; t3++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=IB_W)
#pragma HLS PIPELINE II=1
				input_buffer[t1+TN_MIN][t2][t3] = 0;
		}}}
	}

}

//void memcpy_ifm(hls::stream<float> &ibuf_local,/*float *ibuf_local,*/  float *input, uint32_t burstlen, uint32_t offset_base, uint8_t ifm_lmode, uint16_t Input_w, uint32_t IHxIW, uint16_t px_t,
//				uint16_t t12, uint16_t T2_MIN, uint16_t *t1, uint16_t *t2){
//
//#pragma HLS INLINE off
//	uint32_t offset_ex = offset_base;
//	assert((ifm_lmode >=0)&&(ifm_lmode <=2));
//
//	static uint16_t t1m, t2m;
//
//	if(t12==0){
//		t1m = 0;
//		t2m = 0;
//	}else{
//		t2m++;
//		if(t2m==T2_MIN){
//			t2m = 0;
//			t1m++;
//		}
//	}
//
//	if(ifm_lmode==0){
//		*t1 = 0; *t2 = 0;
//	}else if(ifm_lmode==1)
//	{
//		*t1 = t1m; *t2 = 0;
//		offset_ex += (t1m*IHxIW);
//	}else if(ifm_lmode==2)
//	{
//		*t1 = t1m; *t2 = t2m + px_t;
//		offset_ex += ( t1m*IHxIW + t2m*Input_w);
//	}
//
//	float *ifm_addr = input + offset_ex;
//	assert(burstlen <= TnxIB_HxIB_W);
//
//	for(uint16_t t=0; t<burstlen; t++){
//DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=TnxIB_HxIB_W)
//#pragma HLS PIPELINE II=1
////		ibuf_local[t] = ifm_addr[t];
//		float tmp_f32 = ifm_addr[t];
//		ibuf_local.write(tmp_f32);
//	}
//}
//
//void memcpy_ifm2buf(float input_buffer[Tn][OnChipIB_Height][OnChipIB_Width], hls::stream<float> &ibuf_local,/*float *ibuf_local,*/
//		uint8_t TN_MIN, uint16_t TRow_MIN, uint16_t TCol_MIN, float pad_val, uint8_t ifm_lmode,
//		uint16_t px_l, uint16_t px_t, uint16_t row_max, uint16_t col_max, uint16_t t1_in, uint16_t t2_in, uint16_t T1_MIN2B, uint16_t T2_MIN2B){
//
//#pragma HLS INLINE off
////	uint32_t ifm_offset = 0;
//	if((ifm_lmode==0) || (ifm_lmode==1))
//	{
//		for(uint16_t t1 = 0;t1 < T1_MIN2B; t1++){
//DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tn)
//		for(uint16_t t2 = 0;t2 < T2_MIN2B; t2++){
//DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=IB_H)
//		for(uint16_t t3 = 0;t3 < TCol_MIN; t3++){
//DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=IB_W)
//#pragma HLS PIPELINE II=1
//			bool PX_R_en = (t2 >= px_t) && (t2 < row_max);
//			bool PX_C_en = (t3 >= px_l) && (t3 < col_max);
//			if(PX_R_en&&PX_C_en){
//				input_buffer[t1+t1_in][t2][t3] = ibuf_local.read();
////				input_buffer[t1+t1_in][t2][t3] = ibuf_local[ifm_offset];
////				ifm_offset++;
//			}else{
//				input_buffer[t1+t1_in][t2][t3] = pad_val;
//			}
//		}}}
//	}else
//	{
//		for(uint16_t t3 = 0;t3 < TCol_MIN; t3++){
//DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=IB_W)
//#pragma HLS PIPELINE II=1
//			bool PX_C_en = (t3 >= px_l) && (t3 < col_max);
//			if(PX_C_en){
//				input_buffer[t1_in][t2_in][t3] = ibuf_local.read();
////				input_buffer[t1_in][t2_in][t3] = ibuf_local[ifm_offset];
////				ifm_offset++;
//			}else{
//				input_buffer[t1_in][t2_in][t3] = pad_val;
//			}
//		}
//	}
//}
//
//void load_ifm_df_wrapper(float *input,float input_buffer[Tn][OnChipIB_Height][OnChipIB_Width], uint16_t Input_w, uint8_t TN_MIN, uint32_t IHxIW, float pad_val,
//	uint32_t burstlen, uint32_t offset_base, uint8_t ifm_lmode, uint16_t TRow_MIN, uint16_t TCol_MIN,
//	uint16_t px_l, uint16_t px_t, uint16_t row_max, uint16_t col_max, uint16_t T1_MIN2B, uint16_t T2_MIN2B, uint16_t T1_MIN, uint16_t T2_MIN, uint16_t T12_MIN){
//
//	for(uint16_t t12=0; t12<T12_MIN; t12++){
//DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=TnxIB_H)
//#pragma HLS DATAFLOW
//		uint16_t t1_in, t2_in;
////		float ibuf_local[TnxIB_HxIB_W];
//		hls::stream<float> ibuf_local;
//#pragma HLS STREAM variable=ibuf_local depth=TnxIB_HxIB_W dim=1
//
//		memcpy_ifm(ibuf_local, input, burstlen, offset_base, ifm_lmode, Input_w, IHxIW, px_t, t12, T2_MIN, &t1_in, &t2_in);
//
//		memcpy_ifm2buf(input_buffer, ibuf_local, TN_MIN, TRow_MIN, TCol_MIN, pad_val, ifm_lmode,
//						px_l, px_t, row_max, col_max, t1_in, t2_in, T1_MIN2B, T2_MIN2B);
//	}
//
//}
//
//void input_load(float *input,float input_buffer[Tn][OnChipIB_Height][OnChipIB_Width], uint16_t r, uint16_t c, uint16_t n, uint8_t Kstride, uint8_t ksize, uint8_t Padding,
//		 uint16_t TR_MIN, uint16_t TC_MIN,
//		 uint16_t Input_w, uint16_t Input_h, uint8_t TN_MIN, uint32_t IHxIW, float pad_val, bool IsPad0ExtraFM)
//{
//	const int16_t Coffset = c*Kstride - Padding;
//	const int16_t Roffset = r*Kstride - Padding;
//
//	uint16_t TRow_MIN = (TR_MIN-1)*Kstride + ksize;
//	uint16_t TCol_MIN = (TC_MIN-1)*Kstride + ksize;
//
//	assert(TN_MIN < 65535);
//	assert(TRow_MIN < 65535);
//	assert(TCol_MIN < 65535);
//
//	uint16_t tl_y, tl_x, br_y, br_x;
//	uint16_t px_l, px_r, px_b, px_t;
//	if(Coffset < 0){
//		tl_x = 0;
//		px_l = -Coffset;
//	}
//	else{
//		tl_x = Coffset;
//		px_l = 0;
//	}
//
//	br_x = Coffset + TCol_MIN -1;
//	if(br_x < Input_w){
//		px_r = 0;
//	}else{
//		px_r = br_x - Input_w + 1;
//	}
//
//	if(Roffset < 0){
//		tl_y = 0;
//		px_t = -Roffset;
//	}
//	else{
//		tl_y = Roffset;
//		px_t = 0;
//	}
//	br_y = Roffset + TRow_MIN -1;
//	if(br_y < Input_h){
//		px_b = 0;
//	}else{
//		px_b = br_y - Input_h + 1;
//	}
//
//	uint16_t px_tb = (px_t + px_b);
//	uint16_t px_lr = (px_l + px_r);
//	uint16_t col_len = TCol_MIN - px_lr;
//	uint16_t row_len = TRow_MIN - px_tb;
//	uint16_t col_max = col_len + px_l;
//	uint16_t row_max = row_len + px_t;
//
//	bool IsColCont = (col_len == Input_w);
//	bool IsRowCont = (row_len == Input_h);
//	bool IsAllCont = IsRowCont && IsColCont;
//
//	uint32_t burstlen;
//	uint32_t offset_base;
//	uint16_t T1_MIN, T2_MIN, T1_MIN2B, T2_MIN2B, T12_MIN;
//	uint8_t ifm_lmode;
//	if(IsAllCont){
//		burstlen = TN_MIN*row_len*col_len;
//		offset_base = n*IHxIW;
//		T1_MIN = 1;
//		T2_MIN = 1;
//		ifm_lmode = 0;
//		T1_MIN2B = TN_MIN;
//		T2_MIN2B = TRow_MIN;
//		T12_MIN = 1;
//	}else if(IsColCont){
//		burstlen = row_len*col_len;
//		offset_base = n*IHxIW + tl_y*Input_w;
//		T1_MIN = TN_MIN;
//		T2_MIN = 1;
//		ifm_lmode = 1;
//		T1_MIN2B = 1;
//		T2_MIN2B = TRow_MIN;
//		T12_MIN = TN_MIN;
//	}else{
//		burstlen = col_len;
//		offset_base = n*IHxIW + tl_y*Input_w + tl_x;
//		T1_MIN = TN_MIN;
//		T2_MIN = row_len;
//		ifm_lmode = 2;
//		T1_MIN2B = 1;
//		T2_MIN2B = 1;
//		T12_MIN = TN_MIN*row_len;
//	}
//
//	load_ifm_df_wrapper(input, input_buffer, Input_w, TN_MIN, IHxIW, pad_val, burstlen, offset_base,
//		 ifm_lmode, TRow_MIN, TCol_MIN, px_l, px_t, row_max, col_max, T1_MIN2B, T2_MIN2B, T1_MIN, T2_MIN, T12_MIN);
//
//	if((ifm_lmode==2)&&px_tb){
//		for(uint16_t t1 = 0;t1 < TN_MIN; t1++){
//DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tn)
//		for(uint16_t t2 = 0;t2 < px_tb; t2++){
//DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=IB_H)
//		for(uint16_t t3 = 0;t3 < TCol_MIN; t3++){
//DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=IB_W)
//#pragma HLS PIPELINE II=1
//			if(t2 < px_t)
//				input_buffer[t1][t2][t3] = pad_val;
//			else
//				input_buffer[t1][t2+row_len][t3] = pad_val;
//		}}}
//	}
//
//	uint16_t TN_left = Tn - TN_MIN;
//	if(IsPad0ExtraFM&&TN_left){
//		for(uint16_t t1 = 0;t1 < TN_left; t1++){
//DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tn)
//		for(uint16_t t2 = 0;t2 < TRow_MIN; t2++){
//DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=IB_H)
//		for(uint16_t t3 = 0;t3 < TCol_MIN; t3++){
//DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=IB_W)
//#pragma HLS PIPELINE II=1
//				input_buffer[t1+TN_MIN][t2][t3] = 0;
//		}}}
//	}
//
//}

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

	uint8_t TM_left = Tm - TM_MIN;
	uint8_t TN_left = Tn - TN_MIN;

	if(TM_left || TN_left){

	for(t3 = 0;t3 <Ksize; t3++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=K)
	for(t4 = 0;t4 <Ksize; t4++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=K)
#pragma HLS PIPELINE II=1
	for(t1 = 0;t1 < Tm; t1++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tm)
	for(t2 = 0;t2 < Tn; t2++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tn)
			if((t1 >= TM_MIN)&&(t2 >= TN_MIN))
				weight_buffer[t1][t2][t3][t4] = 0;
	}}}}

	}
}


//void memcpy_w(hls::stream<float> &wbuf_local /*float *wbuf_local*/,  float *weight_addr, uint32_t burstlen){
//
//#pragma HLS INLINE off
//	for(uint16_t t=0; t<burstlen; t++){
//DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=TmxTnxKxK)
//#pragma HLS PIPELINE II=1
////		wbuf_local[t] = weight_addr[t];
//		float tmp_f32 = weight_addr[t];
//		wbuf_local.write(tmp_f32);
//	}
//}
//
//void memcpy_w2buf(float weight_buffer[Tm][Tn][K][K], hls::stream<float> &wbuf_local /*float *wbuf_local*/, uint8_t Ksize, uint8_t TM_MIN, uint8_t TN_MIN){
//
//#pragma HLS INLINE off
//	uint8_t t1,t2,t3,t4;
////	uint32_t woffset_local = 0;
//	for(t3 = 0;t3 <Ksize; t3++){
//DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=K)
//	for(t4 = 0;t4 <Ksize; t4++){
//DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=K)
//	for(t1 = 0;t1 < Tm; t1++){
//DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tm)
//	for(t2 = 0;t2 < Tn; t2++){
//DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tn)
//#pragma HLS PIPELINE II=1
//			if((t1 < TM_MIN)&&(t2 < TN_MIN)){
//				weight_buffer[t1][t2][t3][t4] = wbuf_local.read();
////				weight_buffer[t1][t2][t3][t4] = wbuf_local[woffset_local];
////				woffset_local++;
//			}else
//				weight_buffer[t1][t2][t3][t4] = 0;
//	}}}}
//}
//
//void memcpy_w_df_wrapper(float *weight_addr, uint32_t burstlen,
//					float weight_buffer[Tm][Tn][K][K], uint8_t Ksize, uint8_t TM_MIN, uint8_t TN_MIN){
//
//#pragma HLS DATAFLOW
////	float wbuf_local[TmxTnxKxK];
//	hls::stream<float> wbuf_local;
//#pragma HLS STREAM variable=wbuf_local depth=TmxTnxKxK dim=1
//
//	memcpy_w(wbuf_local,  weight_addr, burstlen);
//	memcpy_w2buf( weight_buffer, wbuf_local, Ksize, TM_MIN, TN_MIN);
//}
//
//void weight_load_reorg(float *Weight,float weight_buffer[Tm][Tn][K][K],bool w_enable,int m,int n, uint32_t IFM_numxKxK, uint8_t KxK, uint8_t Ksize, uint8_t TM_MIN, uint8_t TN_MIN)
//{
//	if(!w_enable)
//		return;
//
//	static int Woffset;
//
//	assert((TM_MIN > 0)&&(TM_MIN <= Tm));
//	assert((TN_MIN > 0)&&(TN_MIN <= Tn));
//	assert((KxK > 0)&&(KxK <= K*K));
//
//	if(m==0&&n==0)
//		Woffset = 0;
//
//	float *weight_addr = Weight + Woffset;
//	uint32_t burstlen = TM_MIN*TN_MIN*KxK;
//	Woffset += burstlen;
//	assert(burstlen <= TmxTnxKxK);
//
//	memcpy_w_df_wrapper(weight_addr, burstlen, weight_buffer, Ksize, TM_MIN, TN_MIN);
//}

void copy_input_weight(float *input,float *Weight,int IFM_num,int Input_w, int Input_h, int Ksize,int Kstride,int r,int c,int m,int n,
		int TM_MIN,int TN,int TR_MIN,int TC_MIN,int Padding,float input_buffer[Tn][OnChipIB_Height][OnChipIB_Width],float weight_buffer[Tm][Tn][K][K],
		bool weight_load_enable,const int IHxIW,const int KxK,const int IFM_numxKxK,const int LayerType, bool enable)
{
	if(!enable)
		return;

	int TN_MIN = MIN(TN, IFM_num-n);
	bool IsPad0ExtraFM = (LayerType==0);
	float pad_val = 0.0f;
	if(!IsPad0ExtraFM)
		pad_val = -1024*1024;

	input_load(input, input_buffer, r, c, n, Kstride, Ksize, Padding, TR_MIN, TC_MIN, Input_w, Input_h, TN_MIN, IHxIW, pad_val, IsPad0ExtraFM);
	weight_load_reorg(Weight,weight_buffer,weight_load_enable,m,n,IFM_numxKxK,KxK,Ksize,TM_MIN,TN_MIN);
}

void conv2d_tile(float input_buffer[Tn][OnChipIB_Height][OnChipIB_Width],float output_buffer[Tm][Tr][Tc],
		float weight_buffer[Tm][Tn][K][K], uint16_t n_next,
		uint8_t Ksize, uint8_t Kstride,
		uint16_t TM_MIN, uint16_t TR_MIN, uint16_t TC_MIN,bool enable)
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
#pragma HLS PIPELINE II=2
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

const int OFM_BLmax = 256;

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

void ofm_mmcpy_cont(float *Output, hls::stream<float> &local_buf /*float *local_buf*/, uint32_t ofm_offset, uint16_t data_num)
{
#pragma HLS INLINE off
	float *OFM_base = Output + ofm_offset;
	for(uint32_t x=0; x< data_num; x++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=OFM_BLmax)
#pragma HLS PIPELINE II=1
//		OFM_base[x] = local_buf[x];
		float tmp_f32 = local_buf.read();
		OFM_base[x] = tmp_f32;
	}
}

void pp_cont(hls::stream<float> &local_buf /*float *local_buf*/, float output_buffer[Tm][Tr][Tc], float *bias_buffer,
				uint8_t TR_MIN, uint8_t TC_MIN, uint32_t OHxOW, uint16_t ofm_w, bool IsNL, bool IsBias,
				uint16_t m, uint32_t offset, uint32_t burstlen, uint32_t blmax_tc, uint16_t T2_MIN,
				uint16_t *TBL_MIN_p, uint32_t *ofm_offset_p, bool IsAllCont, bool IsColCont, bool Init_en){

#pragma HLS INLINE off
	static uint16_t tr,tm,tc;
	static uint16_t t1, t2;
	static uint32_t tbl;

	if(Init_en){
		t1 = 0;
		t2 = 0;
		tbl = 0;
	}

	uint32_t tbl_r = tbl << 8;
	if(tbl_r==0){
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
	uint16_t TBL_MIN = MIN(OFM_BLmax, burstlen-tbl_r);
//	uint8_t bl_cnt = 0;
	for(int tbl_min=0; tbl_min < TBL_MIN; tbl_min++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=OFM_BLmax)
#pragma HLS PIPELINE II=1
		float tmp_out = postproc(output_buffer[tm][tr][tc], bias_buffer[m+tm], IsBias, IsNL);
//		local_buf[bl_cnt] = tmp_out;
//		bl_cnt++;
		local_buf.write(tmp_out);
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

	uint32_t ofm_offset0 = offset + tbl_r;
	uint32_t ofm_offset;
	if(IsAllCont){
		ofm_offset = ofm_offset0;
	}else if(IsColCont){
		ofm_offset = ofm_offset0 + t1*OHxOW;
	}else{
		ofm_offset = ofm_offset0 + t1*OHxOW + t2*ofm_w;
	}

	*TBL_MIN_p = TBL_MIN;
	*ofm_offset_p = ofm_offset;

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

void ofm_df_warpper(float output_buffer[Tm][Tr][Tc], float bias_buffer[MAX_BETA_LENGTH],
		float *Output,uint16_t m, uint16_t ofm_w, bool IsAllCont, bool IsColCont,
		uint32_t Tcomb_TC, uint32_t offset, uint32_t burstlen, uint32_t blmax_tc, uint16_t T2_MIN,
		uint8_t TR_MIN, uint8_t TC_MIN, uint32_t OHxOW, bool IsNL, bool IsBias){

	const uint32_t Tcomb_TC_MAX = Tm*Tr*(Tc/256+1);
	assert(Tcomb_TC <= Tcomb_TC_MAX);
	for(uint32_t tcb = 0; tcb < Tcomb_TC; tcb++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tcomb_TC_MAX)
#pragma HLS DATAFLOW
//	float local_buf[OFM_BLmax];
	hls::stream<float> local_buf;
#pragma HLS STREAM variable=local_buf depth=OFM_BLmax dim=1
	uint16_t TBL_MIN_p[1];
	uint32_t ofm_offset_p[1];

		pp_cont( local_buf, output_buffer, bias_buffer,
			TR_MIN, TC_MIN, OHxOW, ofm_w, IsNL, IsBias,
			m, offset, burstlen, blmax_tc, T2_MIN,
			TBL_MIN_p, ofm_offset_p, IsAllCont, IsColCont, tcb==0);

		ofm_mmcpy_cont(Output, local_buf, ofm_offset_p[0], TBL_MIN_p[0]);
	}

}

void write_back_output_reorg(float output_buffer[Tm][Tr][Tc], float bias_buffer[MAX_BETA_LENGTH], /*float bias_buffer[MAX_BETA_LENGTH],*/
		float *Output, uint16_t r,uint16_t c,uint16_t m,
		uint16_t ofm_num, uint16_t ofm_h, uint16_t ofm_w,
		uint8_t TM_MIN, uint8_t TR_MIN, uint8_t TC_MIN, uint32_t OHxOW, bool IsNL, bool IsBias, bool enable)
{
	if(!enable)
		return;

	assert((TM_MIN >0)&&(TM_MIN <=Tm));
	assert((TR_MIN >0)&&(TR_MIN <=Tr));
	assert((TC_MIN >0)&&(TC_MIN <=Tc));
	assert((OFM_BLmax>=1)&&(OFM_BLmax<=256));

	uint32_t offset = m*OHxOW + r*ofm_w + c;

	bool IsColCont = (ofm_w==TC_MIN);//equal FMCont
	bool IsRowCont = (ofm_h==TR_MIN);
	bool IsAllCont = IsColCont && IsRowCont;

	uint32_t burstlen;
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

	// uint32_t burstlen_a_blmax = burstlen + OFM_BLmax;
	uint8_t add_1b;
	if(burstlen & 0xFF)
		add_1b = 1;
	else
		add_1b = 0;

	uint32_t blmax_tc =  (burstlen >> 8) + add_1b;
	uint32_t Tcomb_TC = T1_MIN*T2_MIN*blmax_tc;

	ofm_df_warpper( output_buffer, bias_buffer, Output, m, ofm_w, IsAllCont, IsColCont,
		Tcomb_TC, offset, burstlen, blmax_tc, T2_MIN, TR_MIN, TC_MIN, OHxOW, IsNL, IsBias);

//	const uint32_t Tcomb_TC_MAX = Tm*Tr*(Tc/256+1);
//	assert(Tcomb_TC <= Tcomb_TC_MAX);
//	for(uint32_t tcb = 0; tcb < Tcomb_TC; tcb++){
//DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tcomb_TC_MAX)
//#pragma HLS DATAFLOW
//	float local_buf[OFM_BLmax];
//	uint16_t TBL_MIN_p[1];
//	uint32_t ofm_offset_p[1];
//
//		pp_cont( local_buf, output_buffer, bias_buffer,
//			TR_MIN, TC_MIN, OHxOW, ofm_w, IsNL, IsBias,
//			m, offset, burstlen, blmax_tc, T2_MIN,
//			TBL_MIN_p, ofm_offset_p, IsAllCont, IsColCont, tcb==0);
//
//		ofm_mmcpy_cont(Output, local_buf, ofm_offset_p[0], TBL_MIN_p[0]);
//	}
}

void pool_yolo2(float Input[Tn][OnChipIB_Height][OnChipIB_Width],float Output[Tm][Tr][Tc],
		uint8_t Ksize, uint8_t Kstride, uint16_t TM_MIN, uint16_t TR_MIN, uint16_t TC_MIN, bool enable)
{
	if(!enable)
		return;

	float tmp[Tn];
#pragma HLS ARRAY_PARTITION variable=tmp complete dim=1

	for(uint16_t tr = 0;tr < TR_MIN;tr++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tr)
	for(uint16_t tc = 0;tc < TC_MIN;tc++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=Tc)
	for(uint16_t i =0;i < Ksize; i++){
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=K)
	for(uint16_t j = 0;j < Ksize; j++)
	{
DO_PRAGMA(HLS LOOP_TRIPCOUNT min=1 max=K)
#pragma HLS PIPELINE II=1
		for(uint16_t of = 0; of < Tn; of++)
		{
			if(i==0&&j==0)
				tmp[of] = -1024*1024;

			if(Input[of][tr*Kstride+i][tc*Kstride+j] > tmp[of])
				tmp[of] = Input[of][tr*Kstride+i][tc*Kstride+j];

			if(i==1&&j==1)
				Output[of][tr][tc] = tmp[of];
		}
	}}}}
}

void load_compute_wrapper(float *ifm, float *weight, float ofm_buffer[Tm][Tr][Tc], int ksize, int kstride, int ifm_num, int ifm_w, int ifm_h,int ofm_num,
	 int ofm_h, int ofm_w, int pad_int, int ltype, int IHW, int KK, int INumxKK, int TM, int TN, int TR, int TC, int tm_r, int tr_r, int tc_r,
	 int tx_n1[3],int TX_MIN_n1[3],bool pp,bool in_flag,bool proc_flag, bool LoadBias, int NTif, uint8_t lmode, bool enable)
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

	static int tx_n[2][3];//tmtrtc, between load2compute
#pragma HLS ARRAY_PARTITION variable=tx_n complete dim=1
#pragma HLS ARRAY_PARTITION variable=tx_n complete dim=2
	static int TX_MIN_n[2][3];
#pragma HLS ARRAY_PARTITION variable=TX_MIN_n complete dim=1
#pragma HLS ARRAY_PARTITION variable=TX_MIN_n complete dim=2

	if(!enable)
		return ;

	int	TR_MIN = MIN(TR,ofm_h-tr_r);
	int	TC_MIN = MIN(TC,ofm_w-tc_r);
	int	TM_MIN = MIN(TM,ofm_num-tm_r);

	if(lmode==0){
		if(pp){
			if(ltype == 0)
			{
				copy_input_weight(ifm,weight,ifm_num,ifm_w,ifm_h,ksize,kstride,tr_r,tc_r,tm_r, 0,
					TM_MIN,TN,TR_MIN,TC_MIN,pad_int,ifm_buffer0,weight_buffer0,1,IHW,KK,INumxKK,ltype,in_flag);
				conv2d_tile(ifm_buffer1, ofm_buffer, weight_buffer1, 0, ksize, kstride, TX_MIN_n[1][0], TX_MIN_n[1][1], TX_MIN_n[1][2], proc_flag);
			}else if(ltype == 1)
			{
				copy_input_weight(ifm,weight,ifm_num,ifm_w,ifm_h,ksize,kstride,tr_r,tc_r,tm_r,tm_r,
										TM_MIN,TM,TR_MIN,TC_MIN,0,ifm_buffer0,weight_buffer0,0,IHW,KK,INumxKK,ltype,in_flag);
				pool_yolo2(ifm_buffer1, ofm_buffer, ksize, kstride, TX_MIN_n[1][0], TX_MIN_n[1][1], TX_MIN_n[1][2], proc_flag);
			}

			tx_n[0][0] = tm_r; tx_n[0][1] = tr_r; tx_n[0][2] = tc_r;
			TX_MIN_n[0][0] = TM_MIN; TX_MIN_n[0][1] = TR_MIN; TX_MIN_n[0][2] = TC_MIN;

			tx_n1[0] = tx_n[1][0]; tx_n1[1] = tx_n[1][1]; tx_n1[2] = tx_n[1][2];
			TX_MIN_n1[0] = TX_MIN_n[1][0]; TX_MIN_n1[1] = TX_MIN_n[1][1]; TX_MIN_n1[2] = TX_MIN_n[1][2];
		}else{
			if(ltype == 0)
			{
				copy_input_weight(ifm,weight,ifm_num,ifm_w,ifm_h,ksize,kstride,tr_r,tc_r,tm_r, 0,
					TM_MIN,TN,TR_MIN,TC_MIN,pad_int,ifm_buffer1,weight_buffer1,1,IHW,KK,INumxKK,ltype,in_flag);
				conv2d_tile(ifm_buffer0, ofm_buffer, weight_buffer0, 0, ksize, kstride, TX_MIN_n[0][0], TX_MIN_n[0][1], TX_MIN_n[0][2], proc_flag);
			}else if(ltype == 1)
			{
				copy_input_weight(ifm,weight,ifm_num,ifm_w,ifm_h,ksize,kstride,tr_r,tc_r,tm_r,tm_r,
										TM_MIN,TM,TR_MIN,TC_MIN,0,ifm_buffer1,weight_buffer1,0,IHW,KK,INumxKK,ltype,in_flag);
				pool_yolo2(ifm_buffer0, ofm_buffer, ksize, kstride, TX_MIN_n[0][0], TX_MIN_n[0][1], TX_MIN_n[0][2], proc_flag);
			}

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
						copy_input_weight(ifm,weight,ifm_num,ifm_w,ifm_h,ksize,kstride,tr_r,tc_r,tm_r, tn_r,
							TM_MIN,TN,TR_MIN,TC_MIN,pad_int,ifm_buffer0,weight_buffer0,1,IHW,KK,INumxKK,ltype,in_flag);
						conv2d_tile(ifm_buffer1, ofm_buffer, weight_buffer1, tn_n[1], ksize, kstride, TX_MIN_n[1][0], TX_MIN_n[1][1], TX_MIN_n[1][2], proc_flag);

						tn_n[0] = tn_r;
						tx_n[0][0] = tm_r; tx_n[0][1] = tr_r; tx_n[0][2] = tc_r;
						TX_MIN_n[0][0] = TM_MIN; TX_MIN_n[0][1] = TR_MIN; TX_MIN_n[0][2] = TC_MIN;
						tn_r += TN;

						pp_tn = 0;
					}else{
						copy_input_weight(ifm,weight,ifm_num,ifm_w,ifm_h,ksize,kstride,tr_r,tc_r,tm_r, tn_r,
							TM_MIN,TN,TR_MIN,TC_MIN,pad_int,ifm_buffer1,weight_buffer1,1,IHW,KK,INumxKK,ltype,in_flag);
						conv2d_tile(ifm_buffer0, ofm_buffer, weight_buffer0, tn_n[0], ksize, kstride, TX_MIN_n[0][0], TX_MIN_n[0][1], TX_MIN_n[0][2], proc_flag);

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

	const int KxK = Ksize*Ksize;
	const int IFM_numxKxK = IFM_num*KxK;

	const int IHxIW   = Input_h*Input_w;
	const int OHxOW = Output_h*Output_w;

	static float bias_buffer[MAX_BETA_LENGTH];
	static float ofm_buffer0[Tm][Tr][Tc];
#pragma HLS ARRAY_PARTITION variable=ofm_buffer0 complete dim=1
	static float ofm_buffer1[Tm][Tr][Tc];
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
				tx_n10, TX_MIN_n10, pp, in_flag, proc_flag, IsBias, NTif, lmode, lc_enable);

			write_back_output_reorg( ofm_buffer1, bias_buffer, Output, tx_n11[1], tx_n11[2], tx_n11[0], OFM_num, Output_h, Output_w, TX_MIN_n11[0], TX_MIN_n11[1], TX_MIN_n11[2], OHxOW, IsNL, IsBias, out_flag);
			pp = 0;
		}else{
			load_compute_wrapper(Input, Weight, ofm_buffer1, Ksize, Kstride, IFM_num, Input_w, Input_h, OFM_num, Output_h, Output_w,
				Padding, LayerType, IHxIW, KxK, IFM_numxKxK, TM, TN, TR, TC, tm_r, tr_r, tc_r,
				tx_n11, TX_MIN_n11, pp, in_flag, proc_flag, IsBias, NTif, lmode, lc_enable);

			write_back_output_reorg( ofm_buffer0, bias_buffer, Output, tx_n10[1], tx_n10[2], tx_n10[0], OFM_num, Output_h, Output_w, TX_MIN_n10[0], TX_MIN_n10[1], TX_MIN_n10[2], OHxOW, IsNL, IsBias, out_flag);
			pp = 1;
		}

		tm++;
		tm_r += TM;
		if(tm==NTof)
		{
			tm=0; tc++;
			tm_r = 0; tc_r += TC;
			if(tc==NTox)
			{
				tc=0; tr++;
				tc_r = 0; tr_r += TR;
		}}
	}
}
