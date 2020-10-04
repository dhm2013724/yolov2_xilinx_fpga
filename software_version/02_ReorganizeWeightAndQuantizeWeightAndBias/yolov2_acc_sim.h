

///////////////////////////////////////////////////////////////////////20181229 anti-reorg start => KxKxTmxTn
#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))
#define S 2
#define K 3

#define Tn 4
#define Tm 32
#define Tr 26
#define Tc 26
#define OnChipIB_Width  ((Tc-1)*S+K)
#define OnChipIB_Height ((Tr-1)*S+K)
#define MAX_BETA_LENGTH (1024)

//#define REORG_GEN
#define REORG_TEST

void input_load(float *input,float input_buffer[Tn][OnChipIB_Height][OnChipIB_Width],int r,int c,int n,int Kstride,int Padding,int TRow,int TCol,int Input_w,int Input_h,int TN_MIN,int IHxIW,int LayerType)
{
	int t1,t2,t3,t4;
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
		for(t2 = 0;t2 < TRow; t2++)
		{
			memcpy((float *)(input_memcpy_buffer + input_mmcpy_offset),(float *)(input + CurrentOffset + t1*IHxIW + t2*Input_w),TCol*sizeof(float));
			input_mmcpy_offset += TCol;
		}

	input_mmcpy_offset = 0;
	for(t1 = 0;t1 < Tn; t1++)
		for(t2 = 0;t2 < TRow; t2++)
			for(t3 = 0;t3 < TCol; t3++)
			{
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

void weight_load(float *Weight,float weight_buffer[Tm][Tn][K][K],bool weight_load_enable,int m,int n,int IFM_numxKxK,int KxK,int Ksize,int TM_MIN,int TN_MIN)
{
	int t1,t2,t3,t4;
	static float weight_memcpy_buffer[Tm*Tn*K*K];
	
	if(!weight_load_enable)
		return;
	
	const int Woffset = m*IFM_numxKxK + n*KxK;
		
	int weight_memcpy_offset = 0;
	for(t1 = 0;t1 < TM_MIN; t1++)
		for(t2 = 0;t2 < TN_MIN; t2++)
		{
			memcpy((float *)(weight_memcpy_buffer + weight_memcpy_offset),(float *)(Weight + Woffset + t1*IFM_numxKxK + t2*KxK),KxK*sizeof(float));
			weight_memcpy_offset += KxK;
		}
	
	weight_memcpy_offset = 0;
	for(t1 = 0;t1 < Tm; t1++)
		for(t2 = 0;t2 < Tn; t2++)
			for(t3 = 0;t3 <Ksize; t3++)
				for(t4 = 0;t4 <Ksize; t4++)
				{
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
	int t1,t2,t3,t4;
	static float weight_memcpy_buffer[Tm*Tn*K*K];
	static int Woffset;

	if(!weight_load_enable)
		return;

	if(m==0&&n==0)
		Woffset = 0;

	memcpy(weight_memcpy_buffer,(float *)(Weight + Woffset),TM_MIN*TN_MIN*KxK*sizeof(float));
	Woffset += TM_MIN*TN_MIN*KxK;
	
	int weight_memcpy_offset = 0;
	for(t3 = 0;t3 <Ksize; t3++)
		for(t4 = 0;t4 <Ksize; t4++)
			for(t1 = 0;t1 < Tm; t1++)
				for(t2 = 0;t2 < Tn; t2++)
				{
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
	//memcpy(local_beta_buffer,(float *)(beta_buffer+m),TM_MIN*sizeof(float));
	int offset;
	int tm;
	for(tm = 0,offset = m;tm < TM_MIN;tm++)
	{
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

	int i,j,tr,tc,tm,tn;
	int n = n_next[0];
	float partial_mul[Tm][Tn];
	float partial_add[Tm];

	for(i =0;i < Ksize; i++)
#pragma HLS LOOP_TRIPCOUNT min=1 max=5
		for(j = 0;j < Ksize; j++)
#pragma HLS LOOP_TRIPCOUNT min=1 max=5
			for(tr = 0;tr < TR_MIN;tr++)
#pragma HLS LOOP_TRIPCOUNT min=14 max=14
				for(tc = 0;tc < TC_MIN;tc++)
				{
#pragma HLS LOOP_TRIPCOUNT min=14 max=14
#pragma HLS PIPELINE
					for(tm = 0;tm < Tm;tm++)
					{
						if(i==0&&j==0&&n==0)
							partial_add[tm] = local_beta_buffer[tm];
						else
							partial_add[tm] = output_buffer[tm][tr][tc];
					}

					for(tm = 0;tm < Tm;tm++)
						for(tn = 0;tn <Tn;tn++)
						{
							partial_mul[tm][tn] = weight_buffer[tm][tn][i][j]*input_buffer[tn][Kstride*tr+i][Kstride*tc+j];
						}

					
					for(tm = 0;tm < Tm;tm++)
					{
						float partial_sum = 0;
						for(tn = 0;tn <Tn;tn++)
						{
							 partial_sum += partial_mul[tm][tn];
						}
						output_buffer[tm][tr][tc] = partial_add[tm] + partial_sum;
					}
				}
}

void nonlinear_leaky(float Input[Tm][Tr][Tc],const int TM_MIN,const int TR_MIN,const int TC_MIN,const bool IsNL)
{
	int tr,tc,tm;

	if(!IsNL)
		return ;
	
	for(tm = 0;tm < TM_MIN;tm++)
#pragma HLS LOOP_TRIPCOUNT min=1 max=1
		for(tr = 0;tr < TR_MIN;tr++)
#pragma HLS LOOP_TRIPCOUNT min=1 max=14
			for(tc = 0;tc < TC_MIN;tc++)
			{
#pragma HLS LOOP_TRIPCOUNT min=14 max=14
#pragma HLS PIPELINE
				float tmp = Input[tm][tr][tc];
				if(tmp < 0)
					Input[tm][tr][tc] = tmp*0.1;
			}

}

void write_back_output_reorg(float output_buffer[Tm][Tr][Tc],float *Output,int r,int c,int m,const int Output_w,const int Output_h,
					   const int TM_MIN,const int TR_MIN,const int TC_MIN,const int OHxOW, bool IsNL, bool write_flag)
{
	if(!write_flag)
		return;

	const int offset = m*OHxOW + r*Output_w + c;
	int tr,tm;

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

void pool_yolo2(float Input[Tn][OnChipIB_Height][OnChipIB_Width],float Output[Tm][Tr][Tc],
		  const int Ksize,const int Kstride,
		  const int TM_MIN,const int TR_MIN,const int TC_MIN,bool enable)
{
	if(!enable)
		return;

	int i,j,tr,tc,of;
	float tmp[Tn];

	for(tr = 0;tr < TR_MIN;tr++)
		for(tc = 0;tc < TC_MIN;tc++)
			for(i =0;i < Ksize; i++)
				for(j = 0;j < Ksize; j++)
				{
#pragma HLS PIPELINE
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
    	for( x = 0; x < TC_MIN; x++)
		for(ky= 0;ky < 2; ky++)
		for(kx = 0;kx < 2; kx++)
		{
#pragma HLS PIPELINE
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
				TM_MIN,TM,TRow,TCol,0,input_buffer0,NULL,NOP,input_flag,0,0,IHxIW,KxK,IFM_numxKxK,LayerType);
			pool_yolo2(input_buffer1,output_buffer,Ksize,Kstride,TX_MIN_next[0],TR_MIN,TC_MIN,process_flag);
		}else
		{
			TMP_X_next[0] = tmp_x;
			TX_MIN_next[0] = tmp_tx_min;
			tmp_x = TMP_M;
			tmp_tx_min = TM_MIN;

			copy_input_weight(Input,Weight,IFM_num,Input_w,Input_h,Ksize,Kstride,TMP_R,TMP_C,TMP_M,TMP_M,
				TM_MIN,TM,TRow,TCol,0,input_buffer1,NULL,NOP,input_flag,0,0,IHxIW,KxK,IFM_numxKxK,LayerType);
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
				TM_MIN,TM,TRow,TCol,0,input_buffer0,NULL,NOP,input_flag,0,0,IHxIW,KxK,IFM_numxKxK,LayerType);
			reorg_yolo2(input_buffer1,output_buffer,Ksize,Kstride,TX_MIN_next[0],TR_MIN,TC_MIN,process_flag);
		}else
		{
			TMP_X_next[0] = tmp_x;
			TX_MIN_next[0] = tmp_tx_min;
			tmp_x = TMP_M;
			tmp_tx_min = TM_MIN;

			copy_input_weight(Input,Weight,IFM_num,Input_w,Input_h,Ksize,Kstride,TMP_R,TMP_C,TMP_M,TMP_M,
				TM_MIN,TM,TRow,TCol,0,input_buffer1,NULL,NOP,input_flag,0,0,IHxIW,KxK,IFM_numxKxK,LayerType);
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
		TR_MIN = MIN(TR,Output_h-r);
		for(c = 0; c < Output_w; c += TC)
		{
			TC_MIN = MIN(TC,Output_w-c);
			pingpongm = 0;
			for(m = 0; m < OFM_num_bound; m += TM)
			{
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
									r, c, m, TM_MIN, TR_MIN, TC_MIN, TN, TRow, TCol, Padding,IHxIW,KxK,IFM_numxKxK,LayerType,TM, 										m1,TM_MIN1, pingpongm, input_flag, process_flag);

					write_back_output_reorg(output_buffer,Output, r, c, m0[0],Output_w,Output_h, TM_MIN0[0], TR_MIN, TC_MIN, OHxOW, IsNL, write_flag);
					pingpongm = 1;
				}else
				{
					intra_pingpong_wrapper(Input,Weight,output_buffer,beta_buffer,input_buffer0,input_buffer1,
									IFM_num, Input_w, Input_h, OFM_num, Ksize, Kstride,
									r, c, m, TM_MIN, TR_MIN, TC_MIN, TN, TRow, TCol, Padding,IHxIW,KxK,IFM_numxKxK,LayerType,TM, 										m0,TM_MIN0, pingpongm, input_flag, process_flag);

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
	int offset_index = 0;

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
	int woffset = 0;
	int boffset = 0;
	int TR,TC,TM,TN;
	int output_w,output_h;
	int mLoops;
	double sum_gop = 0.0;

    	for(i = 0; i < net->n; ++i)
	{
        	layer l = net->layers[i];
		printf("Layer[%2d]: ",i);
		switch(l.type)
		{
			case CONVOLUTIONAL:
				printf("outputMemory:%8d;BN=%d;Activation=%d;conv  %5d %2d x%2d /%2d  %4d x%4d x%4d   ->  %4d x%4d x%4d  %5.3f BFLOPs\n",l.outputs,l.batch_normalize,l.activation, l.n, l.size, l.size, l.stride, l.w, l.h, l.c, l.out_w, l.out_h, l.out_c, (2.0 * l.n * l.size*l.size*l.c/l.groups * l.out_h*l.out_w)/1000000000.);
				sum_gop += (2.0 * l.n * l.size*l.size*l.c/l.groups * l.out_h*l.out_w)/1000000000.;
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

				TR = MIN(((OnChipIB_Height-l.stride)/l.stride+1),Tr);//keep Kstride>=1
				TR = MIN(output_h,TR);
				TC = MIN(((OnChipIB_Width-l.stride)/l.stride+1),Tc);
				TC = MIN(output_w,TC);
				TM = 4;
				mLoops = 1;

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
	printf("SUM_GOP=%g\n",sum_gop);

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
