
#include "cnn.h"
// #include <string.h>
// #include <ap_int.h>

////////////////////////////////////////////2 ok start 100MHZ ???
#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))
#define S 2
#define K 3

#define Tn 2
#define Tm 32
#define Tr 26
#define Tc 26
#define OnChipIB_Width  ((Tc-1)*S+K)
#define OnChipIB_Height ((Tr-1)*S+K)
#define MAX_BETA_LENGTH (1024)
#define INTER_WIDTH (19)

////////////////////////////////////////////////No transfer padding start

void clear_input_buffer(short input_buffer[Tn][OnChipIB_Height][OnChipIB_Width],int TRow,int TCol)
{
#pragma HLS INLINE off
	ap_uint<6> t1,t2,t3;
	ap_uint<6> TRow_6b = TRow;
	ap_uint<6> TCol_6b = TCol;
		for(t2 = 0;t2 < TRow_6b; t2++)
			for(t3 = 0;t3 < TCol_6b; t3++)
			{
#pragma HLS PIPELINE
				for(t1 = 0;t1 < Tn; t1++)
				{
					input_buffer[t1][t2][t3] = 0;
				}
			}
}

void input_pixel_load(int *input,int input_memcpy_buffer[Tn*OnChipIB_Height*(OnChipIB_Width+2)/2],int r,int c,int n,int Kernel_stride,int Padding,int TRow,int TCol,int Input_w,int Input_h,int TN_MIN,int IHxIW,
					  	int row_len_o[1],int col_len_o[1],int RowSub[1],int ColSub[1],ap_uint<1> RowBeginByte[Tn*OnChipIB_Height])
{
	ap_uint<6> t1,t2;

	ap_uint<2> Kernel_stride_2b = Kernel_stride;
	ap_uint<1> Padding_1b = Padding;
	ap_uint<9> Input_w_9b = Input_w;
	ap_uint<9> Input_h_9b = Input_h;
	ap_uint<18> IHxIW_18b = IHxIW;
	ap_uint<6> TN_MIN_6b = TN_MIN;
	ap_uint<9> c_9b = c;
	ap_uint<9> r_9b = r;
	ap_uint<11> n_11b = n;
	ap_uint<6> TCol_6b = TCol;
	ap_uint<6> TRow_6b = TRow;

	ap_uint<1> shiftbit;
	if(Kernel_stride_2b==1)
		shiftbit = 0;
	else
		shiftbit = 1;

	ap_int<13> Coffset = (c_9b << shiftbit) - Padding_1b;
	ap_int<13> Roffset = (r_9b << shiftbit) - Padding_1b;

	ap_uint<9> TRow_bottom,TCol_right;
	ap_int<13> TRow_top,TCol_left;
	ap_uint<6> row_len,col_len;

	if(Coffset > 0)
		TCol_left = Coffset;
	else
		TCol_left = 0;

	if((Coffset + TCol_6b-1)<Input_w_9b)
		TCol_right = Coffset + TCol_6b;
	else
		TCol_right = Input_w_9b;

	col_len = TCol_right - TCol_left;

	if(Roffset > 0)
		TRow_top = Roffset;
	else
		TRow_top = 0;

	if((Roffset + TRow_6b-1)<Input_h_9b)
		TRow_bottom = Roffset + TRow_6b;
	else
		TRow_bottom = Input_h_9b;

	row_len = TRow_bottom - TRow_top;

	int IN_OFFSET = n_11b*IHxIW_18b + TRow_top*Input_w_9b +TCol_left;

	int RowIntNum[Tn*OnChipIB_Height];
	int RowOffset[Tn*OnChipIB_Height];

	int InOffset;
	ap_uint<1> LowBit;// 0 or 1
	int BeginByteNum;

	bool IsContinue = (Input_w_9b <= 26);

//	int ColNum = IsContinue ? (row_len*col_len) : col_len;
//	int RowNum = IsContinue ? 1 : row_len;
	ap_uint<12> ColNum;
	ap_uint<6> RowNum;

	if(IsContinue)
	{
		ColNum = row_len*col_len;
		RowNum = 1;
	}else
	{
		ColNum = col_len;
		RowNum = row_len;
	}

	int RowBeginByteIndex = 0;
	for(t1 = 0;t1 < TN_MIN_6b; t1++)
		for(t2 = 0;t2 < RowNum; t2++)
		{
#pragma HLS PIPELINE
			InOffset = IN_OFFSET + t1*IHxIW_18b + t2*Input_w_9b;//maybe <=0
			LowBit = InOffset&0x1;
			BeginByteNum = ColNum + LowBit;
			RowIntNum[RowBeginByteIndex] = BeginByteNum >> 1;
			if(BeginByteNum&0x1)
				RowIntNum[RowBeginByteIndex]++;
			RowBeginByte[RowBeginByteIndex] = LowBit;
			RowOffset[RowBeginByteIndex] = InOffset >> 1;//must shift bucause of <=0
			RowBeginByteIndex++;
		}

	int RowBeginByteIndex2 = 0;
	int input_mmcpy_offset = 0;
	for(t1 = 0;t1 < TN_MIN_6b; t1++)
		for(t2 = 0;t2 < RowNum; t2++)
		{
			memcpy((int *)(input_memcpy_buffer + input_mmcpy_offset),(int *)(input + RowOffset[RowBeginByteIndex2]),RowIntNum[RowBeginByteIndex2]*sizeof(int));
			input_mmcpy_offset += RowIntNum[RowBeginByteIndex2];
			RowBeginByteIndex2++;
		}

	RowSub[0] = TRow_top - Roffset;
	ColSub[0] = TCol_left - Coffset;
	row_len_o[0] = row_len;
	col_len_o[0] = col_len;
}

void copy_input2buf(short input_buffer[Tn][OnChipIB_Height][OnChipIB_Width],int Input_w,int TN_MIN,int row_len[1],int col_len[1],int RowSub[1],int ColSub[1],
					int input_memcpy_buffer[Tn*OnChipIB_Height*(OnChipIB_Width+2)/2],ap_uint<1> RowBeginByte[Tn*OnChipIB_Height])
{
	ap_uint<6> t1,t2,t3;
	bool IsContinue = (Input_w <= 26);
	bool BarIsContinue = !IsContinue;

	int input_mmcpy_offset = 0;
	int RowBeginByteIndex = 0;
	bool NextInputFlag;
	ap_uint<1> cnt;
	short input_array[2];
#pragma HLS ARRAY_PARTITION variable=input_array complete dim=1

	ap_uint<6> TN_MIN_6b = TN_MIN;
	ap_uint<6> row_len_6b = row_len[0];
	ap_uint<6> col_len_6b = col_len[0];
	ap_uint<2> RowSub_2b = RowSub[0];
	ap_uint<2> ColSub_2b = ColSub[0];

	for(t1 = 0;t1 < TN_MIN_6b; t1++)
		for(t2 = 0;t2 < row_len_6b; t2++)
			for(t3 = 0;t3 < col_len_6b; t3++)
			{
#pragma HLS PIPELINE
				bool T3e0 = (t3==0);
				bool T2e0 = (t2==0);

				if(((IsContinue&&T2e0)||BarIsContinue)&&T3e0)
				{
					cnt = RowBeginByte[RowBeginByteIndex];
					RowBeginByteIndex++;
					NextInputFlag = true;
				}

				if(NextInputFlag)
				{
					input_array[0] = input_memcpy_buffer[input_mmcpy_offset];
					input_array[1] = input_memcpy_buffer[input_mmcpy_offset] >> 16;
					input_mmcpy_offset++;
					NextInputFlag = false;
				}

				input_buffer[t1][t2 + RowSub_2b][t3 + ColSub_2b] = input_array[cnt];

				if(cnt==1)
				{
					NextInputFlag = true;
					cnt = 0;
				}else
				{
					cnt = 1;
				}
			}
}

void input_load(int *input,short input_buffer[Tn][OnChipIB_Height][OnChipIB_Width],int r,int c,int n,int Kernel_stride,int Padding,int TRow,int TCol,int Input_w,int Input_h,int TN_MIN,int IHxIW)
{
	int row_len[1],col_len[1];
	int RowSub[1],ColSub[1];

	static int input_memcpy_buffer[Tn*OnChipIB_Height*(OnChipIB_Width+2)/2];
	ap_uint<1> RowBeginByte[Tn*OnChipIB_Height];//0 ro 1

	clear_input_buffer(input_buffer, TRow, TCol);
	input_pixel_load(input,input_memcpy_buffer, r, c, n, Kernel_stride, Padding, TRow, TCol, Input_w, Input_h, TN_MIN, IHxIW,
					  	row_len,col_len,RowSub,ColSub, RowBeginByte);
	copy_input2buf( input_buffer, Input_w, TN_MIN, row_len, col_len, RowSub, ColSub,input_memcpy_buffer, RowBeginByte);

}
////////////////////////////////////////////////No transfer padding end

void weight_load_reorg(int *Weight,short weight_buffer[Tm][Tn][K][K],bool weight_load_enable,int m,int n,int IFM_numxKxK,int KxK,int Kernel_size,int TM_MIN,int TN_MIN)
{
	ap_uint<6> t1,t2;
	ap_uint<2> t3,t4;
	static int weight_memcpy_buffer[Tm*Tn*K*K/2];
	static int Woffset;

	if(!weight_load_enable)
		return;

	if(m==0&&n==0)
		Woffset = 0;
//	if((TM_MIN*TN_MIN*KxK)%2)
//		printf("weight % error\n");

	ap_uint<6> TM_MIN_6b = TM_MIN;
	ap_uint<6> TN_MIN_6b = TN_MIN;
	ap_uint<4> KxK_4b = KxK;
	ap_uint<2> Kernel_size_2b = Kernel_size;

	int ReadBurstByteNum = (TM_MIN_6b*TN_MIN_6b*KxK_4b) >> 1;
	memcpy(weight_memcpy_buffer,(int *)(Weight + Woffset),ReadBurstByteNum*sizeof(int));
	Woffset += ReadBurstByteNum;

	int weight_memcpy_offset = 0;
	bool ReadNextFlag = true;
	short next_value;
	short input_value;
	for(t1 = 0;t1 < Tm; t1++)
		for(t2 = 0;t2 < Tn; t2++)
			for(t3 = 0;t3 <Kernel_size_2b; t3++)
				for(t4 = 0;t4 <Kernel_size_2b; t4++)
				{
#pragma HLS PIPELINE
					bool Enable = (t1 < TM_MIN_6b)&&(t2 < TN_MIN_6b);
					if(Enable)
					{
						if(ReadNextFlag)
						{
							input_value = weight_memcpy_buffer[weight_memcpy_offset];
							next_value = weight_memcpy_buffer[weight_memcpy_offset] >> 16;
							weight_memcpy_offset++;
							ReadNextFlag = false;
						}else
						{
							input_value = next_value;
							ReadNextFlag = true;
						}

						weight_buffer[t1][t2][t3][t4] =  input_value;
					}
					else
						weight_buffer[t1][t2][t3][t4] = 0;
				}
}

void copy_input_weight(int *input,int *Weight,int InFM_num,int Input_w,int Input_h,int OutFM_num,int Kernel_size,int Kernel_stride,int r,int c,int m,int n,
		int TM_MIN,int TN,int TRow,int TCol,int Padding,short input_buffer[Tn][OnChipIB_Height][OnChipIB_Width],short weight_buffer[Tm][Tn][K][K],int TMP_N_next[1],
		bool enable,bool weight_load_enable,bool initialize,const int IHxIW,const int KxK,const int IFM_numxKxK)
{
	if(!enable)
		return ;

	const int TN_MIN = MIN(TN,InFM_num - n);
	TMP_N_next[0] = n;

	input_load(input, input_buffer, r, c, n, Kernel_stride, Padding, TRow, TCol, Input_w, Input_h, TN_MIN, IHxIW);

	weight_load_reorg(Weight,weight_buffer,weight_load_enable,m,n,IFM_numxKxK,KxK,Kernel_size,TM_MIN,TN_MIN);

}
//////////////////////////////////////////////////T3 end
	void copy_local_beta(short beta_buffer[MAX_BETA_LENGTH],short local_beta_buffer[Tm],const int TM_MIN,int m)
{
	//memcpy(local_beta_buffer,(float *)(beta_buffer+m),TM_MIN*sizeof(float));
	int offset;
	int tm;
	for(tm = 0,offset = m;tm < TM_MIN;tm++)
	{
#pragma HLS PIPELINE
		local_beta_buffer[tm] = beta_buffer[offset];
		offset++;
	}
}

void nonlinear_leaky(int Input[Tm][Tr][Tc],short Output[Tm][Tr][Tc],const int TM_MIN,const int TR_MIN,const int TC_MIN,const bool IsNL,const int InterSubOutput,bool postenable)
{
	if(!postenable)
		return;

	ap_uint<6> tr,tc,tm;
	int tmp_output;
	ap_uint<6> TM_MIN_6b = TM_MIN;
	ap_uint<6> TR_MIN_6b = TR_MIN;
	ap_uint<6> TC_MIN_6b = TC_MIN;

	for(tm = 0;tm < TM_MIN_6b;tm++)
		for(tr = 0;tr < TR_MIN_6b;tr++)
			for(tc = 0;tc < TC_MIN_6b;tc++)
			{
#pragma HLS PIPELINE
				int tmp = Input[tm][tr][tc];
				if(IsNL&&tmp<0)
				{
					tmp_output = (((long long)tmp*0xccc)>>15);
				}else
				{
					tmp_output = tmp;
				}
				Output[tm][tr][tc] = tmp_output>>InterSubOutput;
			}

}

void compute(short input_buffer[Tn][OnChipIB_Height][OnChipIB_Width],short output_buffer[Tm][Tr][Tc],
		short weight_buffer[Tm][Tn][K][K],short beta_buffer[MAX_BETA_LENGTH],int TMP_N_next[1],
		const int Kernel_size,const int Kernel_stride,int TMP_M,
		const int TM_MIN,const int TR_MIN,const int TC_MIN,bool enable,const bool IsNL,const bool reluenable,
		const int InterSubBeta,const int WeightAddInputSubInter,const int InterSubOutput)
{
	static short local_beta_buffer[Tm];
#pragma HLS ARRAY_PARTITION variable=local_beta_buffer complete dim=1

	static int compute_buffer[Tm][Tr][Tc];
#pragma HLS ARRAY_PARTITION variable=compute_buffer complete dim=1

	if(!enable)
	{
		copy_local_beta(beta_buffer,local_beta_buffer,TM_MIN,TMP_M);
		return;
	}

	int partial_mul[Tm][Tn];
#pragma HLS ARRAY_PARTITION variable=partial_mul complete dim=1
#pragma HLS ARRAY_PARTITION variable=partial_mul complete dim=2

	ap_uint<2> Kernel_size_2b = Kernel_size;
	ap_uint<2> i,j;
	ap_uint<6> TR_MIN_6b = TR_MIN;
	ap_uint<6> TC_MIN_6b = TC_MIN;
	ap_uint<6> tr,tc,tm,tn;

	int n = TMP_N_next[0];
	for(i = 0;i < Kernel_size_2b; i++)
		for(j = 0;j < Kernel_size_2b; j++)
			for(tr = 0;tr < TR_MIN_6b;tr++)
				for(tc = 0;tc < TC_MIN_6b;tc++)
				{
#pragma HLS PIPELINE
					for(tm = 0;tm < Tm;tm++)
					{
#pragma HLS DEPENDENCE variable=compute_buffer inter false
						int tmp_add_result;
						if(i==0&&j==0&&n==0)
						{
							tmp_add_result = local_beta_buffer[tm] << InterSubBeta;
						}
						else
							tmp_add_result = compute_buffer[tm][tr][tc];

						partial_mul[tm][0] = (weight_buffer[tm][0][i][j]*input_buffer[0][tr+i][tc+j]) >> WeightAddInputSubInter;//Q1+Q2-Q3
						partial_mul[tm][1] = (weight_buffer[tm][1][i][j]*input_buffer[1][tr+i][tc+j]) >> WeightAddInputSubInter;//Q1+Q2-Q3

						compute_buffer[tm][tr][tc] = tmp_add_result + partial_mul[tm][0] + partial_mul[tm][1];
					}
				}

	nonlinear_leaky(compute_buffer,output_buffer,TM_MIN,TR_MIN,TC_MIN,IsNL,InterSubOutput,reluenable);

}

void write_back_output_reorg(short output_buffer[Tm][Tr][Tc],int *Output,int r,int c,int m,const int Output_w,const int Output_h,
					   const int TM_MIN,const int TR_MIN,const int TC_MIN,const int OHxOW,bool write_flag,const int OutputQ)
{
	static int output_tmp[Tm*Tr*Tc/2];
	ap_uint<6> tr,tm,tc;

	if(!write_flag)
		return;

	ap_uint<6> TM_MIN_6b = TM_MIN;
	ap_uint<6> TR_MIN_6b = TR_MIN;
	ap_uint<6> TC_MIN_6b = TC_MIN;
	ap_uint<9> Output_w_9b = Output_w;
	ap_uint<18> OHxOW_18b = OHxOW;

	const int offset = m*OHxOW_18b + r*Output_w_9b + c;

	ap_uint<6> TM_MIN_g;
	if(TM_MIN_6b==9)
		TM_MIN_g = 12;
	else
		TM_MIN_g = TM_MIN_6b;

	char cnt = 0;
	int outputoffset = 0;
	short ouput_array[2];
#pragma HLS ARRAY_PARTITION variable=ouput_array complete dim=1
	for(tm = 0;tm < TM_MIN_g;tm++)
		for(tr = 0;tr < TR_MIN_6b;tr++)
			for(tc = 0;tc < TC_MIN_6b;tc++)
			{
#pragma HLS PIPELINE
				ouput_array[cnt] = output_buffer[tm][tr][tc];
				cnt++;
				if(cnt==2)
				{
					output_tmp[outputoffset] = (ouput_array[0]       &0x0000FFFF) |
											  ((ouput_array[1] << 16 )&0xFFFF0000);
					outputoffset++;
					cnt = 0;
				}
			}

	int OutputLength;
	ap_uint<6> Loop1,Loop2;
	ap_uint<18> OutputOffset1;
	ap_uint<9> OutputOffset2;
	ap_uint<10> OtuputTmpOffset1;
	ap_uint<6> OtuputTmpOffset2;

	if(TC_MIN_6b == 26)
	{
		Loop1 = TM_MIN_6b;
		Loop2 = 26;
		OutputOffset1 = OHxOW_18b;
		OutputOffset2 = Output_w_9b;
		OtuputTmpOffset1 = 26*26;
		OtuputTmpOffset2 = 26;
		OutputLength = 26/2;
	}else//TMxTRxTC TMx13x13 continues
	{
		Loop1 = 1;
		Loop2 = 1;
		OutputOffset1 = 0;
		OutputOffset2 = 0;
		OtuputTmpOffset1 = 0;
		OtuputTmpOffset2 = 0;
		if(TM_MIN_6b==9)
			OutputLength = 6*13*13;
		else
			OutputLength = (TM_MIN_6b >> 1)*169;
	}

	int OutputOffset1_sum;
	int OutputOffset2_sum;
	int OtuputTmpOffset1_sum;
	int OtuputTmpOffset2_sum;
	for(tm = 0,OutputOffset1_sum = 0,OtuputTmpOffset1_sum = 0;tm < Loop1;tm++,OutputOffset1_sum += OutputOffset1, OtuputTmpOffset1_sum += OtuputTmpOffset1)
		for(tr = 0,OutputOffset2_sum = 0,OtuputTmpOffset2_sum = 0;tr < Loop2;tr++,OutputOffset2_sum += OutputOffset2, OtuputTmpOffset2_sum += OtuputTmpOffset2)
		{
			memcpy((int *)(Output + (OutputOffset1_sum + OutputOffset2_sum + offset)/2),(int *)(output_tmp + (OtuputTmpOffset1_sum + OtuputTmpOffset2_sum)/2),OutputLength*sizeof(int));
		}
}

void pool_yolo2(short Input[Tn][OnChipIB_Height][OnChipIB_Width],short Output[Tm][Tr][Tc],
		  const int Kernel_size,const int Kernel_stride,
		  const int TM_MIN,const int TR_MIN,const int TC_MIN,bool enable)
{
	int x, y, of;
	int kx,ky;

	int Yoffset;
	int Xoffset;
	short tmp[Tn];
#pragma HLS ARRAY_PARTITION variable=tmp complete dim=1
	int tr,tc,tm;

	if(!enable)
		return;

    for( y = 0; y < TR_MIN; y++)
    	for( x = 0; x < TC_MIN; x++)
			for(Yoffset = y << 1,Xoffset = x << 1,ky= 0;ky < 2; ky++)
    			for(kx = 0;kx < 2; kx++)
				{
#pragma HLS PIPELINE
					for( of = 0; of < Tn; of++)
					{
						if(kx==0&&ky==0)
							tmp[of] = 0x8001;

						if(Input[of][Yoffset+ky][Xoffset+kx]>tmp[of])
							tmp[of] = Input[of][Yoffset+ky][Xoffset+kx];

						if(kx==1&&ky==1)
							Output[of][y][x] = tmp[of];
					}
    			}
}

void reorg_yolo2(short Input[Tn][OnChipIB_Height][OnChipIB_Width],short Output[Tm][Tr][Tc],
		  const int Kernel_size,const int Kernel_stride,
		  const int TM_MIN,const int TR_MIN,const int TC_MIN,bool enable)
{
	int x, y,kx,ky;
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

void intra_pingpong_wrapper(int *Input,int *Weight, short output_buffer[Tm][Tr][Tc],short beta_buffer[MAX_BETA_LENGTH],
								 short input_buffer0[Tn][OnChipIB_Height][OnChipIB_Width],short input_buffer1[Tn][OnChipIB_Height][OnChipIB_Width],
								 int InFM_num,int Input_w,int Input_h,int OutFM_num,int Kernel_size,int Kernel_stride,
								 int TMP_R,int TMP_C,int TMP_M,int m,int TM_MIN,int TR_MIN,int TC_MIN,int TN,int TRow,int TCol,int Padding,
								 int IHxIW,int KxK,int IFM_numxKxK,int nLoops,bool IsNL,int LayerType,int TM,int TMP_X_next[1],int TX_MIN_next[1],bool pingpongx,bool input_flag,bool process_flag,
								 const int InterSubBeta,const int WeightAddInputSubInter,const int InterSubOutput)
{
	static short weight_buffer0[Tm][Tn][K][K];
#pragma HLS ARRAY_PARTITION variable=weight_buffer0 complete dim=1
#pragma HLS ARRAY_PARTITION variable=weight_buffer0 complete dim=2

	static short weight_buffer1[Tm][Tn][K][K];
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
		int TMP_N_next0[1];
		int TMP_N_next1[1];
		int n;
		int TMP_N;
		for(TMP_N = 0,n = 0;n < nLoops+1; n++,TMP_N += TN)
		{
			if(pingpong == 1)
			{
				copy_input_weight(Input,Weight,InFM_num,Input_w,Input_h,OutFM_num,Kernel_size,Kernel_stride,TMP_R,TMP_C,TMP_M,TMP_N,
					TM_MIN,TN,TRow,TCol,Padding,input_buffer1,weight_buffer1,TMP_N_next1,n!=nLoops,1,(m==0)&&(n==0),IHxIW,KxK,IFM_numxKxK);
				compute(input_buffer0,output_buffer,weight_buffer0,beta_buffer,TMP_N_next0,Kernel_size,Kernel_stride,TMP_M,TM_MIN,TR_MIN,TC_MIN,n!=0,IsNL,n==nLoops,
					InterSubBeta,WeightAddInputSubInter,InterSubOutput);
				pingpong = 0;
			}else
			{
				copy_input_weight(Input,Weight,InFM_num,Input_w,Input_h,OutFM_num,Kernel_size,Kernel_stride,TMP_R,TMP_C,TMP_M,TMP_N,
					TM_MIN,TN,TRow,TCol,Padding,input_buffer0,weight_buffer0,TMP_N_next0,n!=nLoops,1,(m==0)&&(n==0),IHxIW,KxK,IFM_numxKxK);
				compute(input_buffer1,output_buffer,weight_buffer1,beta_buffer,TMP_N_next1,Kernel_size,Kernel_stride,TMP_M,TM_MIN,TR_MIN,TC_MIN,n!=0,IsNL,n==nLoops,
					InterSubBeta,WeightAddInputSubInter,InterSubOutput);
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

			copy_input_weight(Input,Weight,InFM_num,Input_w,Input_h,OutFM_num,Kernel_size,Kernel_stride,TMP_R,TMP_C,TMP_M,TMP_M,
				TM_MIN,TM,TRow,TCol,0,input_buffer0,weight_buffer0,NOP,input_flag,0,0,IHxIW,KxK,IFM_numxKxK);
			pool_yolo2(input_buffer1,output_buffer,Kernel_size,Kernel_stride,TX_MIN_next[0],TR_MIN,TC_MIN,process_flag);
		}else
		{
			TMP_X_next[0] = tmp_x;
			TX_MIN_next[0] = tmp_tx_min;
			tmp_x = TMP_M;
			tmp_tx_min = TM_MIN;

			copy_input_weight(Input,Weight,InFM_num,Input_w,Input_h,OutFM_num,Kernel_size,Kernel_stride,TMP_R,TMP_C,TMP_M,TMP_M,
				TM_MIN,TM,TRow,TCol,0,input_buffer1,weight_buffer0,NOP,input_flag,0,0,IHxIW,KxK,IFM_numxKxK);
			pool_yolo2(input_buffer0,output_buffer,Kernel_size,Kernel_stride,TX_MIN_next[0],TR_MIN,TC_MIN,process_flag);
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

			copy_input_weight(Input,Weight,InFM_num,Input_w,Input_h,OutFM_num,Kernel_size,Kernel_stride,TMP_R,TMP_C,TMP_M,TMP_M,
				TM_MIN,TM,TRow,TCol,0,input_buffer0,weight_buffer0,NOP,input_flag,0,0,IHxIW,KxK,IFM_numxKxK);
			reorg_yolo2(input_buffer1,output_buffer,Kernel_size,Kernel_stride,TX_MIN_next[0],TR_MIN,TC_MIN,process_flag);
		}else
		{
			TMP_X_next[0] = tmp_x;
			TX_MIN_next[0] = tmp_tx_min;
			tmp_x = TMP_M;
			tmp_tx_min = TM_MIN;

			copy_input_weight(Input,Weight,InFM_num,Input_w,Input_h,OutFM_num,Kernel_size,Kernel_stride,TMP_R,TMP_C,TMP_M,TMP_M,
				TM_MIN,TM,TRow,TCol,0,input_buffer1,weight_buffer0,NOP,input_flag,0,0,IHxIW,KxK,IFM_numxKxK);
			reorg_yolo2(input_buffer0,output_buffer,Kernel_size,Kernel_stride,TX_MIN_next[0],TR_MIN,TC_MIN,process_flag);
		}

	}

}

void copy_beta(short beta_buffer[MAX_BETA_LENGTH],int *Beta,const int OFM_NUM)
{
	static int beta_tmp[MAX_BETA_LENGTH/2];
	int NUM = (OFM_NUM+1)/2;
	memcpy(beta_tmp,(int *)Beta,NUM*sizeof(int));
	int x;
	for(x = 0;x < NUM;x++)
	{
#pragma HLS PIPELINE
		beta_buffer[2*x]   = beta_tmp[x];
		beta_buffer[2*x+1] = beta_tmp[x]>>16;
	}
}

void YOLO2_FPGA(int *Input,int *Output,int *Weight,int *Beta,const int InFM_num,const int OutFM_num,
							  const int Kernel_size,const int Kernel_stride,
							  const int Input_w,const int Input_h,const int Padding,const bool IsNL,const bool IsBN,
							  const int TM,const int TN,const int TR,const int TC,
							  const int mLoops,const int nLoops,const int rLoops,const int cLoops,const int LayerType,
							  const int InputQ,const int OutputQ,const int WeightQ,const int BetaQ)
{
#pragma HLS INTERFACE m_axi depth=512 port=Input offset=slave bundle=DATA_BUS num_read_outstanding=32 num_write_outstanding=32 max_read_burst_length=64 max_write_burst_length=64
#pragma HLS INTERFACE m_axi depth=512 port=Output offset=slave bundle=DATA_BUS num_read_outstanding=32 num_write_outstanding=32 max_read_burst_length=64 max_write_burst_length=64
#pragma HLS INTERFACE m_axi depth=512 port=Weight offset=slave bundle=DATA_BUS2 num_read_outstanding=32 max_read_burst_length=64
#pragma HLS INTERFACE m_axi depth=512 port=Beta offset=slave bundle=DATA_BUS2 num_read_outstanding=32 max_read_burst_length=64

#pragma HLS INTERFACE s_axilite register port=return bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=InFM_num bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=OutFM_num bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=Kernel_size bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=Kernel_stride bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=Input_w bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=Input_h bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=Padding bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=IsNL bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=IsBN bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=TM bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=TN bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=TR bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=TC bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=mLoops bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=nLoops bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=rLoops bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=cLoops bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=LayerType bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=InputQ bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=OutputQ bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=WeightQ bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=BetaQ bundle=CTRL_BUS

#pragma HLS INTERFACE s_axilite register port=Input bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=Output bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=Weight bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite register port=Beta bundle=CTRL_BUS

	const int output_w = (Input_w - Kernel_size + 2*Padding)/Kernel_stride + 1 ;
	const int output_h = (Input_h - Kernel_size + 2*Padding)/Kernel_stride + 1 ;
	const int OHxOW = output_h*output_w;
	const int TRow = (TR-1)*Kernel_stride+Kernel_size;
	const int TCol = (TC-1)*Kernel_stride+Kernel_size;
	const int IHxIW   = Input_h*Input_w;
	const int KxK = Kernel_size*Kernel_size;
	const int IFM_numxKxK = InFM_num*KxK;
	const int mLoops_bound = LayerType ? (mLoops + 2): (mLoops + 1);
	const int InterSubOutput = INTER_WIDTH-OutputQ;
	const int InterSubBeta = INTER_WIDTH - BetaQ;
	const int WeightAddInputSubInter = WeightQ+InputQ-INTER_WIDTH;

	static short input_buffer0[Tn][OnChipIB_Height][OnChipIB_Width];
#pragma HLS ARRAY_PARTITION variable=input_buffer0 complete dim=1

	static short input_buffer1[Tn][OnChipIB_Height][OnChipIB_Width];
#pragma HLS ARRAY_PARTITION variable=input_buffer1 complete dim=1

	static short output_buffer[Tm][Tr][Tc];
#pragma HLS ARRAY_PARTITION variable=output_buffer complete dim=1

	static short output_buffer1[Tm][Tr][Tc];
#pragma HLS ARRAY_PARTITION variable=output_buffer1 complete dim=1

	static short beta_buffer[MAX_BETA_LENGTH];

	int r,c,m;
/////////////////////////////////param
	int TMP_R,TMP_C,TMP_M;
	int TM_MIN,TR_MIN,TC_MIN;
///////////////////////////////////////

	int TMP_M_next0[1];
	int TMP_M_next1[1];
	int TM_MIN_next0[1];
	int TM_MIN_next1[1];
	bool pingpongm;

	if(LayerType==0)
		copy_beta(beta_buffer,Beta,OutFM_num);

	for(TMP_R = 0,r = 0; r < rLoops; r++, TMP_R += TR)
	{
		TR_MIN = MIN(TR,output_h -TMP_R);
		for(TMP_C = 0,c = 0; c < cLoops; c++,TMP_C += TC)
		{
			TC_MIN = MIN(TC,output_w -TMP_C);
			pingpongm = 0;
			for(TMP_M = 0, m = 0; m < mLoops_bound; m++,TMP_M += TM)
			{
				TM_MIN = MIN(TM,OutFM_num-TMP_M);
				bool MneZero = (m!=0);
				bool MneOne = (m!=1);
				bool MnemLoops = (m!=mLoops);
				bool MneMLoopsaddOne = (m!=(mLoops+1));
				bool input_flag = LayerType ? MnemLoops&&MneMLoopsaddOne: MnemLoops;
				bool process_flag = LayerType ? MneZero&&MneMLoopsaddOne : MnemLoops;
				bool write_flag = LayerType ? MneZero&&MneOne : MneZero;

				if(pingpongm==0)
				{
					intra_pingpong_wrapper(Input,Weight,output_buffer1,beta_buffer,input_buffer0,input_buffer1,
									InFM_num, Input_w, Input_h, OutFM_num, Kernel_size, Kernel_stride,
									TMP_R, TMP_C, TMP_M, m, TM_MIN, TR_MIN, TC_MIN, TN, TRow, TCol, Padding,IHxIW,KxK,IFM_numxKxK,nLoops,IsNL,LayerType,TM, TMP_M_next1,TM_MIN_next1, pingpongm, input_flag, process_flag,
									InterSubBeta,WeightAddInputSubInter,InterSubOutput);

					write_back_output_reorg(output_buffer,Output,TMP_R,TMP_C,TMP_M_next0[0],output_w,output_h,TM_MIN_next0[0],TR_MIN,TC_MIN,OHxOW,write_flag,OutputQ);
					pingpongm = 1;
				}else
				{
					intra_pingpong_wrapper(Input,Weight,output_buffer,beta_buffer,input_buffer0,input_buffer1,
									InFM_num, Input_w, Input_h, OutFM_num, Kernel_size, Kernel_stride,
									TMP_R, TMP_C, TMP_M, m, TM_MIN, TR_MIN, TC_MIN, TN, TRow, TCol, Padding,IHxIW,KxK,IFM_numxKxK,nLoops,IsNL,LayerType,TM, TMP_M_next0,TM_MIN_next0, pingpongm, input_flag, process_flag,
									InterSubBeta,WeightAddInputSubInter,InterSubOutput);

					write_back_output_reorg(output_buffer1,Output,TMP_R,TMP_C,TMP_M_next1[0],output_w,output_h,TM_MIN_next1[0],TR_MIN,TC_MIN,OHxOW,write_flag,OutputQ);
					pingpongm = 0;
				}

			}
		}
	}
}
////////////////////////////////////////////////////2 ok end ???


