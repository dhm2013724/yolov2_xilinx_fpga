
void weight_load_tile_buf_reorder(float *weight, float *weight_reorder, int INumxKK, int KK, int TM_MIN, int TN_MIN, int tm, int tn, int ksize, int TN4type)//assume w equal h
{
	static float local_wbuf[Tm*Tn*K*K];
	int t1,t2,t3,t4;
	static int base_offset;

	if((tm==0)&&(tn==0))
		base_offset = 0;	

	int offset = 0;
	for(t1 = 0;t1 < TM_MIN; t1++)
	for(t2 = 0;t2 < TN_MIN; t2++)
	for(t3 = 0;t3 <ksize; t3++)
	for(t4 = 0;t4 <ksize; t4++)
	{
		local_wbuf[offset++] =  weight[(t1+tm)*INumxKK + (t2+tn)*KK + t3*ksize + t4];
	}

	memcpy(weight_reorder + base_offset, local_wbuf, TM_MIN*TN_MIN*KK*sizeof(float));
	base_offset += TM_MIN*TN_MIN*KK;
}

void Reorder_weight(float *weight, float *weight_reorder, int ksize, int ifm_num, int ofm_num, int ltype, int TM, int TN)
{
	const int KK = ksize*ksize;
	const int INumxKK = ifm_num*KK;

	int tm, tn;
	int TM_MIN, TN_MIN;


    for(tm = 0; tm < ofm_num; tm += TM)
    {
	TM_MIN = MIN_diy(TM,ofm_num-tm);
	if(ltype == LT_CONV)
	{
		for(tn = 0; tn < ifm_num; tn += TN)
		{
			TN_MIN = MIN_diy(TN, ifm_num-tn);
			weight_load_tile_buf_reorder(weight, weight_reorder, INumxKK, KK, TM_MIN, TN_MIN, tm, tn, ksize, Tn);
		}
	}else if(ltype == LT_DCONV)
	{	
		weight_load_tile_buf_reorder(weight, weight_reorder,      KK, KK, TM_MIN,      1, tm,  0, ksize,  1);
	}
    }	

}

int Reorder_weight_64(float *Weight, int IFM_NUM,int OFM_NUM,int Ksize,int TM,int TN, int32_t *add_offset, FILE *fp)
{
	const int KxK = Ksize*Ksize;
	const int IFM_NUMxKxK = IFM_NUM*KxK;

	int m,n;
	int tm,tn,tk;

	float weight_buffer[Tm*Tn*K*K+64];
	
	int offset_cur = IFM_NUM*OFM_NUM*KxK;
	printf("Layer[x]; param num=%12d ", offset_cur);
	*add_offset = 0;

	int TM_MIN,TN_MIN;
	for( m = 0; m < OFM_NUM; m += TM)
	{
		TM_MIN = MIN_diy(TM,OFM_NUM - m);
		for(n = 0;n < IFM_NUM; n += TN)
		{
			TN_MIN = MIN_diy(TN,IFM_NUM - n);
			int Woffset = m*IFM_NUMxKxK + n*KxK;
			for(tm = 0;tm < TM_MIN; tm++)
			{
				memcpy((float *)(weight_buffer + tm*TN_MIN*KxK),
					(float *)(Weight + tm*IFM_NUMxKxK + Woffset),TN_MIN*KxK*sizeof(float));
			}
			
			int local_offset = TM_MIN*TN_MIN*KxK;
			if(local_offset & 0x1)
			{
				*add_offset = *add_offset + 1;
			}
			fwrite(weight_buffer,sizeof(float), local_offset + (local_offset & 0x1), fp);
			//memcpy((float *)(Weight_reorg+offset),weight_buffer2,TM_MIN*TN_MIN*KxK*sizeof(float));
		}							
	}
	printf("add_offset = %d\n", *add_offset);

	return 0;
}

#define MIN_VALUE (-1024*1024*1024)
#define MAX_VALUE (1024*1024*1024)

void weight_reorder_tile_f32(float *weight, int INumxKK, int KK, int TM_MIN, int TN_MIN, int tm, int tn, int ksize, 
								int OFM_NUM, int *add_offset, FILE *fp)//assume w equal h
{
	static float local_wbuf[(Tm+3)/4*Tn*K*K*4];
	int t1,t2,t3,t4,t5;

	int offset = 0;
	float tmp_in_float;
	int ma3_d3 = (TM_MIN+3)/4;
	for(t1 = 0;t1 < ma3_d3; t1++)
	for(t2 = 0;t2 < TN_MIN; t2++)
	for(t3 = 0;t3 <ksize; t3++)
	for(t4 = 0;t4 <ksize; t4++)
	for(t5 = 0;t5 <4; t5++)
	{
		int idx_i = (t1*4+t5+tm)*INumxKK + (t2+tn)*KK + t3*ksize + t4;
		int idx_o = t1*TN_MIN*KK*4 + t2*KK*4 + t3*ksize*4 + t4*4 + t5;
		if((t1*4+t5)>=TM_MIN)
			tmp_in_float = 0;
		else
			tmp_in_float = weight[idx_i];
		local_wbuf[idx_o] = tmp_in_float;
	}

	uint16_t val_left = ma3_d3*TN_MIN*KK*4 - TM_MIN*TN_MIN*KK;
	*add_offset = *add_offset + val_left;
	fwrite(local_wbuf,sizeof(float), ma3_d3*TN_MIN*KK*4, fp);

	// fwrite(local_wbuf,sizeof(short), offset, fp);

}

int reorg_weight_ft32(float *Weight, int IFM_NUM,int OFM_NUM,int Ksize,int TM,int TN, int ltype, int *add_offset, FILE *fp)
{
	const int KxK = Ksize*Ksize;
	const int IFM_NUMxKxK = IFM_NUM*KxK;

	int tm,tn,tk;

	int offset_cur = IFM_NUM*OFM_NUM*KxK;
	printf("Layer[x]; param num=%12d ", offset_cur);
	*add_offset = 0;

	double max_error,min_error,sum_error;
	sum_error = 0;
	max_error = MIN_VALUE;
	min_error = MAX_VALUE;

	int TM_MIN,TN_MIN;
    for(tm = 0; tm < OFM_NUM; tm += TM)
    {
		TM_MIN = MIN_diy(TM,OFM_NUM-tm);
		if(ltype == LT_CONV)
		{
			for(tn = 0; tn < IFM_NUM; tn += TN)
			{
				TN_MIN = MIN_diy(TN, IFM_NUM-tn);
				weight_reorder_tile_f32(Weight, IFM_NUMxKxK, KxK, TM_MIN, TN_MIN, tm, tn, Ksize, OFM_NUM, add_offset, fp);
			}
		}else if(ltype == LT_DCONV)
		{	
			weight_reorder_tile_f32(Weight,      KxK, KxK, TM_MIN,      1, tm,  0, Ksize, OFM_NUM, add_offset, fp);
		}
    }	

	printf("add_offset = %d\n", *add_offset);

	return 0;
}

void weight_reorder_tile_i16_c4(float *weight, int INumxKK, int KK, int TM_MIN, int TN_MIN, int tm, int tn, int ksize, 
		int maxQ, int *add_offset, FILE *fp, double *max_error, double *min_error, double *sum_error)//assume w equal h
{
	static short local_wbuf[(Tm+3)/4*Tn*K*K*4];
	int t1,t2,t3,t4, t5;

	double max_error_l, min_error_l, sum_error_l;
	max_error_l = *max_error;
	min_error_l = *min_error;
	sum_error_l = *sum_error;

	// int offset = 0;
	int ma3_d3 = (TM_MIN+3)/4;
	float tmp_in_float;
	for(t1 = 0;t1 < ma3_d3; t1++)
	for(t2 = 0;t2 < TN_MIN; t2++)
	for(t3 = 0;t3 <ksize; t3++)
	for(t4 = 0;t4 <ksize; t4++)
	for(t5 = 0;t5 <4; t5++)
	{
		int idx_i = (t1*4+t5+tm)*INumxKK + (t2+tn)*KK + t3*ksize + t4;
		int idx_o = t1*TN_MIN*KK*4 + t2*KK*4 + t3*ksize*4 + t4*4 + t5;
		if((t1*4+t5)>=TM_MIN)
			tmp_in_float = 0;
		else
			tmp_in_float = weight[idx_i];

		short tmp_fixed = (short)(tmp_in_float*pow(2.0,maxQ));
		float tmp_out_float = (float)tmp_fixed*pow(2.0,-maxQ);
		double error = (tmp_out_float - tmp_in_float)*(tmp_out_float - tmp_in_float);
		error = sqrt(error);
		sum_error_l += error;
		if(error<min_error_l)
			min_error_l = error;
		if(error>max_error_l)
			max_error_l = error;

		// local_wbuf[offset] = tmp_fixed;
		local_wbuf[idx_o] = tmp_fixed;
		// offset++;
	}

	*max_error = max_error_l;
	*min_error = min_error_l;
	*sum_error = sum_error_l;

	uint16_t val_left = ma3_d3*TN_MIN*KK*4 - TM_MIN*TN_MIN*KK;
	*add_offset = *add_offset + val_left;
	fwrite(local_wbuf,sizeof(short), ma3_d3*TN_MIN*KK*4, fp);	

	// uint16_t val_left = 0;
	// if(offset & 0x3)
	// {
	// 	val_left = 4 - (offset & 0x3);
	// 	*add_offset = *add_offset + val_left;
	// 	// local_wbuf[offset] = 0;
	// 	for(uint32_t j=0; j<val_left;j++){
	// 		local_wbuf[offset+j]=0;
	// 	}		
	// }
	// fwrite(local_wbuf,sizeof(short), offset + val_left, fp);

	// fwrite(local_wbuf,sizeof(short), offset, fp);

}

int reorg_quantize_weight_int16_c4(float *Weight, int IFM_NUM,int OFM_NUM,int Ksize,int TM,int TN, int ltype, float *ap16_range,int *maxQ_array, int *add_offset, FILE *fp)
{
	const int KxK = Ksize*Ksize;
	const int IFM_NUMxKxK = IFM_NUM*KxK;

	int tm,tn,tk;

	int offset_cur = IFM_NUM*OFM_NUM*KxK;
	printf("Layer[x]; param num=%12d ", offset_cur);
	*add_offset = 0;

	float min,max;
	min = MAX_VALUE;
	max = MIN_VALUE;
	for(int j=0;j<offset_cur;j++)
	{
		float tmp_in_float = Weight[j];
		if(tmp_in_float<min)
			min = tmp_in_float;
		if(tmp_in_float>max)
			max = tmp_in_float;
	}
	printf("float min=%.7lf,max=%.7lf ",min,max);//find float min max

	int maxQ = -1;
	for(int k=0;k<16;k++)//find maxQ
	{
		if(min>ap16_range[2*k]&&max<ap16_range[2*k+1])
		{
			maxQ = k;
		}
		else if(k==0)
		{
			printf("beyond Q0 min=%.7lf,max=%.7lf ",min,max);
			break;
		}
	}
	printf("maxQ=%d ",maxQ);
	*maxQ_array = maxQ;

	double max_error,min_error,sum_error;
	sum_error = 0;
	max_error = MIN_VALUE;
	min_error = MAX_VALUE;

	int TM_MIN,TN_MIN;
    for(tm = 0; tm < OFM_NUM; tm += TM)
    {
		TM_MIN = MIN_diy(TM,OFM_NUM-tm);
		if(ltype == LT_CONV)
		{
			for(tn = 0; tn < IFM_NUM; tn += TN)
			{
				TN_MIN = MIN_diy(TN, IFM_NUM-tn);
				weight_reorder_tile_i16_c4(Weight, IFM_NUMxKxK, KxK, TM_MIN, TN_MIN, tm, tn, Ksize, maxQ, add_offset, fp,&max_error,&min_error,&sum_error);
			}
		}else if(ltype == LT_DCONV)
		{	
			weight_reorder_tile_i16_c4(Weight,      KxK, KxK, TM_MIN,      1, tm,  0, Ksize, maxQ, add_offset, fp,&max_error,&min_error,&sum_error);
		}
    }	

	printf("sum2_error = %.7lf,min_error=%.7lf,max_error=%.7lf,add_offset = %d",sum_error,min_error,max_error, *add_offset);
	printf("\n");

	return 0;
}



void weight_reorder_tile_i16(float *weight, int INumxKK, int KK, int TM_MIN, int TN_MIN, int tm, int tn, int ksize, 
		int maxQ, int *add_offset, FILE *fp, double *max_error, double *min_error, double *sum_error)//assume w equal h
{
	static short local_wbuf[Tm*Tn*K*K+8];
	int t1,t2,t3,t4;

	double max_error_l, min_error_l, sum_error_l;
	max_error_l = *max_error;
	min_error_l = *min_error;
	sum_error_l = *sum_error;

	int offset = 0;
	for(t1 = 0;t1 < TM_MIN; t1++)
	for(t2 = 0;t2 < TN_MIN; t2++)
	for(t3 = 0;t3 <ksize; t3++)
	for(t4 = 0;t4 <ksize; t4++)
	{
		float tmp_in_float = weight[(t1+tm)*INumxKK + (t2+tn)*KK + t3*ksize + t4];
		short tmp_fixed = (short)(tmp_in_float*pow(2.0,maxQ));
		float tmp_out_float = (float)tmp_fixed*pow(2.0,-maxQ);
		double error = (tmp_out_float - tmp_in_float)*(tmp_out_float - tmp_in_float);
		error = sqrt(error);
		sum_error_l += error;
		if(error<min_error_l)
			min_error_l = error;
		if(error>max_error_l)
			max_error_l = error;

		local_wbuf[offset] = tmp_fixed;
		offset++;
	}

	*max_error = max_error_l;
	*min_error = min_error_l;
	*sum_error = sum_error_l;	

	uint16_t val_left = 0;
	if(offset & 0x3)
	{
		val_left = 4 - (offset & 0x3);
		*add_offset = *add_offset + val_left;
		// local_wbuf[offset] = 0;
		for(uint32_t j=0; j<val_left;j++){
			local_wbuf[offset+j]=0;
		}		
	}
	fwrite(local_wbuf,sizeof(short), offset + val_left, fp);

	// fwrite(local_wbuf,sizeof(short), offset, fp);

}

int reorg_quantize_weight_int16(float *Weight, int IFM_NUM,int OFM_NUM,int Ksize,int TM,int TN, int ltype, float *ap16_range,int *maxQ_array, int *add_offset, FILE *fp)
{
	const int KxK = Ksize*Ksize;
	const int IFM_NUMxKxK = IFM_NUM*KxK;

	int tm,tn,tk;

	int offset_cur = IFM_NUM*OFM_NUM*KxK;
	printf("Layer[x]; param num=%12d ", offset_cur);
	*add_offset = 0;

	float min,max;
	min = MAX_VALUE;
	max = MIN_VALUE;
	for(int j=0;j<offset_cur;j++)
	{
		float tmp_in_float = Weight[j];
		if(tmp_in_float<min)
			min = tmp_in_float;
		if(tmp_in_float>max)
			max = tmp_in_float;
	}
	printf("float min=%.7lf,max=%.7lf ",min,max);//find float min max

	int maxQ = -1;
	for(int k=0;k<16;k++)//find maxQ
	{
		if(min>ap16_range[2*k]&&max<ap16_range[2*k+1])
		{
			maxQ = k;
		}
		else if(k==0)
		{
			printf("beyond Q0 min=%.7lf,max=%.7lf ",min,max);
			break;
		}
	}
	printf("maxQ=%d ",maxQ);
	*maxQ_array = maxQ;

	double max_error,min_error,sum_error;
	sum_error = 0;
	max_error = MIN_VALUE;
	min_error = MAX_VALUE;

	int TM_MIN,TN_MIN;
    for(tm = 0; tm < OFM_NUM; tm += TM)
    {
		TM_MIN = MIN_diy(TM,OFM_NUM-tm);
		if(ltype == LT_CONV)
		{
			for(tn = 0; tn < IFM_NUM; tn += TN)
			{
				TN_MIN = MIN_diy(TN, IFM_NUM-tn);
				weight_reorder_tile_i16(Weight, IFM_NUMxKxK, KxK, TM_MIN, TN_MIN, tm, tn, Ksize, maxQ, add_offset, fp,&max_error,&min_error,&sum_error);
			}
		}else if(ltype == LT_DCONV)
		{	
			weight_reorder_tile_i16(Weight,      KxK, KxK, TM_MIN,      1, tm,  0, Ksize, maxQ, add_offset, fp,&max_error,&min_error,&sum_error);
		}
    }	

	printf("sum2_error = %.7lf,min_error=%.7lf,max_error=%.7lf,add_offset = %d",sum_error,min_error,max_error, *add_offset);
	printf("\n");

	return 0;
}

// int reorg_quantize_weight_int16(float *Weight, int IFM_NUM,int OFM_NUM,int Ksize,int TM,int TN,float *ap16_range,int *maxQ_array, int *add_offset, FILE *fp)
// {
// 	const int KxK = Ksize*Ksize;
// 	const int IFM_NUMxKxK = IFM_NUM*KxK;

// 	int m,n;
// 	int tm,tn,tk;

// 	float weight_buffer[Tm*Tn*K*K];
// 	short weight_buffer2[Tm*Tn*K*K+32];
	
// 	int offset_cur = IFM_NUM*OFM_NUM*KxK;
// 	printf("Layer[x]; param num=%12d ", offset_cur);
// 	*add_offset = 0;

// 	int j;
// 	float min,max;
// 	min = MAX_VALUE;
// 	max = MIN_VALUE;
// 	for(j=0;j<offset_cur;j++)
// 	{
// 		float tmp_in_float = Weight[j];
// 		if(tmp_in_float<min)
// 			min = tmp_in_float;
// 		if(tmp_in_float>max)
// 			max = tmp_in_float;
// 	}
// 	printf("float min=%.7lf,max=%.7lf ",min,max);//find float min max
// 	int k;
// 	int maxQ = -1;
// 	for(k=0;k<16;k++)//find maxQ
// 	{
// 		if(min>ap16_range[2*k]&&max<ap16_range[2*k+1])
// 		{
// 			maxQ = k;
// 		}
// 		else if(k==0)
// 		{
// 			printf("beyond Q0 min=%.7lf,max=%.7lf ",min,max);
// 			break;
// 		}
// 	}
// 	printf("maxQ=%d ",maxQ);
// 	*maxQ_array = maxQ;

// 	double max_error,min_error,sum_error;
// 	sum_error = 0;
// 	max_error = MIN_VALUE;
// 	min_error = MAX_VALUE;

// 	int TM_MIN,TN_MIN;
// 	for( m = 0; m < OFM_NUM; m += TM)
// 	{
// 		TM_MIN = MIN_diy(TM,OFM_NUM - m);
// 		for(n = 0;n < IFM_NUM; n += TN)
// 		{
// 			TN_MIN = MIN_diy(TN,IFM_NUM - n);
// 			int Woffset = m*IFM_NUMxKxK + n*KxK;
// 			for(tm = 0;tm < TM_MIN; tm++)
// 			{
// 				memcpy((float *)(weight_buffer + tm*TN_MIN*KxK),
// 					(float *)(Weight + tm*IFM_NUMxKxK + Woffset),TN_MIN*KxK*sizeof(float));
// 			}

// 			int TN_MINxTM_MIN = TN_MIN*TM_MIN;
// 			for(tk = 0;tk < KxK; tk++)
// 				for(tm = 0;tm < TM_MIN; tm++)
// 					for(tn = 0;tn < TN_MIN;tn++)
// 					{
// 						//weight_buffer2[tk*TN_MINxTM_MIN + tm*TN_MIN + tn] = weight_buffer[tm*TN_MIN*KxK + tn*KxK + tk];
// 						float tmp_in_float = weight_buffer[tm*TN_MIN*KxK + tn*KxK + tk];
// 						short tmp_fixed = (short)(tmp_in_float*pow(2.0,maxQ));
// 						float tmp_out_float = (float)tmp_fixed*pow(2.0,-maxQ);
// 						double error = (tmp_out_float - tmp_in_float)*(tmp_out_float - tmp_in_float);
// 						error = sqrt(error);
// 						sum_error += error;
// 						if(error<min_error)
// 							min_error = error;
// 						if(error>max_error)
// 							max_error = error;

// 						weight_buffer2[tk*TN_MINxTM_MIN + tm*TN_MIN + tn] = tmp_fixed;
// 					}
			
// 			int local_offset = TM_MIN*TN_MIN*KxK;
// 			if(local_offset & 0x1)
// 			{
// 				*add_offset = *add_offset + 1;
// 			}
// 			fwrite(weight_buffer2,sizeof(short), local_offset + (local_offset & 0x1), fp);
// 			//memcpy((float *)(Weight_reorg+offset),weight_buffer2,TM_MIN*TN_MIN*KxK*sizeof(float));
// 		}							
// 	}

// 	printf("sum2_error = %.7lf,min_error=%.7lf,max_error=%.7lf,add_offset = %d",sum_error,min_error,max_error, *add_offset);
// 	printf("\n");

// 	return 0;
// }

void weight_reorder_tile_i16c_scale(float *weight, int16_t *local_wbuf, int INumxKK, int KK, int TM_MIN, int TN_MIN, int tm, int tn, int ksize, int lane_num, 
		 int maxQ, int *add_offset, FILE *fp, double *max_error, double *min_error, double *sum_error, float *scale_array, int scale_maxQ)
{
	int t1,t2,t3,t4, t5;

	double max_error_l, min_error_l, sum_error_l;
	max_error_l = *max_error;
	min_error_l = *min_error;
	sum_error_l = *sum_error;

	int ma3_d3 = (TM_MIN+lane_num-1)/lane_num;
	float tmp_in_float;

	float tmp_exp_w = pow(2.0,maxQ);
	float tmp_exp_w_d = pow(2.0,-maxQ);
	float tmp_exp_scale = pow(2.0,scale_maxQ);
	float tmp_exp_scale_d = pow(2.0,-scale_maxQ);

	for(t1 = 0;t1 < ma3_d3; t1++)
	for(t2 = 0;t2 < TN_MIN; t2++)
	for(t3 = 0;t3 <ksize; t3++)
	for(t4 = 0;t4 <ksize; t4++)
	for(t5 = 0;t5 <lane_num; t5++)
	{
		int m = (t1*lane_num+t5+tm);
		int idx_i = (t1*lane_num+t5+tm)*INumxKK + (t2+tn)*KK + t3*ksize + t4;
		int idx_o = t1*TN_MIN*KK*lane_num + t2*KK*lane_num + t3*ksize*lane_num + t4*lane_num + t5;

		float scale_f32;
		if((t1*lane_num+t5)>=TM_MIN){
			scale_f32 = 0;
			tmp_in_float = 0;
		}
		else{
			scale_f32 = scale_array[m];
			float w_in = weight[idx_i];
			float w_in_abs = abs(w_in); 
			if(w_in_abs>=scale_f32){
				int16_t scale_i16 = scale_f32*tmp_exp_scale;
				scale_i16 = scale_i16 -1;
				float scale_f32_1 = scale_i16*tmp_exp_scale_d;
				if(w_in<0){
					w_in = -scale_f32_1;
				}
				else{
					w_in = scale_f32_1;
				}
			}
			tmp_in_float = w_in/scale_f32;
		}

		int16_t tmp_fixed = (int16_t)(tmp_in_float*tmp_exp_w);
		float tmp_out_float = (float)tmp_fixed*tmp_exp_w_d;
		double error = (tmp_out_float - tmp_in_float)*(tmp_out_float - tmp_in_float)*scale_f32*scale_f32;
		error = sqrt(error);
		sum_error_l += error;
		if(error<min_error_l)
			min_error_l = error;
		if(error>max_error_l)
			max_error_l = error;

		local_wbuf[idx_o] = tmp_fixed;
	}

	uint16_t val_left = ma3_d3*TN_MIN*KK*lane_num - TM_MIN*TN_MIN*KK;
	*add_offset = *add_offset + val_left;
	fwrite(local_wbuf,sizeof(int16_t), ma3_d3*TN_MIN*KK*lane_num, fp);

	*max_error = max_error_l;
	*min_error = min_error_l;
	*sum_error = sum_error_l;
}

int weight_reorder_i16c_scale(float *Weight, int IFM_NUM,int OFM_NUM,int Ksize,int TM,int TN, int ltype, float *ap16_range,int *scaleQ_array, 
							int lane_num, int *add_offset, float *scale_array, FILE *fp)
{
	const int KxK = Ksize*Ksize;
	const int IFM_NUMxKxK = IFM_NUM*KxK;

	// float local_wbuf[(Tm+3)/4*Tn*K*K*4];
	int lbuf_size = ((Tm+lane_num-1)/lane_num)*Tn*K*K*lane_num;
	int16_t *local_wbuf = (int16_t *)malloc(lbuf_size*sizeof(int16_t));

	int tm,tn,tk;
	// int offset_cur = IFM_NUM*OFM_NUM*KxK;
	// printf("Layer[x]; param num=%12d ", offset_cur);
	*add_offset = 0;
	// int maxQ = 14;
	int maxQ = 15;
	printf("maxQ=%d ",maxQ);
	// *maxQ_array = maxQ;

	float scale_max = 0;
	float scale_min = 1024*1024;
	for(tm = 0; tm < OFM_NUM; tm++)//calc sacle param
    {
		float max = 0;
		int m_offset = tm*IFM_NUMxKxK;
		for(int tnkk = 0; tnkk < IFM_NUMxKxK; tnkk++){
			float tmp_f32 = Weight[m_offset + tnkk];
			tmp_f32 = abs(tmp_f32);
			if(tmp_f32 > max)
				max = tmp_f32;
		}
		scale_array[tm] = max;
		if(max > scale_max)
			scale_max = max;
		if(max < scale_min)
			scale_min = max;
    }
	printf("float scale abs min=%.7lf,max=%.7lf ", scale_min, scale_max);//find float min max

	int scale_maxQ = -1;
	for(int k=0;k<16;k++)//find maxQ
	{
		if(scale_min>ap16_range[2*k]&&scale_max<ap16_range[2*k+1])//here make sure that all scale_maxQ belong to (min, max)
		{
			scale_maxQ = k;
		}
		else if(k==0)
		{
			printf("beyond Q0 min=%.7lf,max=%.7lf ",scale_min,scale_max);
			break;
		}
	}
	printf("scale_maxQ=%d ,", scale_maxQ);
	*scaleQ_array = scale_maxQ;

	float tmp_exp = pow(2.0, scale_maxQ);
	float tmp_exp_d = pow(2.0, -scale_maxQ);
	for(tm = 0; tm < OFM_NUM; tm++)//calc sacle param
    {
		float scale_l = scale_array[tm];
		int16_t scale_l_16b = scale_l*tmp_exp;
		scale_array[tm] = (scale_l_16b+1)*tmp_exp_d;//to make all params in single channel (-1, 1), +1 will not beyond the bound [min, max]
    }				

	double max_error,min_error,sum_error;
	sum_error = 0;
	max_error = MIN_VALUE;
	min_error = MAX_VALUE;	

    for(tm = 0; tm < OFM_NUM; tm += TM)
    {
		int TM_MIN = MIN_diy(TM,OFM_NUM-tm);
		if(ltype == LT_CONV)
		{
			for(tn = 0; tn < IFM_NUM; tn += TN)
			{
				int TN_MIN = MIN_diy(TN, IFM_NUM-tn);
				weight_reorder_tile_i16c_scale(Weight, local_wbuf, IFM_NUMxKxK, KxK, TM_MIN, TN_MIN, tm, tn, Ksize, lane_num, 
											maxQ, add_offset, fp, &max_error,&min_error,&sum_error, scale_array, scale_maxQ);
			}
		}else if(ltype == LT_DCONV)
		{	
			weight_reorder_tile_i16c_scale(Weight, local_wbuf, KxK, KxK, TM_MIN,      1, tm,  0, Ksize, lane_num,
										 maxQ, add_offset, fp, &max_error,&min_error,&sum_error, scale_array, scale_maxQ);
		}
    }	

	printf("sum2_error = %.7lf,min_error=%.7lf,max_error=%.7lf,add_offset = %d",sum_error,min_error,max_error, *add_offset);
	printf("\n");

	free(local_wbuf);

	return 0;
}

// void weight_reorder_tile_i16c_scale(float *weight, int16_t *local_wbuf, int INumxKK, int KK, int TM_MIN, int TN_MIN, int tm, int tn, int ksize, int lane_num, 
// 		 int maxQ, int *add_offset, FILE *fp, double *max_error, double *min_error, double *sum_error, float *scale_array)
// {
// 	int t1,t2,t3,t4, t5;

// 	double max_error_l, min_error_l, sum_error_l;
// 	max_error_l = *max_error;
// 	min_error_l = *min_error;
// 	sum_error_l = *sum_error;

// 	int ma3_d3 = (TM_MIN+lane_num-1)/lane_num;
// 	float tmp_in_float;
// 	for(t1 = 0;t1 < ma3_d3; t1++)
// 	for(t2 = 0;t2 < TN_MIN; t2++)
// 	for(t3 = 0;t3 <ksize; t3++)
// 	for(t4 = 0;t4 <ksize; t4++)
// 	for(t5 = 0;t5 <lane_num; t5++)
// 	{
// 		int m = (t1*lane_num+t5+tm);
// 		int idx_i = (t1*lane_num+t5+tm)*INumxKK + (t2+tn)*KK + t3*ksize + t4;
// 		int idx_o = t1*TN_MIN*KK*lane_num + t2*KK*lane_num + t3*ksize*lane_num + t4*lane_num + t5;

// 		float scale_f32;
// 		if((t1*lane_num+t5)>=TM_MIN){
// 			scale_f32 = 0;
// 			tmp_in_float = 0;
// 		}
// 		else{
// 			scale_f32 = scale_array[m];
// 			tmp_in_float = weight[idx_i]/scale_array[m];
// 		}

// 		int16_t tmp_fixed = (int16_t)(tmp_in_float*pow(2.0,maxQ));
// 		float tmp_out_float = (float)tmp_fixed*pow(2.0,-maxQ);
// 		double error = (tmp_out_float - tmp_in_float)*(tmp_out_float - tmp_in_float)*scale_f32*scale_f32;
// 		error = sqrt(error);
// 		sum_error_l += error;
// 		if(error<min_error_l)
// 			min_error_l = error;
// 		if(error>max_error_l)
// 			max_error_l = error;

// 		local_wbuf[idx_o] = tmp_fixed;
// 	}

// 	uint16_t val_left = ma3_d3*TN_MIN*KK*lane_num - TM_MIN*TN_MIN*KK;
// 	*add_offset = *add_offset + val_left;
// 	fwrite(local_wbuf,sizeof(int16_t), ma3_d3*TN_MIN*KK*lane_num, fp);

// 	*max_error = max_error_l;
// 	*min_error = min_error_l;
// 	*sum_error = sum_error_l;
// }

// int weight_reorder_i16c_scale(float *Weight, int IFM_NUM,int OFM_NUM,int Ksize,int TM,int TN, int ltype, float *ap16_range,int *maxQ_array, 
// 							int lane_num, int *add_offset, float *scale_array, FILE *fp)
// {
// 	const int KxK = Ksize*Ksize;
// 	const int IFM_NUMxKxK = IFM_NUM*KxK;

// 	// float local_wbuf[(Tm+3)/4*Tn*K*K*4];
// 	int lbuf_size = ((Tm+lane_num-1)/lane_num)*Tn*K*K*lane_num;
// 	int16_t *local_wbuf = (int16_t *)malloc(lbuf_size*sizeof(int16_t));

// 	int tm,tn,tk;
// 	// int offset_cur = IFM_NUM*OFM_NUM*KxK;
// 	// printf("Layer[x]; param num=%12d ", offset_cur);
// 	*add_offset = 0;
// 	int maxQ = 14;
// 	printf("maxQ=%d ",maxQ);
// 	*maxQ_array = maxQ;

// 	for(tm = 0; tm < OFM_NUM; tm++)//calc sacle param
//     {
// 		float max = 0;
// 		int m_offset = tm*IFM_NUMxKxK;
// 		for(int tnkk = 0; tnkk < IFM_NUMxKxK; tnkk++){
// 			float tmp_f32 = Weight[m_offset + tnkk];
// 			tmp_f32 = abs(tmp_f32);
// 			if(tmp_f32 > max)
// 				max = tmp_f32;
// 		}
// 		scale_array[tm] = max;
//     }		

// 	double max_error,min_error,sum_error;
// 	sum_error = 0;
// 	max_error = MIN_VALUE;
// 	min_error = MAX_VALUE;	

//     for(tm = 0; tm < OFM_NUM; tm += TM)
//     {
// 		int TM_MIN = MIN_diy(TM,OFM_NUM-tm);
// 		if(ltype == LT_CONV)
// 		{
// 			for(tn = 0; tn < IFM_NUM; tn += TN)
// 			{
// 				int TN_MIN = MIN_diy(TN, IFM_NUM-tn);
// 				weight_reorder_tile_i16c_scale(Weight, local_wbuf, IFM_NUMxKxK, KxK, TM_MIN, TN_MIN, tm, tn, Ksize, lane_num, 
// 											maxQ, add_offset, fp, &max_error,&min_error,&sum_error, scale_array);
// 			}
// 		}else if(ltype == LT_DCONV)
// 		{	
// 			weight_reorder_tile_i16c_scale(Weight, local_wbuf, KxK, KxK, TM_MIN,      1, tm,  0, Ksize, lane_num,
// 										 maxQ, add_offset, fp, &max_error,&min_error,&sum_error, scale_array);
// 		}
//     }	

// 	printf("sum2_error = %.7lf,min_error=%.7lf,max_error=%.7lf,add_offset = %d",sum_error,min_error,max_error, *add_offset);
// 	printf("\n");

// 	free(local_wbuf);

// 	return 0;
// }

void weight_reorder_tile_i16c(float *weight, int16_t *local_wbuf, int INumxKK, int KK, int TM_MIN, int TN_MIN, int tm, int tn, int ksize, int lane_num, 
		 int maxQ, int *add_offset, FILE *fp, double *max_error, double *min_error, double *sum_error)
{
	int t1,t2,t3,t4, t5;

	double max_error_l, min_error_l, sum_error_l;
	max_error_l = *max_error;
	min_error_l = *min_error;
	sum_error_l = *sum_error;

	int ma3_d3 = (TM_MIN+lane_num-1)/lane_num;
	float tmp_in_float;
	for(t1 = 0;t1 < ma3_d3; t1++)
	for(t2 = 0;t2 < TN_MIN; t2++)
	for(t3 = 0;t3 <ksize; t3++)
	for(t4 = 0;t4 <ksize; t4++)
	for(t5 = 0;t5 <lane_num; t5++)
	{
		int idx_i = (t1*lane_num+t5+tm)*INumxKK + (t2+tn)*KK + t3*ksize + t4;
		int idx_o = t1*TN_MIN*KK*lane_num + t2*KK*lane_num + t3*ksize*lane_num + t4*lane_num + t5;
		if((t1*lane_num+t5)>=TM_MIN)
			tmp_in_float = 0;
		else
			tmp_in_float = weight[idx_i];

		int16_t tmp_fixed = (int16_t)(tmp_in_float*pow(2.0,maxQ));
		float tmp_out_float = (float)tmp_fixed*pow(2.0,-maxQ);
		double error = (tmp_out_float - tmp_in_float)*(tmp_out_float - tmp_in_float);
		error = sqrt(error);
		sum_error_l += error;
		if(error<min_error_l)
			min_error_l = error;
		if(error>max_error_l)
			max_error_l = error;

		local_wbuf[idx_o] = tmp_fixed;
	}

	uint16_t val_left = ma3_d3*TN_MIN*KK*lane_num - TM_MIN*TN_MIN*KK;
	*add_offset = *add_offset + val_left;
	fwrite(local_wbuf,sizeof(int16_t), ma3_d3*TN_MIN*KK*lane_num, fp);

	*max_error = max_error_l;
	*min_error = min_error_l;
	*sum_error = sum_error_l;
}

int weight_reorder_i16c(float *Weight, int IFM_NUM,int OFM_NUM,int Ksize,int TM,int TN, int ltype, float *ap16_range,int *maxQ_array, 
							int lane_num, int *add_offset, FILE *fp)
{
	const int KxK = Ksize*Ksize;
	const int IFM_NUMxKxK = IFM_NUM*KxK;

	// float local_wbuf[(Tm+3)/4*Tn*K*K*4];
	int lbuf_size = ((Tm+lane_num-1)/lane_num)*Tn*K*K*lane_num;
	int16_t *local_wbuf = (int16_t *)malloc(lbuf_size*sizeof(int16_t));

	int tm,tn,tk;

	int offset_cur = IFM_NUM*OFM_NUM*KxK;
	printf("Layer[x]; param num=%12d ", offset_cur);
	*add_offset = 0;

	float min,max;
	min = MAX_VALUE;
	max = MIN_VALUE;
	for(int j=0;j<offset_cur;j++)
	{
		float tmp_in_float = Weight[j];
		if(tmp_in_float<min)
			min = tmp_in_float;
		if(tmp_in_float>max)
			max = tmp_in_float;
	}
	printf("float min=%.7lf,max=%.7lf ",min,max);//find float min max

	int maxQ = -1;
	for(int k=0;k<16;k++)//find maxQ
	{
		if(min>ap16_range[2*k]&&max<ap16_range[2*k+1])
		{
			maxQ = k;
		}
		else if(k==0)
		{
			printf("beyond Q0 min=%.7lf,max=%.7lf ",min,max);
			break;
		}
	}
	printf("maxQ=%d ",maxQ);
	*maxQ_array = maxQ;

	double max_error,min_error,sum_error;
	sum_error = 0;
	max_error = MIN_VALUE;
	min_error = MAX_VALUE;	

	int TM_MIN,TN_MIN;
    for(tm = 0; tm < OFM_NUM; tm += TM)
    {
		TM_MIN = MIN_diy(TM,OFM_NUM-tm);
		if(ltype == LT_CONV)
		{
			for(tn = 0; tn < IFM_NUM; tn += TN)
			{
				TN_MIN = MIN_diy(TN, IFM_NUM-tn);
				weight_reorder_tile_i16c(Weight, local_wbuf, IFM_NUMxKxK, KxK, TM_MIN, TN_MIN, tm, tn, Ksize, lane_num, 
											maxQ, add_offset, fp, &max_error,&min_error,&sum_error);
			}
		}else if(ltype == LT_DCONV)
		{	
			weight_reorder_tile_i16c(Weight, local_wbuf, KxK, KxK, TM_MIN,      1, tm,  0, Ksize, lane_num,
										 maxQ, add_offset, fp, &max_error,&min_error,&sum_error);
		}
    }	

	printf("sum2_error = %.7lf,min_error=%.7lf,max_error=%.7lf,add_offset = %d",sum_error,min_error,max_error, *add_offset);
	printf("\n");

	free(local_wbuf);

	return 0;
}

int bias_reorder_i16c(float *in,int *offset,int layer_num, int *ltype_set, float *ap16_range,int *maxQ_array, int lane_num, int *add_offset, FILE *fp)
{
	static int16_t bias_buf_c[MAX_BETA_LENGTH+128];

	int woffset = 0;
	int offset_cur = 0;
	for(int i=0;i<layer_num;i++)
	{
		if((ltype_set[i]==LT_AVGPOOL)||(ltype_set[i] == LT_MAXPOOL)){
			continue;
		}
		offset_cur = offset[i];
		printf("Layer %2d;param num=%12d ",i, offset_cur);

		float min,max;
		min = MAX_VALUE;
		max = MIN_VALUE;
		for(int j=0;j<offset_cur;j++)
		{
			float tmp_in_float = in[woffset+j];
			if(tmp_in_float<min)
				min = tmp_in_float;
			if(tmp_in_float>max)
				max = tmp_in_float;
		}
		printf("float min=%.7lf,max=%.7lf ",min,max);//find float min max

		int k;
		int maxQ = -1;
		for(k=0;k<16;k++)//find maxQ
		{
			if(min>ap16_range[2*k]&&max<ap16_range[2*k+1])
			{
				maxQ = k;
			}
			else if(k==0)
			{
				printf("beyond Q0 min=%.7lf,max=%.7lf ",min,max);
				break;
			}
		}
		printf("maxQ=%d ",maxQ);
		maxQ_array[i] = maxQ;		

		double max_error,min_error,sum_error;
		sum_error = 0;
		max_error = MIN_VALUE;
		min_error = MAX_VALUE;

		int data_size = ((offset_cur + lane_num -1)/lane_num)*lane_num;
		for(int j=0;j<data_size;j++)
		{
			float tmp_in_float;
			if(j < offset_cur)
				tmp_in_float = in[woffset+j];
			else
				tmp_in_float = 0;

			int16_t tmp_fixed = (int16_t)(tmp_in_float*pow(2.0,maxQ));
			float tmp_out_float = (float)tmp_fixed*pow(2.0,-maxQ);
			double error = (tmp_out_float - tmp_in_float)*(tmp_out_float - tmp_in_float);
			error = sqrt(error);
			sum_error += error;
			if(error<min_error)
				min_error = error;
			if(error>max_error)
				max_error = error;

			bias_buf_c[j] = tmp_fixed;
		}
		printf("sum2_error = %.7lf,min_error=%.7lf,max_error=%.7lf\n",sum_error,min_error,max_error);

		if(offset_cur%lane_num)
		{
			add_offset[i] = lane_num - (offset_cur%lane_num);
			printf("add %d\n", add_offset[i]);
		}

		fwrite(bias_buf_c,sizeof(int16_t), data_size, fp);

		woffset += offset_cur;
	}

	return 0;
}

int bias_reorder_i16c(float *in,int *offset,int layer_num, float *ap16_range,int *maxQ_array, int lane_num, int *add_offset, FILE *fp)
{
	static int16_t bias_buf_c[MAX_BETA_LENGTH+128];

	int woffset = 0;
	int offset_cur = 0;
	for(int i=0;i<layer_num;i++)
	{
		offset_cur = offset[i];
		printf("Layer %2d;param num=%12d ",i, offset_cur);

		float min,max;
		min = MAX_VALUE;
		max = MIN_VALUE;
		for(int j=0;j<offset_cur;j++)
		{
			float tmp_in_float = in[woffset+j];
			if(tmp_in_float<min)
				min = tmp_in_float;
			if(tmp_in_float>max)
				max = tmp_in_float;
		}
		printf("float min=%.7lf,max=%.7lf ",min,max);//find float min max

		int k;
		int maxQ = -1;
		for(k=0;k<16;k++)//find maxQ
		{
			if(min>ap16_range[2*k]&&max<ap16_range[2*k+1])
			{
				maxQ = k;
			}
			else if(k==0)
			{
				printf("beyond Q0 min=%.7lf,max=%.7lf ",min,max);
				break;
			}
		}
		printf("maxQ=%d ",maxQ);
		maxQ_array[i] = maxQ;		

		double max_error,min_error,sum_error;
		sum_error = 0;
		max_error = MIN_VALUE;
		min_error = MAX_VALUE;

		int data_size = ((offset_cur + lane_num -1)/lane_num)*lane_num;
		for(int j=0;j<data_size;j++)
		{
			float tmp_in_float;
			if(j < offset_cur)
				tmp_in_float = in[woffset+j];
			else
				tmp_in_float = 0;

			int16_t tmp_fixed = (int16_t)(tmp_in_float*pow(2.0,maxQ));
			float tmp_out_float = (float)tmp_fixed*pow(2.0,-maxQ);
			double error = (tmp_out_float - tmp_in_float)*(tmp_out_float - tmp_in_float);
			error = sqrt(error);
			sum_error += error;
			if(error<min_error)
				min_error = error;
			if(error>max_error)
				max_error = error;

			bias_buf_c[j] = tmp_fixed;
		}
		printf("sum2_error = %.7lf,min_error=%.7lf,max_error=%.7lf\n",sum_error,min_error,max_error);

		if(offset_cur%lane_num)
		{
			add_offset[i] = lane_num - (offset_cur%lane_num);
			printf("add %d\n", add_offset[i]);
		}

		fwrite(bias_buf_c,sizeof(int16_t), data_size, fp);

		woffset += offset_cur;
	}

	return 0;
}

int quantize_bias_int16(float *in, short *tmp_out,int *offset,int layer_num, int *ltype_set, float *ap16_range,int *maxQ_array, int *add_offset, FILE *fp)
{
	int i;
	int woffset = 0;
	int offset_cur = 0;
	for(i=0;i<layer_num;i++)
	{
		if(ltype_set[i]==LT_AVGPOOL){
			maxQ_array[i] = 0;
			continue;
		}

		offset_cur = offset[i];
		printf("Layer %2d;param num=%12d ",i, offset_cur);
		int j;
		float min,max;
		min = MAX_VALUE;
		max = MIN_VALUE;
		for(j=0;j<offset_cur;j++)
		{
			float tmp_in_float = in[woffset+j];
			if(tmp_in_float<min)
				min = tmp_in_float;
			if(tmp_in_float>max)
				max = tmp_in_float;
		}
		printf("float min=%.7lf,max=%.7lf ",min,max);//find float min max

		int k;
		int maxQ = -1;
		for(k=0;k<16;k++)//find maxQ
		{
			if(min>ap16_range[2*k]&&max<ap16_range[2*k+1])
			{
				maxQ = k;
			}
			else if(k==0)
			{
				printf("beyond Q0 min=%.7lf,max=%.7lf ",min,max);
				break;
			}
		}
		printf("maxQ=%d ",maxQ);
		maxQ_array[i] = maxQ;

		double max_error,min_error,sum_error;
		sum_error = 0;
		max_error = MIN_VALUE;
		min_error = MAX_VALUE;
		for(j=0;j<offset_cur;j++)
		{
			float tmp_in_float = in[woffset+j];
			short tmp_fixed = (short)(tmp_in_float*pow(2.0,maxQ));
			float tmp_out_float = (float)tmp_fixed*pow(2.0,-maxQ);
			double error = (tmp_out_float - tmp_in_float)*(tmp_out_float - tmp_in_float);
			error = sqrt(error);
			sum_error += error;
			if(error<min_error)
				min_error = error;
			if(error>max_error)
				max_error = error;

			tmp_out[j] = tmp_fixed;
		}
		printf("sum2_error = %.7lf,min_error=%.7lf,max_error=%.7lf",sum_error,min_error,max_error);
		printf("\n");

		if(offset_cur & 0x3)
		{
			add_offset[i] = 4 - (offset_cur & 0x3);
			printf("add %d\n", add_offset[i]);
			for(j=0; j<add_offset[i];j++){
				tmp_out[offset_cur+j]=0;
			}
		}

		fwrite(tmp_out,sizeof(short), offset_cur + add_offset[i], fp);

		woffset += offset_cur;
	}

	return 0;
}

int quantize_ifm_int16(float *in,float *out,int offset,float *ap16_range,int *maxQ_array)
{
	printf("Layer[x]'s ofm; num=%12d ", offset);
	int j;
	float min,max;
	min = MAX_VALUE;
	max = MIN_VALUE;
	for(j=0;j<offset;j++)
	{
		float tmp_in_float = in[j];
		if(tmp_in_float<min)
			min = tmp_in_float;
		if(tmp_in_float>max)
			max = tmp_in_float;
	}
	printf("float min=%.7lf,max=%.7lf ",min,max);//find float min max

	int maxQ = -1;
	for(j=0;j<16;j++)//find maxQ
	{
		if(min>ap16_range[2*j]&&max<ap16_range[2*j+1])
		{
			maxQ = j;
		}
		else if(j==0)
		{
			printf("beyond Q0 min=%.7lf,max=%.7lf ",min,max);
			break;
		}
	}
	
	// int old_maxQ = *maxQ_array;
	// if(old_maxQ!=0)
	// {
	// 	printf("old_maxQ=%d, maxQ=%d ", old_maxQ, maxQ);
	// 	if(maxQ < old_maxQ)
	// 		*maxQ_array = maxQ;
	// }
	// else
	// {
	*maxQ_array = maxQ;
	printf("maxQ=%d ", maxQ);
	// }
	double max_error,min_error,sum_error;
	sum_error = 0;
	max_error = MIN_VALUE;
	min_error = MAX_VALUE;
	for(j=0;j<offset;j++)
	{
		float tmp_in_float = in[j];
		short tmp_fixed = (short)(tmp_in_float*pow(2.0,maxQ));
		float tmp_out_float = (float)tmp_fixed*pow(2.0,-maxQ);
		double error = (tmp_out_float - tmp_in_float)*(tmp_out_float - tmp_in_float);
		error = sqrt(error);
		sum_error += error;
		if(error<min_error)
			min_error = error;
		if(error>max_error)
			max_error = error;

		out[j] = tmp_out_float;
	}
	printf("sum2_error = %.7lf,min_error=%.7lf,max_error=%.7lf",sum_error,min_error,max_error);
	printf("\n");


	return 0;
}


//////////////////////////////////////////quantize_reorder f32_c start

void weight_reorder_tile_f32c(float *weight, float *local_wbuf, int INumxKK, int KK, int TM_MIN, int TN_MIN,
		 int tm, int tn, int ksize, int lane_num, int *add_offset, FILE *fp)
{
	int t1,t2,t3,t4, t5;

	int ma3_d3 = (TM_MIN+lane_num-1)/lane_num;
	float tmp_in_float;
	for(t1 = 0;t1 < ma3_d3; t1++)
	for(t2 = 0;t2 < TN_MIN; t2++)
	for(t3 = 0;t3 <ksize; t3++)
	for(t4 = 0;t4 <ksize; t4++)
	for(t5 = 0;t5 <lane_num; t5++)
	{
		int idx_i = (t1*lane_num+t5+tm)*INumxKK + (t2+tn)*KK + t3*ksize + t4;
		int idx_o = t1*TN_MIN*KK*lane_num + t2*KK*lane_num + t3*ksize*lane_num + t4*lane_num + t5;
		if((t1*lane_num+t5)>=TM_MIN)
			tmp_in_float = 0;
		else
			tmp_in_float = weight[idx_i];

		local_wbuf[idx_o] = tmp_in_float;

	}

	uint16_t val_left = ma3_d3*TN_MIN*KK*lane_num - TM_MIN*TN_MIN*KK;
	*add_offset = *add_offset + val_left;
	fwrite(local_wbuf,sizeof(float), ma3_d3*TN_MIN*KK*lane_num, fp);
}

int weight_reorder_f32c(float *Weight, int IFM_NUM,int OFM_NUM,int Ksize,int TM,int TN, int ltype, int lane_num, int *add_offset, FILE *fp)
{
	const int KxK = Ksize*Ksize;
	const int IFM_NUMxKxK = IFM_NUM*KxK;

	// float local_wbuf[(Tm+3)/4*Tn*K*K*4];
	int lbuf_size = ((Tm+lane_num-1)/lane_num)*Tn*K*K*lane_num;
	float *local_wbuf = (float *)malloc(lbuf_size*sizeof(float));

	int tm,tn,tk;

	int offset_cur = IFM_NUM*OFM_NUM*KxK;
	printf("Layer[x]; param num=%12d ", offset_cur);
	*add_offset = 0;

	int TM_MIN,TN_MIN;
    for(tm = 0; tm < OFM_NUM; tm += TM)
    {
		TM_MIN = MIN_diy(TM,OFM_NUM-tm);
		if(ltype == LT_CONV)
		{
			for(tn = 0; tn < IFM_NUM; tn += TN)
			{
				TN_MIN = MIN_diy(TN, IFM_NUM-tn);
				weight_reorder_tile_f32c(Weight, local_wbuf, IFM_NUMxKxK, KxK, TM_MIN, TN_MIN, tm, tn, Ksize, lane_num, add_offset, fp);
			}
		}else if(ltype == LT_DCONV)
		{	
			weight_reorder_tile_f32c(Weight, local_wbuf, KxK, KxK, TM_MIN,      1, tm,  0, Ksize, lane_num, add_offset, fp);
		}
    }	

	printf("add_offset = %d", *add_offset);
	printf("\n");
	free(local_wbuf);

	return 0;
}

int bias_reorder_i16c(float *in,int *offset,int layer_num, int *ltype_set, float *ap16_range,int *maxQ_array, int *remap_order, int lane_num, int *add_offset, FILE *fp)
{
	static int16_t bias_buf_c[MAX_BETA_LENGTH+128];

	int woffset = 0;
	int offset_cur = 0;
	for(int i=0;i<layer_num;i++)
	{
		int lnum = remap_order[i];
		if((ltype_set[lnum]==LT_AVGPOOL)||(ltype_set[lnum] == LT_MAXPOOL)){
			continue;
		}
		offset_cur = offset[lnum];
		printf("Layer %2d;param num=%12d ",lnum, offset_cur);

		float min,max;
		min = MAX_VALUE;
		max = MIN_VALUE;
		for(int j=0;j<offset_cur;j++)
		{
			float tmp_in_float = in[woffset+j];
			if(tmp_in_float<min)
				min = tmp_in_float;
			if(tmp_in_float>max)
				max = tmp_in_float;
		}
		printf("float min=%.7lf,max=%.7lf ",min,max);//find float min max

		int k;
		int maxQ = -1;
		for(k=0;k<16;k++)//find maxQ
		{
			if(min>ap16_range[2*k]&&max<ap16_range[2*k+1])
			{
				maxQ = k;
			}
			else if(k==0)
			{
				printf("beyond Q0 min=%.7lf,max=%.7lf ",min,max);
				break;
			}
		}
		printf("maxQ=%d ",maxQ);
		maxQ_array[lnum] = maxQ;		

		double max_error,min_error,sum_error;
		sum_error = 0;
		max_error = MIN_VALUE;
		min_error = MAX_VALUE;

		int data_size = ((offset_cur + lane_num -1)/lane_num)*lane_num;
		for(int j=0;j<data_size;j++)
		{
			float tmp_in_float;
			if(j < offset_cur)
				tmp_in_float = in[woffset+j];
			else
				tmp_in_float = 0;

			int16_t tmp_fixed = (int16_t)(tmp_in_float*pow(2.0,maxQ));
			float tmp_out_float = (float)tmp_fixed*pow(2.0,-maxQ);
			double error = (tmp_out_float - tmp_in_float)*(tmp_out_float - tmp_in_float);
			error = sqrt(error);
			sum_error += error;
			if(error<min_error)
				min_error = error;
			if(error>max_error)
				max_error = error;

			bias_buf_c[j] = tmp_fixed;
		}
		printf("sum2_error = %.7lf,min_error=%.7lf,max_error=%.7lf\n",sum_error,min_error,max_error);

		if(offset_cur%lane_num)
		{
			add_offset[lnum] = lane_num - (offset_cur%lane_num);
			printf("add %d\n", add_offset[lnum]);
		}

		fwrite(bias_buf_c,sizeof(int16_t), data_size, fp);

		woffset += offset_cur;
	}

	return 0;
}

int bias_reorder_f32c(float *in,int *offset,int layer_num, int *ltype_set, int *remap_order, int lane_num, int *add_offset, FILE *fp)
{
	static float bias_buf_c[MAX_BETA_LENGTH+128];

	int woffset = 0;
	int offset_cur = 0;
	for(int i=0;i<layer_num;i++)
	{
		int lnum = remap_order[i];
		if((ltype_set[lnum]==LT_AVGPOOL)||(ltype_set[lnum] == LT_MAXPOOL)){
			continue;
		}
		offset_cur = offset[lnum];
		printf("Layer %2d;param num=%12d ",lnum, offset_cur);

		int data_size = ((offset_cur + lane_num -1)/lane_num)*lane_num;
		for(int j=0;j<data_size;j++)
		{
			float tmp_in_float;
			if(j < offset_cur)
				tmp_in_float = in[woffset+j];
			else
				tmp_in_float = 0;

			bias_buf_c[j] = tmp_in_float;
		}

		if(offset_cur%lane_num)
		{
			add_offset[lnum] = lane_num - (offset_cur%lane_num);
			printf("add %d\n", add_offset[lnum]);
		}else{
			printf("\n");
		}

		fwrite(bias_buf_c,sizeof(float), data_size, fp);

		woffset += offset_cur;
	}

	return 0;
}

int bias_reorder_f32c(float *in,int *offset,int layer_num, int *ltype_set, int lane_num, int *add_offset, FILE *fp)
{
	static float bias_buf_c[MAX_BETA_LENGTH+128];

	int woffset = 0;
	int offset_cur = 0;
	for(int i=0;i<layer_num;i++)
	{
		if((ltype_set[i]==LT_AVGPOOL)||(ltype_set[i] == LT_MAXPOOL)){
			continue;
		}
		offset_cur = offset[i];
		printf("Layer %2d;param num=%12d ",i, offset_cur);

		int data_size = ((offset_cur + lane_num -1)/lane_num)*lane_num;
		for(int j=0;j<data_size;j++)
		{
			float tmp_in_float;
			if(j < offset_cur)
				tmp_in_float = in[woffset+j];
			else
				tmp_in_float = 0;

			bias_buf_c[j] = tmp_in_float;
		}

		if(offset_cur%lane_num)
		{
			add_offset[i] = lane_num - (offset_cur%lane_num);
			printf("add %d\n", add_offset[i]);
		}else{
			printf("\n");
		}

		fwrite(bias_buf_c,sizeof(float), data_size, fp);

		woffset += offset_cur;
	}

	return 0;
}


int bias_reorder_f32c(float *in,int *offset,int layer_num, int lane_num, int *add_offset, FILE *fp)
{
	static float bias_buf_c[MAX_BETA_LENGTH+128];

	int woffset = 0;
	int offset_cur = 0;
	for(int i=0;i<layer_num;i++)
	{
		offset_cur = offset[i];
		printf("Layer %2d;param num=%12d ",i, offset_cur);

		int data_size = ((offset_cur + lane_num -1)/lane_num)*lane_num;
		for(int j=0;j<data_size;j++)
		{
			float tmp_in_float;
			if(j < offset_cur)
				tmp_in_float = in[woffset+j];
			else
				tmp_in_float = 0;

			bias_buf_c[j] = tmp_in_float;
		}

		if(offset_cur%lane_num)
		{
			add_offset[i] = lane_num - (offset_cur%lane_num);
			printf("add %d\n", add_offset[i]);
		}else{
			printf("\n");
		}

		fwrite(bias_buf_c,sizeof(float), data_size, fp);

		woffset += offset_cur;
	}

	return 0;
}
//////////////////////////////////////////quantize_reorder f32_c end