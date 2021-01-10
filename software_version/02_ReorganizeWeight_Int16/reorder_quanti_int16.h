
#ifdef REORDER_GEN

#define MIN_VALUE (-1024*1024*1024)
#define MAX_VALUE (1024*1024*1024)

int reorg_quantize_weight_int16(float *Weight, int IFM_NUM,int OFM_NUM,int Ksize,int TM,int TN,float *ap16_range,int *maxQ_array, int *add_offset, FILE *fp)
{
	const int KxK = Ksize*Ksize;
	const int IFM_NUMxKxK = IFM_NUM*KxK;

	int m,n;
	int tm,tn,tk;

	float weight_buffer[Tm*Tn*K*K];
	short weight_buffer2[Tm*Tn*K*K+32];
	
	int offset_cur = IFM_NUM*OFM_NUM*KxK;
	printf("Layer[x]; param num=%12d ", offset_cur);
	*add_offset = 0;

	int j;
	float min,max;
	min = MAX_VALUE;
	max = MIN_VALUE;
	for(j=0;j<offset_cur;j++)
	{
		float tmp_in_float = Weight[j];
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
	*maxQ_array = maxQ;

	double max_error,min_error,sum_error;
	sum_error = 0;
	max_error = MIN_VALUE;
	min_error = MAX_VALUE;

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

			int TN_MINxTM_MIN = TN_MIN*TM_MIN;
			for(tk = 0;tk < KxK; tk++)
				for(tm = 0;tm < TM_MIN; tm++)
					for(tn = 0;tn < TN_MIN;tn++)
					{
						//weight_buffer2[tk*TN_MINxTM_MIN + tm*TN_MIN + tn] = weight_buffer[tm*TN_MIN*KxK + tn*KxK + tk];
						float tmp_in_float = weight_buffer[tm*TN_MIN*KxK + tn*KxK + tk];
						short tmp_fixed = (short)(tmp_in_float*pow(2.0,maxQ));
						float tmp_out_float = (float)tmp_fixed*pow(2.0,-maxQ);
						double error = (tmp_out_float - tmp_in_float)*(tmp_out_float - tmp_in_float);
						error = sqrt(error);
						sum_error += error;
						if(error<min_error)
							min_error = error;
						if(error>max_error)
							max_error = error;

						weight_buffer2[tk*TN_MINxTM_MIN + tm*TN_MIN + tn] = tmp_fixed;
					}
			
			int local_offset = TM_MIN*TN_MIN*KxK;
			if(local_offset & 0x1)
			{
				*add_offset = *add_offset + 1;
			}
			fwrite(weight_buffer2,sizeof(short), local_offset + (local_offset & 0x1), fp);
			//memcpy((float *)(Weight_reorg+offset),weight_buffer2,TM_MIN*TN_MIN*KxK*sizeof(float));
		}							
	}

	printf("sum2_error = %.7lf,min_error=%.7lf,max_error=%.7lf,add_offset = %d",sum_error,min_error,max_error, *add_offset);
	printf("\n");

	return 0;
}

int quantize_bias_int16(float *in, short *tmp_out,int *offset,int layer_num,float *ap16_range,int *maxQ_array, int *add_offset, FILE *fp)
{
	int i;
	int offset_index = 0;
	int woffset = 0;
	int offset_cur = 0;
	for(i=0;i<layer_num;i++)
	{
		offset_cur = offset[offset_index];
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

		if(offset_cur & 0x1)
		{
			printf("add 1\n");
			tmp_out[offset_cur]=0;
			add_offset[i] = 1;
		}else
			add_offset[i] = 0;

		fwrite(tmp_out,sizeof(short), offset_cur + add_offset[i], fp);

		woffset += offset_cur;
		offset_index++;
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
	
	int old_maxQ = *maxQ_array;
	if(old_maxQ!=0)
	{
		printf("old_maxQ=%d, maxQ=%d ", old_maxQ, maxQ);
		if(maxQ < old_maxQ)
			*maxQ_array = maxQ;
	}
	else
	{
		*maxQ_array = maxQ;
		printf("maxQ=%d ", maxQ);
	}
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


#endif

