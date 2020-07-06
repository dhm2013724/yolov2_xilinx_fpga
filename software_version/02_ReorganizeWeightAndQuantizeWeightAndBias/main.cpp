
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <fcntl.h>
#include <string.h>
#include <time.h>
#include "yolov2.h"

#define MIN_VALUE (-1024*1024*1024)
#define MAX_VALUE (1024*1024*1024)

#define QUANTI

#ifndef QUANTI
int main( int argc, char *argv[])
{
	//freopen("result.txt","w",stdout);
	printf("YOLOv2 TEST Begin\n");
    	char **names = get_labels("coco.names");
	int x;
	for(x=0;x<80;x++)//80 classe labels
	{
		printf("[%d]%s\n",x,names[x]);
	}
    	image **alphabet = load_alphabet();
    	network *net = load_network("yolov2.cfg");
	set_batch_network(net, 1);

////////////////////load img resize img begin
	char img_buff[256];
	char *input_imgfn = img_buff;
	if(argc==1)
		strncpy(input_imgfn, "../test_imgs/dog.jpg", 256);
	else
		strncpy(input_imgfn, argv[1], 256);
	image im = load_image_stb(input_imgfn, 3);//3 channel img
	printf("Input img:%s\n w=%d,h=%d,c=%d\n", input_imgfn, im.w, im.h, im.c);
	image sized = letterbox_image(im, 416, 416);
	save_image_png(sized, "sized");// convert to yolov3 net input size 416x416x3
////////////////////load img resize img end

	time_t first, second;       
	layer l = net->layers[net->n-1];
    	float *X = sized.data;

	first=time(NULL); 
	yolov2_hls_ps(net, X);
	second=time(NULL); 
	printf("%s: Predicted in %f seconds.\n", input_imgfn, difftime(second,first));

    int nboxes = 0;
    float nms=.45;
	float thresh = .5;
	float hier_thresh = .5;
    detection *dets = get_network_boxes(net, im.w, im.h, thresh, hier_thresh, 0, 1, &nboxes);
    printf("%d\n", nboxes);
	//for(x=0;x<nboxes;x++)
	//{
	//	printf("[%3d]:h=%f,w=%f,x=%f,y=%f,objectness=%f\n",x,dets[x].bbox.h,dets[x].bbox.w,dets[x].bbox.x,dets[x].bbox.y,dets[x].objectness);
	//}

    if (nms) do_nms_sort(dets, nboxes, l.classes, nms);
    draw_detections(im, dets, nboxes, thresh, names, alphabet, l.classes);

    free_detections(dets, nboxes);
	
///////////////////write predictions img
	save_image_png(im, "predictions");// output

	free_image(im);
    free_image(sized);

	printf("YOLOv2 TEST End\n");

    return 0;
}

#else

int quantize_short16(float *in,short *out,int *offset,int layer_num,float *ap16_range,int *maxQ_array)
{
	int i;
	int offset_index = 0;
	int woffset = 0;
	for(i=0;i<layer_num;i++)
	{
		if(offset[offset_index]==0)
			return i;
		printf("Layer %2d;weight num=%12d ",i,offset[offset_index]);
		int j;
		float min,max;
		min = MAX_VALUE;
		max = MIN_VALUE;
		for(j=0;j<offset[offset_index];j++)
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
		for(j=0;j<offset[offset_index];j++)
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

			out[woffset+j] = tmp_fixed;
		}
		printf("sum2_error = %.7lf,min_error=%.7lf,max_error=%.7lf",sum_error,min_error,max_error);
		printf("\n");

		woffset += offset[offset_index];
		offset_index++;
	}

	return 0;
}

int main(int argc,char *argv[])
{
	int i;
	printf("Test fixed-point\n");

	int weight_offset[32] = {864, 18432, 73728, 8192, 73728,
		294912, 32768, 294912, 1179648, 131072, 1179648, 131072,
		1179648, 4718592, 524288, 4718592, 524288, 4718592, 9437184,
		9437184, 32768, 11796480, 435200, 0, 0, 0, 0, 0, 0, 0, 0, 0};

	int beta_offset[32] = {32, 64, 128, 64, 128, 256, 128, 256, 512, 256, 512, 256, 512, 1024,
		512, 1024, 512, 1024, 1024, 1024, 64, 1024, 425, 0, 0, 0, 0, 0, 0, 0, 0, 0};

	short *Weight_fixed_buf = (short *)calloc(203767168/4,sizeof(short));
	float *Weight_buf = (float *)calloc(203767168/4,sizeof(float));
	float *Beta_buf   = (float *)calloc(43044/4,sizeof(float));

	FILE *fp_w = fopen("weights_reorg.bin", "rb");
    if(!fp_w) printf("fopen weights_reorg.bin error\n");
	FILE *fp_b = fopen("bias.bin", "rb");
    if(!fp_b) printf("fopen bias.bin error\n");

	fread(Weight_buf, sizeof(float), 203767168/4, fp_w);
	fread(Beta_buf, sizeof(float), 43044/4, fp_b);

	fclose(fp_w);
	fclose(fp_b);
////////////////////////////////

	short ap16_min = 0x8000;
	short ap16_max = 0x7fff;
	printf("ap16_min = %d \nap16_max = %d\n",ap16_min,ap16_max);
	float ap16_range[16*2];
	for(i=0;i<16;i++)
	{
		printf("Q%2d:",i);
		ap16_range[2*i]   = (float)ap16_min*pow((float)2,-i);//min
		ap16_range[2*i+1] = (float)ap16_max*pow((float)2,-i);//max
		printf("min=%.7lf,max=%.7lf\n",ap16_range[2*i],ap16_range[2*i+1]);
	}
////////////////////////////////
	int maxQ_array[32];
	int layer_num;
	FILE* fout;
	char layer_num_string[256];
	char s[256];

	printf("weight quantize begin\n");
	layer_num = quantize_short16(Weight_buf,Weight_fixed_buf,weight_offset,32,ap16_range,maxQ_array);
	for(i=0;i<layer_num;i++)
	{
		printf("[%d]=%d\n",i,maxQ_array[i]);
	}
	sprintf(s,"weights_reorg_ap16_maxQ_%d.bin", layer_num);
	printf("%s\n",s);

	fout = fopen(s,"wb");
    if(!fout) printf("fopen %s error\n",s);
	fwrite(maxQ_array,sizeof(int), layer_num,fout);
	fclose(fout);

	fout = fopen("weights_reorg_ap16.bin","wb");
    if(!fout) printf("fopen weights_reorg_ap16.bin error\n");
	fwrite(Weight_fixed_buf,sizeof(short), 203767168/4,fout);
	fclose(fout);
	printf("weight quantize end\n");

	printf("beta quantize begin\n");
	layer_num = quantize_short16(Beta_buf,Weight_fixed_buf,beta_offset,32,ap16_range,maxQ_array);
	for(i=0;i<layer_num;i++)
	{
		printf("[%d]=%d\n",i,maxQ_array[i]);
	}
	sprintf(s,"bias_ap16_maxQ_%d.bin", layer_num);
	printf("%s\n",s);

	fout = fopen(s,"wb");
    if(!fout) printf("fopen %s error\n",s);
	fwrite(maxQ_array,sizeof(int), layer_num,fout);
	fclose(fout);

	fout = fopen("bias_ap16.bin","wb");
    if(!fout) printf("fopen bias_ap16.bin error\n");
	fwrite(Weight_fixed_buf,sizeof(short), 43044/4+1,fout);
	fclose(fout);
	printf("beta quantize end\n");

	free(Weight_fixed_buf);
	free(Weight_buf);
	free(Beta_buf);
	
	printf("0.1=%.10lf\n",((short)(0.1*pow(2.0,15)))*pow(2.0,-15));
	printf("0.1=%x\n",(short)(0.1*pow(2.0,15)));
	return 0;
}

#endif


