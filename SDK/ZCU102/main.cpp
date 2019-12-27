
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <fcntl.h>
#include <string.h>
#include <time.h>
#include "yolov2.h"

int main( int argc, char *argv[])
{
	unsigned int WEIGHT_BASE = 0x60000000;//0x 0612 9EC0
    unsigned int BETA_BASE = 0x66400000;// 0x 5414
    unsigned int MEM_BASE  = 0x66800000;// 49770*1024*4 = 203857920 = C26A000

	printf("YOLOv2 TEST Begin\n");
    char **names = get_labels("coco.names");
	int x;
	for(x=0;x<80;x++)//80 classe labels
	{
		printf("[%d]%s\n",x,names[x]);
	}
    image **alphabet = load_alphabet();
    network *net = load_network("yolov2.cfg", "yolov2.weights", 0);
    set_batch_network(net, 1);

	//load_weights_hls(net, "yolov2.weights", 0, net->n);//write weight to file
////////////////////load img resize img begin
	char buff[256];
    char *input_imgfn = buff;
    if(argc==1)
    {
    	strncpy(input_imgfn, "person.jpg", 256);
//        printf("WEIGHT_BASE = %d|%x\n",WEIGHT_BASE,WEIGHT_BASE);
//        printf("BETA_BASE = %d|%x\n",BETA_BASE,BETA_BASE);
//        printf("MEM_BASE = %d|%x\n",MEM_BASE,MEM_BASE);
    }
    else
    {
    	strncpy(input_imgfn, argv[1], 256);
//    	WEIGHT_BASE = atoi(argv[2]);
//        BETA_BASE = atoi(argv[3]);
//        MEM_BASE  = atoi(argv[4]);
//        printf("WEIGHT_BASE = %d|%x\n",WEIGHT_BASE,WEIGHT_BASE);
//        printf("BETA_BASE = %d|%x\n",BETA_BASE,BETA_BASE);
//        printf("MEM_BASE = %d|%x\n",MEM_BASE,MEM_BASE);
    }
	printf("Input img:%s\n",input_imgfn);
	image im = load_image_stb(input_imgfn, 3);//3 channel img
	printf("img w=%d,h=%d,c=%d\n",im.w,im.h,im.c);
	image sized = letterbox_image(im, 416, 416);
	save_image_png(sized, "sized");// convert to yolov3 net input size 416x416x3
////////////////////load img resize img end

//	time_t first, second;
	double time;
	layer l = net->layers[net->n-1];
    float *X = sized.data;
//	first=time(NULL);
    time = what_time_is_it_now();
    //network_predict(net, X);
	yolov2_hls_ps(net, X,WEIGHT_BASE,BETA_BASE,MEM_BASE);
//	second=time(NULL);
//	printf("%s: Predicted in %f seconds.\n", input_imgfn, difftime(second,first));
	printf("Predicted in %f seconds.\n",what_time_is_it_now()-time);

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
