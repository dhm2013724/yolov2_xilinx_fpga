
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
	printf("YOLOv2 TEST Begin\n");
    char **names = get_labels("coco.names");
	int x;

    image **alphabet = load_alphabet();
    network *net = load_network("yolov2.cfg");
	set_batch_network(net, 1);

	char img_buff[256];
	char *input_imgfn = img_buff;
	if(argc==1)
		strncpy(input_imgfn, "../test_imgs/000000000139.jpg", 256);
	else
		strncpy(input_imgfn, argv[1], 256);
	image im = load_image_stb(input_imgfn, 3);//3 channel img
	printf("Input img:%s\n w=%d,h=%d,c=%d\n", input_imgfn, im.w, im.h, im.c);
	image sized = letterbox_image(im, 416, 416);
	// save_image_png(sized, "sized");// convert to yolov3 net input size 416x416x3

	time_t first, second;
	first=time(NULL); 
	yolov2_hls_ps(net, sized.data);
	second=time(NULL); 
	printf("%s: Predicted in %f seconds.\n", input_imgfn, difftime(second,first));

	int nboxes = 0;
    float nms=.45;
	float thresh = .5;
	float hier_thresh = .5;
	detection *dets = get_network_boxes(net, im.w, im.h, thresh, hier_thresh, 0, 1, &nboxes);
	printf("%d\n", nboxes);
	for(x=0;x<nboxes;x++)
	{
		if(dets[x].objectness > 0.0)
			printf("[%3d]:h=%f,w=%f,x=%f,y=%f,objectness=%f\n",x,dets[x].bbox.h,dets[x].bbox.w,dets[x].bbox.x,dets[x].bbox.y,dets[x].objectness);
	}
	// printf("class_num = %d\n", l.classes);
	int class_num = 80;
	if (nms) do_nms_sort(dets, nboxes, class_num, nms);
	draw_detections(im, dets, nboxes, thresh, names, alphabet, class_num);

	free_detections(dets, nboxes);	
///////////////////write predictions img
	save_image_png(im, "predictions");// output

	free_image(im);
	free_image(sized);
	printf("YOLOv2 TEST End\n");

	return 0;
}

