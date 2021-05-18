# Extract weight and bias from Darknet and combined with Batch Normalization
This folder is about how to extract YOLOv2's weight and bias from darknet. 

First, you should download the darknet source from [https://github.com/pjreddie/darknet](https://github.com/pjreddie/darknet) and yolov2.weights from [https://pjreddie.com/media/files/yolov2.weights](https://pjreddie.com/media/files/yolov2.weights). 

Second, add the code segment to function __void load_convolutional_weights(layer l, FILE * fp)__ in __src/parse.c__. This code segment here combines the batch normalization with weights and biases.

Third, make the whole project and run the command: ./darknet detect cfg/yolov2.cfg yolov2.weights data/dog.jpg. 

Then, you will find two bin files (__weights.bin__ and __bias.bin__) that used for next step.(modify the yolov2.cfg width and height __608-> 416__ in next step.)

Currently, we have included the darknet framework and modified this code. You only need to download the yolov2.weights and follow previous step.(2021-05-19)
