# Extract weight and bias from Darknet and combined with Batch Normalization
This reop is about how to extract YOLOv2's weight and bias from darknet. 

First, you should download the darknet source from [https://github.com/pjreddie/darknet](https://github.com/pjreddie/darknet) and yolov2.weights from [https://pjreddie.com/media/files/yolov2-voc.weights](https://pjreddie.com/media/files/yolov2-voc.weights). 

Second, add the code segment to function _void load_convolutional_weights(layer l, FILE *fp)_ in _src/parse.c_. This code segment here combines the batch normalization with weights and biases.

Third, make the whole project and run the command:./darknet detect cfg/yolov2.cfg yolov2.weights data/dog.jpg. Then, you will find two bin files (weights.bin and bias.bin), and rename this two files to weightsv2_comb.bin and biasv2_comb.bin for next step.
![two bin files](https://github.com/dhm2013724/yolov2_xilinx_fpga/blob/150MHzTn4Tm32Tr26Tc26Cin4Cout2/Software%20version/01_ExtractWeightAndBiasFromDarknet/s3.jpg)
