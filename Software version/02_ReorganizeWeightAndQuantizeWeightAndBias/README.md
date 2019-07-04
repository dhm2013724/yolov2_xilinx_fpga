# Reorganize Weight and Quantize Weight and Bias
This reop is about how to reorganize YOLOv2's weight and quantize weight and bias from float32 to fixed16. (I ran this project in visual studio 2012. If you want to run in other environments, please modify some include headers by yourself.)

First, make sure the marco define __REORG_GEN__ in yolov2.h is used, and annotate the marco define __REORG_TEST__ in yolov2.h and __QUANTI__ in main.cpp. Run this code to reorganize weight in the order of memory access. This step will save the reorganized weight to "weightsv2_comb_reorg.bin".

Second(__Optional__), if you want to test the reorganized weight, you can annotate the marco define __REORG_GEN__ and anti-annotate marco define __REORG_TEST__ in yolov2.h. Run this code to compare the output with First step.


