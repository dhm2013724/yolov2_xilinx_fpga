# Reorganize/Reorder Weights. FLOAT32_VERSION
This folder is about how to reorganize YOLOv2's weights, and perform the software simulation of YOLOv2 Accelerator in float32 percision. 

Firstly, just copy two files(weights.bin and bias.bin) here from step 1.

Secondly, type __make clean; make gen;__ in shell to gen the executable file. Then, type __./test__ to run the bin file 'test' and generate the __weights_reorg.bin__; 
Off course, you can change the hardware design parameters, like: Tn, Tm, Tr, and Tc. But maybe, there still some bugs for reorg layers to limit these choices. I will solve this problem in recent days!

3rdly, type __make clean; make test;__ in shell to gen the executable file. Then, type __./test__ to test the design;
If you get the same results with above step, you can turn to hls design, just modify some hardware parameters to get your own design. Good Luck.
		

