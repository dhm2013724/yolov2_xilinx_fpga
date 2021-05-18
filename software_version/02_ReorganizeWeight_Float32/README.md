# Reorganize/Reorder Weights. FLOAT32_VERSION
This folder is about how to reorganize YOLOv2's weights, and perform the software simulation of YOLOv2 Accelerator in float32 percision. 

Firstly, just copy two files(weights.bin and bias.bin) here from step 1.

Then, type __make clean; make gen; ./test__ to run the bin file 'test' and generate __weights_reorg.bin__; 
If you want to change the hardware design parameters(Tn, Tm...), you need to change variables in __hw_params_gen.py__.

Next, type __make clean; make test; ./test__ to test the design;
See Makefile, this step will copy related weights and bias bin files to hls repo. So, you __dont need to copy anything to hls design__. Just checkout the results.

Now, you can turn to hls design, just modify some hardware parameters to get your own design. Good Luck.
		

