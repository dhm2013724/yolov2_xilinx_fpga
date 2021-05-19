# Reorganize/Reorder Weights and Quantize Weights and Biases
This folder is about how to reorganize YOLOv2's weights, and quantize weights and biases from float32 to fixed16. (The number of biases is too small, and we can load each layer's biases just once before the computing phase.)

First, just copy two files(weights.bin and bias.bin) here from step 1.

Then, __make gen_i16; ./test_layers or ./test_layers ../test_imgs/dog.jpg__ This step will generate reorganized weight file(16bit) and some related cfg files.

Last, __make test_i16; ./test_layers__ (The last layer's output is generated from the functon __forward_region_layer__) in yolov2.h. Check the result, and you can turn to hls/src_int16 for hls design. __Dont need to copy any file.__


