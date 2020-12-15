# Reorganize/Reorder Weights and Quantize Weights and Biases
This folder is about how to reorganize YOLOv2's weights, and quantize weights and biases from float32 to fixed16. (The number of biases is too small, and we can load each layer's biases just once before the computing phase.)

Firstly, just copy two files(weights.bin and bias.bin) here from step 1.

Secondly, __make gen; ./test_layers or ./test_layers ../test_imgs/dog.jpg__ This step will generate reorganized weight file: weights_reorg.bin.

3rdly(__Optional__), if you want to test the reorganized weight, you can __make test; ./test_layers or ./test_layers ../test_imgs/dog.jpg__ And, compare the output with 2nd's output. (The last layer's output is generated from the functon __forward_region_layer__) in yolov2.h.

4thly, if you want to quantize the wieghts and bias, __make quanti; ./test_layers or ./test_layers ../test_imgs/dog.jpg__. This step will generate four files: "bias_ap16.bin",  "bias_ap16_maxQ_23.bin", "weights_reorg_ap16.bin", and "weights_reorg_ap16_maxQ_23.bin". Two bin files(*_ap16.bin) are quantized weights and bias. The other two files(* _ maxQ_23.bin) are fraction bit for weights and bias in 23 convolutional layers.

