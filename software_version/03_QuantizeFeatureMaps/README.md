## Quantize Input/Output Feature Maps
Compile and run this code, it will quantize each layer's input feature maps. For each layer, the previous layer's output fearure maps are current layer's input feature maps. Besides of 23 convolutional layers, 5 MaxPool layers and 1 Reorg layer will not affect the feature maps' percisions. And for layer 0, considering that resize operation in preprocess phase maps RGB pixel values from 0-255(uint8) to 0-1(float), I just set input feature maps' fraction bit to be 14.

Besides of the files included in this folder, other files can be generated from previous step(just 4 generated files in step 2). This step will generate one bin file("yolov2_ap16_inout_maxQ_24.bin") about fraction bit for 23 convolutional layers.(Each MaxPool Layer and Reorg Layer's input/output feature maps' fraction bit will share previous layer's output feature maps' fraction bit.)

About the QNUM's selection, here just for simple, you can make related codes to change the QNUM for percision. I just changed QNUM's value to test the percision, and choose 23. QNUM would only affect the percision of CONV layer and the design of CONV compute module.

Defaultly, quantization stage would perfer to use only one pic __kite.jpg__ in test_imgs. Here just use one pic to get IFM/OFMs' quantized factors, you can also use more pictures to make this stage more logical.
