 Software Version For YOLOv2
 Considering that Vivado HLS C Simulation is very slow, I usually simulate the code in visual studio 2012 or vscode. From this code, you can also quantize the data precison from float-32 to fixed-16. 
 In addition, you can also change the hardware parameters to design different accelerators.
 
 1.You should find the code with float-32 precison, change the hardware parameters to generate different weights and biases after weight-reorganization.
 2.Find the Macro definition about QUANTI in main.cpp, Undo annotation macro definition and run it to generate fixed-16 weight and biases files.
 3.Annotation macro definition QUANTI, find the code with fixed-16 precison in yolov2.h. Undo annotation and change the hardware parameters to the same that you did. Be careful that if you change Tn, you should change the realted code in compute module and input module.
 
Two files about weight and bias are available [YOLOv2 Float-32 Weight & BIAS in BaiDu CloudDisk](https://pan.baidu.com/s/10XW-u79hx_e-8a7kj9EvAA) Extraction code：f501. 

I have just completed the master's thesis test and will continue to update the entire project after June 5.



