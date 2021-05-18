## YOLOv2 Accelerator
This repo is about YOLOv2 accelerator implemented in vivado HLS 2019.2. 

__Target Device:  xc7z020clg484-1__ 

__Target clock: 5.20__ (Estimated: 5.675 Uncertainty: 0.65; __I just ran it in 150MHz__) 
![overview](https://github.com/dhm2013724/yolov2_xilinx_fpga/blob/150MHzTn4Tm32Tr26Tc26Cin4Cout2/hls/c-syn.PNG)

Current testbench can just pass C-simulation, tb for C-RTL cosimulation maybe need several days. 
Some other related files are available from software version foler.

![overview](https://github.com/dhm2013724/yolov2_xilinx_fpga/blob/150MHzTn4Tm32Tr26Tc26Cin4Cout2/hls/c_sim.PNG)

(I just ran this in vivado HLS 2019.2 in 6th July 2020, but the C Simulation completed after nearly one and half hour; I will optimize the code in the future.)

This step will use five bin files that generated from software version's step 2&3.(just some weights, bias, and related quantized params)

above descriptions are related to src_int16. I recently tested the accelerator in float32 percision on Zedboard, and code is in src__float32.
__Target Device:  xc7z020clg484-1__ 
__Target clock: 4.00__

But it can only ran maximun clock in 140MHz, if you dont further optimize this design. (Off course, you can decrease some useless bitwidth for some variables, and so on.)
