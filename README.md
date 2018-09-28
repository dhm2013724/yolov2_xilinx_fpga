# yolov2_xilinx_fpga
A demo for accelerating YOLOv2 in xilinx's fpga PYNQ  
You can follow the step: HLS -> VIVADO -> PYNQ or just jump to PYNQ
Every repo has some steps to help further evaluate or study.  

 
  |  Resource     |  DSP      | BRAM      | LUT        |  FF        | Freq   |
  |  -----        |   -----   | -----     | -----      |  -----     | -----  |
  |Fixed-16(n2m32)| 106(48%)  | 100(72%)  | 27495(52%) | 30118(28%) |	130MHz |
  
  
| Performance              |        |
|  -----                   | -----  |
|CNN models	               |YOLO v2 |
|Board                     | PYNQ   |                
|Clock(MHz)		              |    130 |
|Precision		               |Fixed-16|
|Power (W)		               |   2.71 |
|Operations (GOP)		        |29.47   |
|Performance(GOP/s)		      |11.39   |
|Power Efficiency(GOP/s/W)	|	4.20   |

Result as:  
![image1](https://github.com/dhm2013724/yolov2_xilinx_fpga/blob/master/pynq/result.jpg)

The current design mainly refers to the following papers:  
[1] Maximizing CNN Accelerator Efficiency Through Resource Partitioning  
[2] PLACID: A Platform for FPGA-Based Accelerator Creation for DCNNs  
[3] Going Deeper with Embedded FPGA Platform for Convolutional Neural Network  
[4] DianNao A Small-Footprint High-Throughput Accelerator for Ubiquitous Machine-Learning  
[5] An Automatic RTL Compiler for High-Throughput FPGA Implementation of Diverse Deep Convolutional Neural Networks  
[6] A Dynamic Multi-precision Fixed-Point Data Quantization Strategy for Convolutional Neural Network  
[7] Caffeine: Towards Uniformed Representation and Acceleration for Deep Convolutional Neural Networks  
[8] Optimizing FPGA-based Accelerator Design for Deep Convolutional Neural Networks  


  
  


