# YOLOv2 Accelerator in Xilinx's Zynq-7000 Soc(PYNQ-z2, Zedboard and ZCU102)
A Demo for accelerating YOLOv2 in Xilinx's FPGA PYNQ-z2, Zedboard and ZCU102
__I have graduated from Jiangnan University, China in July 1, 2019. Related papers are available now.__  
Master thesis ["Research of Scalability on FPGA-based Neural Network Accelerator"](https://kns.cnki.net/KCMS/detail/detail.aspx?dbcode=CMFD&dbname=CMFDTEMP&filename=1019228234.nh&uid=WEEvREcwSlJHSldRa1FhdXNXaEhoOGhUTzA5T0tESzdFZ2pyR1NJR1ZBaz0=$9A4hF_YAuvQ5obgVAqNKPCYcEjKensW4IQMovwHtwkF4VYPoHbKxJw!!&v=MjE5NTN5dmdXN3JBVkYyNkY3RzZGdFBQcTVFYlBJUjhlWDFMdXhZUzdEaDFUM3FUcldNMUZyQ1VSTE9lWnVkdUY=)  
Journal article ["Design and implementation of FPGA-based deep learning object detection system"](https://kns.cnki.net/KCMS/detail/detail.aspx?dbcode=CJFQ&dbname=CJFDLAST2019&filename=DZJY201908009&uid=WEEvREcwSlJHSldRa1FhdXNXaEhoOGhUTzA5T0tESzdFZ2pyR1NJR1ZBaz0=$9A4hF_YAuvQ5obgVAqNKPCYcEjKensW4IQMovwHtwkF4VYPoHbKxJw!!&v=MDU0NDJDVVJMT2VadWR1Rnl2Z1c3ck1JVGZCZDdHNEg5ak1wNDlGYllSOGVYMUx1eFlTN0RoMVQzcVRyV00xRnI=)   
Journal article ["Design and Implementation of YOLOv2 Accelerator Based on Zynq7000 FPGA Heterogeneous Platform"](https://kns.cnki.net/KCMS/detail/detail.aspx?dbcode=CJFQ&dbname=CJFDTEMP&filename=KXTS201910005&uid=WEEvREcwSlJHSldRa1FhdXNXaEhoOGhUTzA5T0tESzdFZ2pyR1NJR1ZBaz0=$9A4hF_YAuvQ5obgVAqNKPCYcEjKensW4IQMovwHtwkF4VYPoHbKxJw!!&v=MjkwNzdXTTFGckNVUkxPZVp1ZHVGeXZnVzdyT0xqWGZmYkc0SDlqTnI0OUZZWVI4ZVgxTHV4WVM3RGgxVDNxVHI=)

For PYNQ-z2 and Zedboard, in addition to final Linux application( For PYNQ, turn to PYNQ directory; For Zedboard and ZCU102, turn to SDK and PetaLinux), other steps are almost same:
## (1)Software Simulation
Firstly, you should download the darknet source from [https://github.com/pjreddie/darknet](https://github.com/pjreddie/darknet) and yolov2.weights from [https://pjreddie.com/media/files/yolov2.weights](https://pjreddie.com/media/files/yolov2.weights). 


Secondly, modify the darknet's weight load function to get the weights and biases that we want(Here, considering that batcn normalizaton can be combined with weight and bias).

Thirdly, considering that multiple and add operations that implemented in hardware logic will cost too high resources in FPGA[3][6], we should use lower percision operation instead of float-32. Here, I just follow [3] and [6] to quantize the input/output feature maps, weights and biases to dynamic fixed-16. And use fixed-16 operation to replace multiple, add and relu operations in float-32 percision. 

## (2)HLS Accelerator and Simulation
Oh, this part is too complicated to brightly intoduce. __Current design didnt implemment C/RTL simulation, because that testbench always overflow!__ if anyone can solve it, please tell me and upload it, Thanks!
## (3)Vivado Block Design
Just connect the YOLOv2 IP in Vivado Block Design. Only the clock wizzard configuraiotn should be careful. I remembered that input clock is 100MHz, output's clock is 150MHz, __Reset pin active low__, that's all.
## (4)Vivado SDK for Zedboard
This step just wants to get the executable file to driver and control YOLOv2 Acceleraotr in PL. Here, I reserved 0x1000_0000 bytes memories for accelerator to read/wirte feature maps and read weights.
## (5)PetaLinux
Related steps have been updated in Petalinux direcotry. Just use two files(__.hdf file and .bit file__) that generated from Vivado project to create one Peatalinux. Then test yolov2 acclerator in it.

Every directory has some steps to help further implement or study this accelerator.


# Design and Optimization of YOLOv2 Accelerator Based on FPGA  
According to the analysis of the YOLOv2 network, most layers are serially processed, except for the routing layer. The routing layer can be implemented by setting a specific address in advance.   
From an accelerator perspective, the work required is to interact with memory in order (reading memory data, processing data, and then writing back memory data). Since the amount of data input and output is very large, loop tiling technique is always applied to reuse data and reduce memory access times, which tiles the convolution loop R, C, M, N to Tr, Tc, Tm ,Tn[8].  
The overall architecture of the accelerator is shown below:  
![overview](https://github.com/dhm2013724/yolov2_xilinx_fpga/blob/master/overview.png)

Similar to [4,5,8], the accelerator has two AXI4 master interfaces and one AXI4-Lite slave interface. AXI-Lite slave interface is responsible for reading and writing control, data and status register sets. The input feature maps and weights are read concurrently by two master interfaces, and the output feature maps are written back simultaneously through write channel.   
The Data Scatter module is designed to generate the corresponding write address and distribute the data read from the DRAM to the on-chip buffers. The Data Gather module is designed to generate the DRAM write-back address and write the data in the output buffer back to the DRAM. The other red modules are responsible for the processing of the convolutional layer (Conv and Leaky ReLU), the maximum pooling layer (Pool) and the reorg layer (Reorg).  

## Weight Arrangement   
The effective FPGA bandwidth goes up with the increase of burst length and finally flattens out above some burst length threshold[7]. The data tiling technique usually results in a discontinuous DRAM access for the row-major data layout in DRAM. To reduce the number of memory accesses and increase the effective memory bandwidth, we arrange the kernel weights for an entire tile to a continuous block to ensure a high utilization of the bandwidth of external memory [3].  

## Parallel Convolution Engine  
The acceleration strategy of convolutional layer is similar to [5][6], which utilizes input and output parallelism to accelerate the computation. By designing multiple parallel multiplication units and add trees to achieve input parallelism (Tn parallelism) and output parallelism (Tm parallelism) in convolution calculation. The Tm*Tn multiplication units are calculated in parallel. The add trees of Log2 (Tn) depth are accumulated by pipeline, and generate the partial sums.  

## Ping-Pong operation  
Similar to [8], the design implements ping-pong buffers to overlap the delay of reading input feature maps and weights, writing output feature maps and calculation, which greatly improves the dynamic utilization of the computing engines.  

# Evaulate  
Experiments show that floating point addition in HLS requires three DSP resources, floating point multiplication requires two DSPs; fixed point 16-bit multiplication requires one DSP, and fixed-point 16-bit addition can be implemented only using LUT. After placing and routing, resource consumptions of fixed-16 (Tn=2, Tm=32, Tr=26, Tc=26) are shown as follows:     

  |  Resource     |  DSP      | BRAM      | LUT        |  FF        | Freq   | Dev    |
  |  -----        |   -----   | -----     | -----      |  -----     | -----  |-----   |
  |Fixed-16(n4m32)| 153(69%)  | 88(63%)   | 35977(68%) | 36247(34%) |	150MHz |Zedboard|
  |Float-32(n4m23)| 209(95%)  |115(82%)   | 36348(68%) | 64077(60%) |	140MHz |Zedboard|
  |Fixed-16(n4m32)| 147(6%)   | 88(10%)   | 36759(13%) | 30447(6%)  |	180MHz |ZCU102  |  
  
According to the current design, DSP and BRAM are more expensive. The cost of DSP can be further reduced (there are many bit-width redundant multiplications), and the BRAM cost can be reduced. (As Shen [1] said, BRAM allocates an exponential size of 2 in HLS. Actually, many BRAMs are redundant. ).  
The performance comparison in the two cases is shown in the following table:  
  
| Performance              |        |        |        |        |
|  -----                   | -----  | -----  | -----  | -----  |
|CNN models	           |YOLO v2 |YOLO v2 | YOLO v2| YOLO v2|
|Board                     | PYNQ   |Zedboard| ZCU102 |Zedboard|
|Clock(MHz)                |  150   |  150   |  180   |  140   |
|Precision	           |Fixed-16|Fixed-16|Fixed-16|Float-32|
|Power (W)		   |   2.98 |   1.20 | ?      | ?      |
|Operations (GOP)          |29.47   |29.47   | 29.47  | 29.47  |
|Performance(GOP/s)	   |25.98   |30.15   | 36.13  |  6.63  |
|Power Efficiency(GOP/s/W) | 4.20   | 6.02   | ?      | ?      |

# Result  
![image1](https://github.com/dhm2013724/yolov2_xilinx_fpga/blob/150MHzTn4Tm32Tr26Tc26Cin4Cout2/pynq/result2.jpg)

# References:  
[1] Maximizing CNN Accelerator Efficiency Through Resource Partitioning  
[2] PLACID: A Platform for FPGA-Based Accelerator Creation for DCNNs  
[3] Going Deeper with Embedded FPGA Platform for Convolutional Neural Network  
[4] DianNao A Small-Footprint High-Throughput Accelerator for Ubiquitous Machine-Learning  
[5] An Automatic RTL Compiler for High-Throughput FPGA Implementation of Diverse Deep Convolutional Neural Networks  
[6] A Dynamic Multi-precision Fixed-Point Data Quantization Strategy for Convolutional Neural Network  
[7] Caffeine: Towards Uniformed Representation and Acceleration for Deep Convolutional Neural Networks  
[8] Optimizing FPGA-based Accelerator Design for Deep Convolutional Neural Networks  


  
  


