# yolov2_xilinx_fpga
A demo for accelerating YOLOv2 in xilinx's fpga PYNQ  
You can follow the step: HLS -> VIVADO -> PYNQ or just jump to PYNQ
Every repo has some steps to help further evaluate or study.  

Current implementation:  
 
  |               |  DSP      | BRAM      | LUT        |  FF        | Freq   |
  |  -----        |   -----   | -----     | -----      |  -----     | -----  |
  |Fixed-16(n2m32)| 106(48%)  | 100(72%)  | 27495(52%) | 30118(28%) |	130MHz |
  

|  -----                   | -----  |
|CNN models	               |YOLO v2 |
|Board                     | PYNQ   |                
|Clock(MHz)		              |    130 |
|Precision		               |Fixed-16|
|Power (W)		               |   2.71 |
|Operations (GOP)		        |29.47   |
|Performance(GOP/s)		      |11.39   |
|Power Efficiency(GOP/s/W)	|	4.20   |

  
  


