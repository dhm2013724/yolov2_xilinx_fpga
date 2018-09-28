# yolov2_xilinx_fpga
A demo for accelerating YOLOv2 in xilinx's fpga PYNQ  
You can follow the step: HLS -> VIVADO -> PYNQ or just jump to PYNQ
Every repo has some steps to help further evaluate or study.  

Current implementation:  

  表头  | 表头
  -- | --
 单元格内容  | 单元格内容
 单元格内容l  | 单元格内容
 
   DSP  | BRAM
  ------------- | -------------
 106(48%)  | 100(72%)
 207(94%)  | 120(86%)
 
    DSP  | BRAM  | LUT
  ------------- | -------------
 106(48%)  | 100(72%)  | 27495(52%)
 207(94%)  | 120(86%)  | 27816(52%)
 
 
    | DSP	 | BRAM  | LUT  | FF  | Freq
 ------------- | ------------- | ------------- | ------------- | ------------- | -------------
  |Fixed-16(n2m32) | 106(48%) |	100(72%) |	27495(52%) |  30118(28%) |	130MHz
  |Float-32(n3m9)	 | 207(94%)	| 120(86%) |	27816(52%) |	27816(26%) |	100MHz

