This reop is about creating vivado project to get block_design.tcl and bitstream.
Follow below step:
1. import ip from ..\hls\yolov2_ap16_n2m32_inburst\solution1\impl\ip and add ip Yolo2_fpga
2. add ip ps7.0 apply configuration pynq_revC.tcl
   add constraint PYNQ-Z1_C.xdc
3. create clock_wizard, set output clock 150MHz,set reset type active low
4. connect ip follow the picture vivado_bd.jpg
![image]( https://github.com/dhm2013724/yolov2_xilinx_fpga/blob/150MHzTn4Tm32Tr26Tc26Cin4Cout2/vivado/vivado_bd2.jpg)
5. generate bitstream
6. export block_design
7. get the block_design design_1.tcl and bitstream design_1_wrapper.bit for pynq (existed in generated_demo)
8. change two files to the same name, like yolov2.tcl and yolov2.bit for pynq

__please check the vivado 2019.2 prj.(2021-05-19)__
__because int16's design is similar to float32's, here, I didnt copy vivado bd fot int16.__
