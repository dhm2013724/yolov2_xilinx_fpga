import ip from ..\hls\yolov2_ap16_n2m32_inburst\solution1\impl\ip
add ip Yolo2_fpga
add ip ps7.0 apply configuration pynq_revC.tcl
add constraint PYNQ-Z1_C.xdc
create clock_wizard, set output clock 130MHz,set reset type active low
connect ip follow the picture vivado_bd.jpg
generate bitstream
export block_design
get the block_design design_1.tcl and bitstream design_1_wrapper.bit for pynq 
