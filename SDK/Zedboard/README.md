## SDK
These codes are cpmpiled in __Vivado SDK 2018.2__ in __release__ mode use __-O3__ optimization. And you should add '-static' string in Liunx g++ linker cmd
![Release -static](https://github.com/dhm2013724/yolov2_xilinx_fpga/blob/150MHzTn4Tm32Tr26Tc26Cin4Cout2/SDK/release.PNG)

This project is tested in Zedboard, and reserve 256MB DDR memory for PL.(in this design, reserved memory base address:0x1000_0000, memory size 0x1000_0000). More informations about Linux Reserved Meomory are available in Xilinx's Wiki [Linux Reserved Memory](https://xilinx-wiki.atlassian.net/wiki/spaces/A/pages/18841683/Linux+Reserved+Memory).

