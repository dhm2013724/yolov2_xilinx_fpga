## SDK
These codes are cpmpiled in __Vivado SDK 2018.2__ in __release__ mode use __-O3__ optimization. You should add __-static__ in Liunx g++ linker cmd, and __m__ in Linker->Libraries.
![Release -static](https://github.com/dhm2013724/yolov2_xilinx_fpga/blob/150MHzTn4Tm32Tr26Tc26Cin4Cout2/SDK/release.PNG)

This project is tested in Zedboard and ZCU102 rev1.1. For Zedboard and ZCU102, both reserve 256MB DDR memory for PL.(Zedboard:reserved memory base address:0x1000_0000, memory size 0x1000_0000; ZCU102 start addr 0x6000_0000, size 0x1000_0000). More informations about Linux Reserved Meomory are available in Xilinx's Wiki [Linux Reserved Memory](https://xilinx-wiki.atlassian.net/wiki/spaces/A/pages/18841683/Linux+Reserved+Memory).

