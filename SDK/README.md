## SDK
These codes are cpmpiled in __Vivado SDK 2019.1__ in __release__ mode use __-O2__ optimization. You should add __-static__ in Liunx g++ linker cmd, and __m__ in Linker->Libraries.
![Release -static](release.PNG)

This project is tested in EdgeBoard. For Zynq-7000 and Zynq Ultrascale+ MPSoC, both reserve continous DDR memory region for PL.(Zedboard:reserved memory base address:0x1000_0000, memory size 0x1000_0000; ZCU102/EdgeBoard start addr 0x6000_0000, size 0x1000_0000). More informations about Linux Reserved Meomory are available in Xilinx's Wiki [Linux Reserved Memory](https://xilinx-wiki.atlassian.net/wiki/spaces/A/pages/18841683/Linux+Reserved+Memory).

