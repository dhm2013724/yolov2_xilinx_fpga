## PetaLinux
Recently, I have tried to implement YOLOv2 Accelerator in PetaLinux v2018.2. And here, I just want to share the experience and PetaLinux create steps to help other people re-implement this accelerator.(In Ubuntu 16.04 with Vivado v2018.2 and PetaLinux 2018.2)

# Prepare for Petalinux
First, make sure that you have installed PetaLinux and Vivado in Ubuntu 16.04(See the official Xilinx manual and User Guide or just google or bing it) . After you have installed PetaLinux and Vivado, you should source related settings.sh to use petalinux and cross-compiler tool-chain in shell. like these two cmds: __source /opt/pkg/petalinux/settings.sh__ and __source /opt/Xilinx/Vivado/2018.2/settings64.sh__. Then, copy two files (hardware description file __design_1_wrapper.hdf__ and bitstream file that generated from Vivado project)to one directory.

# 

