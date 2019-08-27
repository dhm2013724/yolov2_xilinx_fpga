# PetaLinux
Recently, I have tried to implement YOLOv2 Accelerator in PetaLinux v2018.2. And here, I just want to share the experience and PetaLinux create steps to help other people re-implement this accelerator.(In Ubuntu 16.04 with Vivado v2018.2 and PetaLinux 2018.2)

## Prepare for Petalinux
First, make sure that you have installed PetaLinux and Vivado in Ubuntu 16.04(See the official Xilinx manual and User Guide or just google or bing it) . After you have installed PetaLinux and Vivado, you should source related settings.sh in order to use petalinux and cross-compiler tool-chain in shell. like these two cmds: __source /opt/pkg/petalinux/settings.sh__ and __source /opt/Xilinx/Vivado/2018.2/settings64.sh__. Then, copy two files (hardware description file __design_1_wrapper.hdf__ and bitstream file that generated from Vivado project)to one directory.

## Create PetaLinux Project and modify the Device-tree to reserve memory
Use petalinux-create cmd to create one petalinux project. If you dont know how to create one project, you can type __petalinux-create -h or --help__ to get help and more informations. Here, I just type __petalinux-create -t project -n yolov2 --template zynq__ to create a petalinux project named yolov2. Then, __cd yolov2/project-spec/meta-user/recipes-bsp/device-tree/files/__ to find the __system-user.dtsi__ and modify it like below:

```
/include/ "system-conf.dtsi"
/ {
	reserved-memory {
		#address-cells = <1>;
		#size-cells = <1>;
		ranges;

		reserved: buffer@0x10000000 {
			 no-map;
			 reg = <0x10000000 0x10000000>;
		};
	};

	reserved-driver@0 {
		compatible = "xlnx,reserved-memory";
		memory-region = <&reserved>;
	};
	
};
```
Here, I just follow the wiki to reserve one continued memories for yolov2(base address 0x1000_0000, size 0x1000_0000 bytes). [https://xilinx-wiki.atlassian.net/wiki/spaces/A/pages/18841683/Linux+Reserved+Memory](https://xilinx-wiki.atlassian.net/wiki/spaces/A/pages/18841683/Linux+Reserved+Memory)

