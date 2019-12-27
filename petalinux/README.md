# PetaLinux
Recently, I have tried to implement YOLOv2 Accelerator in PetaLinux v2018.2. And here, I just want to share the experience and PetaLinux create steps to help other people re-implement this accelerator.(In Ubuntu 16.04 with Vivado v2018.2 and PetaLinux 2018.2)

## Prepare for Petalinux
First, make sure that you have installed PetaLinux and Vivado in Ubuntu 16.04(See the official Xilinx manual and User Guide or just google or bing it) . After you have installed PetaLinux and Vivado, you should source related settings.sh in order to use petalinux and cross-compiler tool-chain in shell. like these two cmds: __source /opt/pkg/petalinux/settings.sh__ and __source /opt/Xilinx/Vivado/2018.2/settings64.sh__. Then, copy two files (hardware description file __design_1_wrapper.hdf__ and bitstream file  __design_1_wrapper.bit__ that generated from Vivado project)to one directory.

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

```
/include/ "system-conf.dtsi"
/ {

reserved-memory {
   #address-cells = <2>;
   #size-cells = <2>;
   ranges;
 
   reserved: buffer@0 {
      no-map;
      reg = <0x0 0x60000000 0x0 0x10000000>;
   };
};
 
reserved-driver@0 {
   compatible = "xlnx,reserved-memory";
   memory-region = <&reserved>;
};

};
```
Here, I just follow the wiki to reserve one continued memories for yolov2(Zedboard:base address 0x1000_0000, size 0x1000_0000 bytes, ZCU102: 0x6000_0000, size 0x1000_0000).  
[https://xilinx-wiki.atlassian.net/wiki/spaces/A/pages/18841683/Linux+Reserved+Memory](https://xilinx-wiki.atlassian.net/wiki/spaces/A/pages/18841683/Linux+Reserved+Memory)

## Config project with .hdf file
Third, back to the project directory __yolov2/__, and use type this cmd to initilize your project: __petalinux-config --get-hw-description DIR_where_you_put_the_design_1_wrapper.hdf__. Then, it will come to one menu like below:
![hdf_config.jpg](https://github.com/dhm2013724/yolov2_xilinx_fpga/blob/150MHzTn4Tm32Tr26Tc26Cin4Cout2/petalinux/hdf_config.jpg)

Here, I like to set rootfs from SD card(Just set __Image packaing configration-->file system type or rootfs(I forgot it...)--> SD card__). So that, the files that you used can be saved in SD card when you power-off it. And, then save configuraiton and exit. You will have to wait for a long time(I think that it will donwloads somethings from Internet or just do some configs. Make sure that your computer can access the Internet).

## Further config or build project
If you have other configuraitons, you can use __petalinux-config or petalinux-config -c XXX__ to further config this project. More informations are available in __petalinux-config -h__. After configuraiotn, you can type __petalinux-build__ to build the whole project. After that, type __petalinux-package -boot --fsbl image/linux/zynq_fsbl.elf --fpga --u-boot --force__ to generate BOOT.BIN. Finally, you will get three files for petalinux: __BOOT.BIN, image.ub and rootfs.cpio__.
![package.jpg](https://github.com/dhm2013724/yolov2_xilinx_fpga/blob/150MHzTn4Tm32Tr26Tc26Cin4Cout2/petalinux/package.jpg)

## Partition SD Card and unzip rootfs
Just use Disks tool in ubuntu to partition SD card into two file systems: one FAT and one EXT4. Copy BOOT.BIN and image.ub into FAT file system and copy rootfs.cpio into EXT4 fs.(If you meet permission denied, use sudo). In EXT4 fs, type cmd __sudo pax -rvf rootfs.cpio__ to unzip the rootfs. Then, umount EXT4 and FAT. Here, the Petalinux has been implemented.

## Download related files and Test YOLOv2 Accelerator
Use ssh or other protocol to transmit related files into EXT4 fs, files that yolov2 accelerator needs are available in previous steps. like below:
![yolov2_files.jpg](https://github.com/dhm2013724/yolov2_xilinx_fpga/blob/150MHzTn4Tm32Tr26Tc26Cin4Cout2/petalinux/files_that_yolov2_need.jpg)
type cmd __chmod 777 yolov2_4port_n4m32.elf__ to make it can be executed. Then, type __./yolov2_4port_n4m32.elf one_pic_name.jpg__ to test it. Then, you will get the prediction results and saved as one picture named __predictions.png__.
![output.jpg](https://github.com/dhm2013724/yolov2_xilinx_fpga/blob/150MHzTn4Tm32Tr26Tc26Cin4Cout2/petalinux/output.jpg)
