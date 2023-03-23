# PetaLinux
If you want to create your own Petaliunx verison based on specfic hardware design, pls follow below steps.

## Prepare for Petalinux
First, make sure that you have installed PetaLinux and Vivado in Linux env (See Xilinx's PetaLinux Tools Documentation: Reference Guide (UG1144)) . After that, you should source related settings.sh in order to use petalinux and cross-compiler tool-chain in shell. like these cmds: __source 2019.2/petalinux/settings.sh__ and __source 2019.2/Vivado/2018.2/settings64.sh__. Then, copy two files (hardware description file __design_1_wrapper.hdf__ and bitstream file  __design_1_wrapper.bit__ that generated from Vivado project)to one directory (suggest include the bitstream within hdf file).

## Create PetaLinux Project and modify the Device-tree to reserve memory
Use petalinux-create cmd to create one petalinux project. If you dont know how to create one project, you can type __petalinux-create -h or --help__ to get help and more informations. Here, I just type __petalinux-create -t project -n yolov2 --template zynqMP__ to create a petalinux project named yolov2. Then, __cd yolov2/project-spec/meta-user/recipes-bsp/device-tree/files/__ to find the __system-user.dtsi__ and modify it like below:

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
Third, back to the project directory __yolov2/__, and use type this cmd to initilize your project: __petalinux-config --get-hw-description DIR_where_you_put_the_design_1_wrapper.hdf__. In menu __Subsystem AUTO Hardware Settings__, check your board's UART and SD/SDIO config (some boards have multiple UARTS).
THen, in menu __Image Packaging Configuration__, select __INITRAMFS__ option, and set INITRAMFS/INITRD Image name as __petalinux-image-minimal__. (Here, just follow UG1144).

## SSH config
For SSH, just mutually select openSSH or dropbear, otherwise the ssh's speed would be very slow.

OpenSSH RootFS Configuration.‌‌
$ petalinux-config -c rootfs

Filesystem Packages - console - network - openssh - openssh‌‌ →  Y

Filesystem Packages - console - network - openssh - openssh-*‌‌ →  Y

Image Features - "imagefeature-ssh-server-openssh"‌‌ →  Y

Filesystem Packages - console - network - dropbear‌‌ →  N

Image Features - "imagefeature-ssh-server-dropbear"‌‌ →  N

Filesystem Packages - misc - "packagegroup-core-ssh-dropbear" - "packagegroup-core-ssh-dropbear"‌‌ →  N


DropBear RootFS Configuratoin.
$petalinux-config -c rootfs

Filesystem Packages - console - network - dropbear‌‌ →  Y

Image Features - imagefeature-ssh-server-dropbear‌‌ →  Y

Filesystem Packages - misc - packagegroup-core-ssh-dropbear - packagegroup-core-ssh-dropbear‌‌ →  Y

Filesystem Packages - console - network - openssh - openssh‌‌ →  N

Filesystem Packages - console - network - openssh - openssh-*‌‌ → N

Image Features - imagefeature-ssh-server-openssh‌‌ →  N

## Root login
Then, for some new Petalinux veriosn (2021 or later), default login is not root, if you want to modify it, just follow:
75610 - PetaLinux: How to SSH and SCP to Xilinx Evaluation Boards As The root User Using Dropbear
2021年9月23日•Knowledge
TITLE
75610 - PetaLinux: How to SSH and SCP to Xilinx Evaluation Boards As The root User Using Dropbear
DESCRIPTION
In PetaLinux or Yocto, how can I SSH and SCP to Xilinx Evaluation boards as the root user using dropbear?

SOLUTION
In Yocto or PetaLinux, root login is disabled by default for SSH or SCP.

In order to allow SSH as the root user you need to enable the "debug-tweaks" feature in PetaLinux or Yocto using any of the below methods.

In 2019.2 and prior releases, enabling debug-tweaks also enables auto-login.
This procedure applies to both Xilinx Evaluation and Custom boards.
PetaLinux:

Method 1:

In 2019.2 and prior releases, you can enable debug-tweaks from the petalinux-config options as shown below.


$ petalinux-config ---> Yocto Settings ---> [*] Enable Debug Tweaks

In 2020.1 and later releases, you can enable debug-tweaks from the petalinux-config options as shown below.

$ petalinux-config -c rootfs ---> Image Features ---> [*] debug-tweaks

$ petalinux-config -c rootfs ---> Image Features ---> [*] auto-login

[75610 - PetaLinux: How to SSH and SCP to Xilinx Evaluation Boards As The root User Using Dropbear](https://support.xilinx.com/s/article/75610?language=en_US#:~:text=In%20Yocto%20or%20PetaLinux%2C%20root%20login%20is%20disabled,and%20prior%20releases%2C%20enabling%20debug-tweaks%20also%20enables%20auto-login.)

## build project
After configuraiotn, you can type __petalinux-build__ to build the whole project. After that, type petalinux-build
petalinux-package --boot --format BIN --fsbl images/linux/zynqmp_fsbl.elf --fpga images/linux/system.bit --u-boot --force to generate BOOT.BIN. Finally, you will get three files for petalinux: __BOOT.BIN, image.ub and rootfs.cpio__.

## Partition SD Card and unzip rootfs
Just use Disks tool in ubuntu to partition SD card into two file systems: one FAT and one EXT4. Copy BOOT.BIN and image.ub into FAT file system. (because INITRAM use memory ram as rootfs included in BOOT.BIN)
