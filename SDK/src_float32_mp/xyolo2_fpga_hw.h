// ==============================================================
// Vivado(TM) HLS - High-Level Synthesis from C, C++ and SystemC v2019.2 (64-bit)
// Copyright 1986-2019 Xilinx, Inc. All Rights Reserved.
// ==============================================================
// CTRL_BUS
// 0x00 : Control signals
//        bit 0  - ap_start (Read/Write/COH)
//        bit 1  - ap_done (Read/COR)
//        bit 2  - ap_idle (Read)
//        bit 3  - ap_ready (Read)
//        bit 7  - auto_restart (Read/Write)
//        others - reserved
// 0x04 : Global Interrupt Enable Register
//        bit 0  - Global Interrupt Enable (Read/Write)
//        others - reserved
// 0x08 : IP Interrupt Enable Register (Read/Write)
//        bit 0  - Channel 0 (ap_done)
//        bit 1  - Channel 1 (ap_ready)
//        others - reserved
// 0x0c : IP Interrupt Status Register (Read/TOW)
//        bit 0  - Channel 0 (ap_done)
//        bit 1  - Channel 1 (ap_ready)
//        others - reserved
// 0x10 : Data signal of Input_r
//        bit 31~0 - Input_r[31:0] (Read/Write)
// 0x14 : reserved
// 0x18 : Data signal of Output_r
//        bit 31~0 - Output_r[31:0] (Read/Write)
// 0x1c : reserved
// 0x20 : Data signal of Weight
//        bit 31~0 - Weight[31:0] (Read/Write)
// 0x24 : reserved
// 0x28 : Data signal of Beta
//        bit 31~0 - Beta[31:0] (Read/Write)
// 0x2c : reserved
// 0x30 : Data signal of IFM_num
//        bit 31~0 - IFM_num[31:0] (Read/Write)
// 0x34 : reserved
// 0x38 : Data signal of OFM_num
//        bit 31~0 - OFM_num[31:0] (Read/Write)
// 0x3c : reserved
// 0x40 : Data signal of Ksize
//        bit 31~0 - Ksize[31:0] (Read/Write)
// 0x44 : reserved
// 0x48 : Data signal of Kstride
//        bit 31~0 - Kstride[31:0] (Read/Write)
// 0x4c : reserved
// 0x50 : Data signal of Input_w
//        bit 31~0 - Input_w[31:0] (Read/Write)
// 0x54 : reserved
// 0x58 : Data signal of Input_h
//        bit 31~0 - Input_h[31:0] (Read/Write)
// 0x5c : reserved
// 0x60 : Data signal of Output_w
//        bit 31~0 - Output_w[31:0] (Read/Write)
// 0x64 : reserved
// 0x68 : Data signal of Output_h
//        bit 31~0 - Output_h[31:0] (Read/Write)
// 0x6c : reserved
// 0x70 : Data signal of Padding
//        bit 31~0 - Padding[31:0] (Read/Write)
// 0x74 : reserved
// 0x78 : Data signal of IsNL
//        bit 0  - IsNL[0] (Read/Write)
//        others - reserved
// 0x7c : reserved
// 0x80 : Data signal of IsBN
//        bit 0  - IsBN[0] (Read/Write)
//        others - reserved
// 0x84 : reserved
// 0x88 : Data signal of TM
//        bit 31~0 - TM[31:0] (Read/Write)
// 0x8c : reserved
// 0x90 : Data signal of TN
//        bit 31~0 - TN[31:0] (Read/Write)
// 0x94 : reserved
// 0x98 : Data signal of TR
//        bit 31~0 - TR[31:0] (Read/Write)
// 0x9c : reserved
// 0xa0 : Data signal of TC
//        bit 31~0 - TC[31:0] (Read/Write)
// 0xa4 : reserved
// 0xa8 : Data signal of OFM_num_bound
//        bit 31~0 - OFM_num_bound[31:0] (Read/Write)
// 0xac : reserved
// 0xb0 : Data signal of mLoopsxTM
//        bit 31~0 - mLoopsxTM[31:0] (Read/Write)
// 0xb4 : reserved
// 0xb8 : Data signal of mLoops_a1xTM
//        bit 31~0 - mLoops_a1xTM[31:0] (Read/Write)
// 0xbc : reserved
// 0xc0 : Data signal of LayerType
//        bit 31~0 - LayerType[31:0] (Read/Write)
// 0xc4 : reserved
// (SC = Self Clear, COR = Clear on Read, TOW = Toggle on Write, COH = Clear on Handshake)

#ifndef _YOLOv2_HW_H
#define _YOLOv2_HW_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <sys/ioctl.h>
#include <unistd.h>
#include <linux/fb.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <sys/time.h>

#define XYOLO2_FPGA_CTRL_BUS_ADDR_AP_CTRL        0x000
#define XYOLO2_FPGA_CTRL_BUS_ADDR_GIE            0x004
#define XYOLO2_FPGA_CTRL_BUS_ADDR_IER            0x008
#define XYOLO2_FPGA_CTRL_BUS_ADDR_ISR            0x00c
#define XYOLO2_FPGA_CTRL_BUS_ADDR_INPUT0_DATA    0x010
#define XYOLO2_FPGA_CTRL_BUS_BITS_INPUT0_DATA    32
#define XYOLO2_FPGA_CTRL_BUS_ADDR_INPUT1_DATA    0x018
#define XYOLO2_FPGA_CTRL_BUS_BITS_INPUT1_DATA    32
#define XYOLO2_FPGA_CTRL_BUS_ADDR_INPUT2_DATA    0x020
#define XYOLO2_FPGA_CTRL_BUS_BITS_INPUT2_DATA    32
#define XYOLO2_FPGA_CTRL_BUS_ADDR_INPUT3_DATA    0x028
#define XYOLO2_FPGA_CTRL_BUS_BITS_INPUT3_DATA    32
#define XYOLO2_FPGA_CTRL_BUS_ADDR_OUTPUT0_DATA   0x030
#define XYOLO2_FPGA_CTRL_BUS_BITS_OUTPUT0_DATA   32
#define XYOLO2_FPGA_CTRL_BUS_ADDR_OUTPUT1_DATA   0x038
#define XYOLO2_FPGA_CTRL_BUS_BITS_OUTPUT1_DATA   32
#define XYOLO2_FPGA_CTRL_BUS_ADDR_WEIGHT0_DATA   0x040
#define XYOLO2_FPGA_CTRL_BUS_BITS_WEIGHT0_DATA   32
#define XYOLO2_FPGA_CTRL_BUS_ADDR_WEIGHT1_DATA   0x048
#define XYOLO2_FPGA_CTRL_BUS_BITS_WEIGHT1_DATA   32
#define XYOLO2_FPGA_CTRL_BUS_ADDR_WEIGHT2_DATA   0x050
#define XYOLO2_FPGA_CTRL_BUS_BITS_WEIGHT2_DATA   32
#define XYOLO2_FPGA_CTRL_BUS_ADDR_WEIGHT3_DATA   0x058
#define XYOLO2_FPGA_CTRL_BUS_BITS_WEIGHT3_DATA   32
#define XYOLO2_FPGA_CTRL_BUS_ADDR_BETA_DATA      0x060
#define XYOLO2_FPGA_CTRL_BUS_BITS_BETA_DATA      32
#define XYOLO2_FPGA_CTRL_BUS_ADDR_IFM_NUM_DATA   0x068
#define XYOLO2_FPGA_CTRL_BUS_BITS_IFM_NUM_DATA   32
#define XYOLO2_FPGA_CTRL_BUS_ADDR_OFM_NUM_DATA   0x070
#define XYOLO2_FPGA_CTRL_BUS_BITS_OFM_NUM_DATA   32
#define XYOLO2_FPGA_CTRL_BUS_ADDR_KSIZE_DATA     0x078
#define XYOLO2_FPGA_CTRL_BUS_BITS_KSIZE_DATA     32
#define XYOLO2_FPGA_CTRL_BUS_ADDR_KSTRIDE_DATA   0x080
#define XYOLO2_FPGA_CTRL_BUS_BITS_KSTRIDE_DATA   32
#define XYOLO2_FPGA_CTRL_BUS_ADDR_INPUT_W_DATA   0x088
#define XYOLO2_FPGA_CTRL_BUS_BITS_INPUT_W_DATA   32
#define XYOLO2_FPGA_CTRL_BUS_ADDR_INPUT_H_DATA   0x090
#define XYOLO2_FPGA_CTRL_BUS_BITS_INPUT_H_DATA   32
#define XYOLO2_FPGA_CTRL_BUS_ADDR_OUTPUT_W_DATA  0x098
#define XYOLO2_FPGA_CTRL_BUS_BITS_OUTPUT_W_DATA  32
#define XYOLO2_FPGA_CTRL_BUS_ADDR_OUTPUT_H_DATA  0x0a0
#define XYOLO2_FPGA_CTRL_BUS_BITS_OUTPUT_H_DATA  32
#define XYOLO2_FPGA_CTRL_BUS_ADDR_PADDING_DATA   0x0a8
#define XYOLO2_FPGA_CTRL_BUS_BITS_PADDING_DATA   32
#define XYOLO2_FPGA_CTRL_BUS_ADDR_ISNL_DATA      0x0b0
#define XYOLO2_FPGA_CTRL_BUS_BITS_ISNL_DATA      1
#define XYOLO2_FPGA_CTRL_BUS_ADDR_TM_DATA        0x0b8
#define XYOLO2_FPGA_CTRL_BUS_BITS_TM_DATA        32
#define XYOLO2_FPGA_CTRL_BUS_ADDR_TN_DATA        0x0c0
#define XYOLO2_FPGA_CTRL_BUS_BITS_TN_DATA        32
#define XYOLO2_FPGA_CTRL_BUS_ADDR_TR_DATA        0x0c8
#define XYOLO2_FPGA_CTRL_BUS_BITS_TR_DATA        32
#define XYOLO2_FPGA_CTRL_BUS_ADDR_TC_DATA        0x0d0
#define XYOLO2_FPGA_CTRL_BUS_BITS_TC_DATA        32
#define XYOLO2_FPGA_CTRL_BUS_ADDR_NTOY_DATA      0x0d8
#define XYOLO2_FPGA_CTRL_BUS_BITS_NTOY_DATA      32
#define XYOLO2_FPGA_CTRL_BUS_ADDR_NTOX_DATA      0x0e0
#define XYOLO2_FPGA_CTRL_BUS_BITS_NTOX_DATA      32
#define XYOLO2_FPGA_CTRL_BUS_ADDR_NTOF_DATA      0x0e8
#define XYOLO2_FPGA_CTRL_BUS_BITS_NTOF_DATA      32
#define XYOLO2_FPGA_CTRL_BUS_ADDR_NTCOMB_DATA    0x0f0
#define XYOLO2_FPGA_CTRL_BUS_BITS_NTCOMB_DATA    32
#define XYOLO2_FPGA_CTRL_BUS_ADDR_NTIF_DATA      0x0f8
#define XYOLO2_FPGA_CTRL_BUS_BITS_NTIF_DATA      32
#define XYOLO2_FPGA_CTRL_BUS_ADDR_LMODE_DATA     0x100
#define XYOLO2_FPGA_CTRL_BUS_BITS_LMODE_DATA     8
#define XYOLO2_FPGA_CTRL_BUS_ADDR_NTCOMB_L_DATA  0x108
#define XYOLO2_FPGA_CTRL_BUS_BITS_NTCOMB_L_DATA  32
#define XYOLO2_FPGA_CTRL_BUS_ADDR_LAYERTYPE_DATA 0x110
#define XYOLO2_FPGA_CTRL_BUS_BITS_LAYERTYPE_DATA 32


//#define XYOLO2_FPGA_CTRL_BUS_ADDR_AP_CTRL        0x000
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_GIE            0x004
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_IER            0x008
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_ISR            0x00c
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_INPUT0_DATA    0x010
//#define XYOLO2_FPGA_CTRL_BUS_BITS_INPUT0_DATA    32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_INPUT1_DATA    0x018
//#define XYOLO2_FPGA_CTRL_BUS_BITS_INPUT1_DATA    32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_INPUT2_DATA    0x020
//#define XYOLO2_FPGA_CTRL_BUS_BITS_INPUT2_DATA    32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_INPUT3_DATA    0x028
//#define XYOLO2_FPGA_CTRL_BUS_BITS_INPUT3_DATA    32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_OUTPUT_R_DATA  0x030
//#define XYOLO2_FPGA_CTRL_BUS_BITS_OUTPUT_R_DATA  32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_WEIGHT0_DATA   0x038
//#define XYOLO2_FPGA_CTRL_BUS_BITS_WEIGHT0_DATA   32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_WEIGHT1_DATA   0x040
//#define XYOLO2_FPGA_CTRL_BUS_BITS_WEIGHT1_DATA   32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_WEIGHT2_DATA   0x048
//#define XYOLO2_FPGA_CTRL_BUS_BITS_WEIGHT2_DATA   32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_WEIGHT3_DATA   0x050
//#define XYOLO2_FPGA_CTRL_BUS_BITS_WEIGHT3_DATA   32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_BETA_DATA      0x058
//#define XYOLO2_FPGA_CTRL_BUS_BITS_BETA_DATA      32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_IFM_NUM_DATA   0x060
//#define XYOLO2_FPGA_CTRL_BUS_BITS_IFM_NUM_DATA   32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_OFM_NUM_DATA   0x068
//#define XYOLO2_FPGA_CTRL_BUS_BITS_OFM_NUM_DATA   32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_KSIZE_DATA     0x070
//#define XYOLO2_FPGA_CTRL_BUS_BITS_KSIZE_DATA     32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_KSTRIDE_DATA   0x078
//#define XYOLO2_FPGA_CTRL_BUS_BITS_KSTRIDE_DATA   32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_INPUT_W_DATA   0x080
//#define XYOLO2_FPGA_CTRL_BUS_BITS_INPUT_W_DATA   32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_INPUT_H_DATA   0x088
//#define XYOLO2_FPGA_CTRL_BUS_BITS_INPUT_H_DATA   32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_OUTPUT_W_DATA  0x090
//#define XYOLO2_FPGA_CTRL_BUS_BITS_OUTPUT_W_DATA  32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_OUTPUT_H_DATA  0x098
//#define XYOLO2_FPGA_CTRL_BUS_BITS_OUTPUT_H_DATA  32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_PADDING_DATA   0x0a0
//#define XYOLO2_FPGA_CTRL_BUS_BITS_PADDING_DATA   32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_ISNL_DATA      0x0a8
//#define XYOLO2_FPGA_CTRL_BUS_BITS_ISNL_DATA      1
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_TM_DATA        0x0b0
//#define XYOLO2_FPGA_CTRL_BUS_BITS_TM_DATA        32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_TN_DATA        0x0b8
//#define XYOLO2_FPGA_CTRL_BUS_BITS_TN_DATA        32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_TR_DATA        0x0c0
//#define XYOLO2_FPGA_CTRL_BUS_BITS_TR_DATA        32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_TC_DATA        0x0c8
//#define XYOLO2_FPGA_CTRL_BUS_BITS_TC_DATA        32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_NTOY_DATA      0x0d0
//#define XYOLO2_FPGA_CTRL_BUS_BITS_NTOY_DATA      32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_NTOX_DATA      0x0d8
//#define XYOLO2_FPGA_CTRL_BUS_BITS_NTOX_DATA      32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_NTOF_DATA      0x0e0
//#define XYOLO2_FPGA_CTRL_BUS_BITS_NTOF_DATA      32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_NTCOMB_DATA    0x0e8
//#define XYOLO2_FPGA_CTRL_BUS_BITS_NTCOMB_DATA    32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_NTIF_DATA      0x0f0
//#define XYOLO2_FPGA_CTRL_BUS_BITS_NTIF_DATA      32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_LMODE_DATA     0x0f8
//#define XYOLO2_FPGA_CTRL_BUS_BITS_LMODE_DATA     8
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_NTCOMB_L_DATA  0x100
//#define XYOLO2_FPGA_CTRL_BUS_BITS_NTCOMB_L_DATA  32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_LAYERTYPE_DATA 0x108
//#define XYOLO2_FPGA_CTRL_BUS_BITS_LAYERTYPE_DATA 32


//#define XYOLO2_FPGA_CTRL_BUS_ADDR_AP_CTRL        0x00
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_GIE            0x04
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_IER            0x08
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_ISR            0x0c
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_INPUT_R_DATA   0x10
//#define XYOLO2_FPGA_CTRL_BUS_BITS_INPUT_R_DATA   32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_OUTPUT_R_DATA  0x18
//#define XYOLO2_FPGA_CTRL_BUS_BITS_OUTPUT_R_DATA  32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_WEIGHT0_DATA   0x20
//#define XYOLO2_FPGA_CTRL_BUS_BITS_WEIGHT0_DATA   32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_WEIGHT1_DATA   0x28
//#define XYOLO2_FPGA_CTRL_BUS_BITS_WEIGHT1_DATA   32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_WEIGHT2_DATA   0x30
//#define XYOLO2_FPGA_CTRL_BUS_BITS_WEIGHT2_DATA   32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_WEIGHT3_DATA   0x38
//#define XYOLO2_FPGA_CTRL_BUS_BITS_WEIGHT3_DATA   32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_BETA_DATA      0x40
//#define XYOLO2_FPGA_CTRL_BUS_BITS_BETA_DATA      32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_IFM_NUM_DATA   0x48
//#define XYOLO2_FPGA_CTRL_BUS_BITS_IFM_NUM_DATA   32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_OFM_NUM_DATA   0x50
//#define XYOLO2_FPGA_CTRL_BUS_BITS_OFM_NUM_DATA   32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_KSIZE_DATA     0x58
//#define XYOLO2_FPGA_CTRL_BUS_BITS_KSIZE_DATA     32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_KSTRIDE_DATA   0x60
//#define XYOLO2_FPGA_CTRL_BUS_BITS_KSTRIDE_DATA   32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_INPUT_W_DATA   0x68
//#define XYOLO2_FPGA_CTRL_BUS_BITS_INPUT_W_DATA   32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_INPUT_H_DATA   0x70
//#define XYOLO2_FPGA_CTRL_BUS_BITS_INPUT_H_DATA   32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_OUTPUT_W_DATA  0x78
//#define XYOLO2_FPGA_CTRL_BUS_BITS_OUTPUT_W_DATA  32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_OUTPUT_H_DATA  0x80
//#define XYOLO2_FPGA_CTRL_BUS_BITS_OUTPUT_H_DATA  32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_PADDING_DATA   0x88
//#define XYOLO2_FPGA_CTRL_BUS_BITS_PADDING_DATA   32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_ISNL_DATA      0x90
//#define XYOLO2_FPGA_CTRL_BUS_BITS_ISNL_DATA      1
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_TM_DATA        0x98
//#define XYOLO2_FPGA_CTRL_BUS_BITS_TM_DATA        32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_TN_DATA        0xa0
//#define XYOLO2_FPGA_CTRL_BUS_BITS_TN_DATA        32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_TR_DATA        0xa8
//#define XYOLO2_FPGA_CTRL_BUS_BITS_TR_DATA        32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_TC_DATA        0xb0
//#define XYOLO2_FPGA_CTRL_BUS_BITS_TC_DATA        32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_NTOY_DATA      0xb8
//#define XYOLO2_FPGA_CTRL_BUS_BITS_NTOY_DATA      32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_NTOX_DATA      0xc0
//#define XYOLO2_FPGA_CTRL_BUS_BITS_NTOX_DATA      32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_NTOF_DATA      0xc8
//#define XYOLO2_FPGA_CTRL_BUS_BITS_NTOF_DATA      32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_NTCOMB_DATA    0xd0
//#define XYOLO2_FPGA_CTRL_BUS_BITS_NTCOMB_DATA    32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_NTIF_DATA      0xd8
//#define XYOLO2_FPGA_CTRL_BUS_BITS_NTIF_DATA      32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_LMODE_DATA     0xe0
//#define XYOLO2_FPGA_CTRL_BUS_BITS_LMODE_DATA     8
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_NTCOMB_L_DATA  0xe8
//#define XYOLO2_FPGA_CTRL_BUS_BITS_NTCOMB_L_DATA  32
//#define XYOLO2_FPGA_CTRL_BUS_ADDR_LAYERTYPE_DATA 0xf0
//#define XYOLO2_FPGA_CTRL_BUS_BITS_LAYERTYPE_DATA 32


//#define YOLO2_BASEADDR 0x43c00000
//#define WEIGHT_BASE (0x10000000) //203767168 = C253D80
//#define BETA_BASE (0x1C25F000) //43044 = 0xA824
//#define MEM_BASE (0x1C26A000) //416*416*32*4+208*208*32*4=173,056+43,264= 216,320*128 = 0x1A6_8000

#define YOLO2_BASEADDR 0xA0000000
#define WEIGHT_BASE (0x60000000) //203767168 = C253D80
#define BETA_BASE (0x6C25F000) //43044 = 0xA824
#define MEM_BASE (0x6C400000) //416*416*32*4+208*208*32*4=173,056+43,264= 216,320*128 = 0x1A6_8000

#define Tn 4
#define Tm 36
#define Tr 26
#define Tc 32

#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))
#define S 2
#define K 3

#define OnChipIB_Width  ((Tc-1)*S+K)
#define OnChipIB_Height ((Tr-1)*S+K)
#define MAX_BETA_LENGTH (1024)

const uint32_t IB_W = OnChipIB_Width;
const uint32_t IB_H = OnChipIB_Height;
const uint32_t IB_HxW = IB_H*IB_W;
const uint32_t TnxIB_H = Tn*IB_H;
const uint32_t TnxIB_HxIB_W = Tn*IB_H*IB_W;
const uint32_t TrxTc = Tr*Tc;


#define WriteReg(BaseAddress, RegOffset, Data) *(volatile unsigned int*)((BaseAddress) + (RegOffset)) = (Data)
#define ReadReg(BaseAddress, RegOffset) *(volatile unsigned int*)((BaseAddress) + (RegOffset))

#define HPAGESIZE (4*1024)

void copy_mem2dev(uint8_t *orig,uint32_t byte_num, unsigned long in_buffer);

void copy_dev2mem(uint8_t *dst,uint32_t byte_num, unsigned long in_buffer);

int copy_file2mem(char *bin_file,uint32_t byte_num,unsigned long in_buffer);

int copy_mem2file(char *bin_file,uint32_t byte_num,unsigned long in_buffer);

void copy_mem2dev(uint8_t *orig,uint32_t byte_num, unsigned long in_buffer)
{
	int fd = open("/dev/mem", O_RDWR);
	unsigned char *virtual_addr;
	uint32_t RequestByteNum;// must page
	if(byte_num%(HPAGESIZE)==0)
		RequestByteNum = byte_num;
	else
	{
		RequestByteNum = ceil(byte_num/(HPAGESIZE*1.0))*(HPAGESIZE);
	}
	virtual_addr = (unsigned char *)mmap(NULL, RequestByteNum, PROT_READ | PROT_WRITE, MAP_SHARED, fd, (off_t)in_buffer);
	if(virtual_addr == MAP_FAILED)
	{
		perror("Virtual_addr_in mappong for absolute memory access failed!\n");
		return;
	}
	memcpy(virtual_addr,orig,byte_num);

	munmap((void *)virtual_addr, byte_num);
	close(fd);
}

void copy_dev2mem(uint8_t *dst,uint32_t byte_num, unsigned long in_buffer)
{
	int fd = open("/dev/mem", O_RDWR);
	unsigned char *virtual_addr;
	uint32_t RequestByteNum;// must page
	if(byte_num%(HPAGESIZE)==0)
		RequestByteNum = byte_num;
	else
	{
		RequestByteNum = ceil(byte_num/(HPAGESIZE*1.0))*(HPAGESIZE);
	}
		virtual_addr = (unsigned char *)mmap(NULL, RequestByteNum, PROT_READ | PROT_WRITE, MAP_SHARED, fd, (off_t)in_buffer);
	if(virtual_addr == MAP_FAILED)
	{
		perror("Virtual_addr_in mappong for absolute memory access failed!\n");
		return;
	}
	printf("copy start-----byte_num=%d\n",byte_num);
	memcpy((uint8_t *)dst,virtual_addr,byte_num);
	printf("copy ok!\n");

	munmap((void *)virtual_addr, byte_num);
	close(fd);
}

int copy_file2mem(char *bin_file,uint32_t byte_num,unsigned long in_buffer)
{
	unsigned char *buffer = (unsigned char *)malloc(HPAGESIZE);
	if(buffer==NULL){
		printf("cannot malloc buffer %d byte\n", HPAGESIZE);
		return -1;
	}
	printf("Total Byte Num = %d\n Address 0x%X\n", byte_num, in_buffer);
	FILE *fp;
	if( (fp = fopen(bin_file, "rb")) == NULL)fprintf(stderr,"CANNOT OPEN bin_file\n");
	int rd_num;
	unsigned long offset = 0;
	while(rd_num = fread(buffer, sizeof(unsigned char), HPAGESIZE, fp))
	{
		if(rd_num < HPAGESIZE)
			rd_num = HPAGESIZE;
		copy_mem2dev(buffer,rd_num, in_buffer+offset);
//		printf("rd_num=%d, offset=%d\n", rd_num, offset);
		offset += rd_num;
	}
	printf("copy_file2mem offset=%d\n",offset);
	fclose(fp);

	free(buffer);


	return 0;
}

int copy_mem2file(char *bin_file,uint32_t byte_num,unsigned long in_buffer)
{
	void *buffer = malloc(HPAGESIZE);
	if(buffer==NULL){
		printf("cannot malloc buffer %d byte\n", HPAGESIZE);
		return -1;
	}

	FILE *fp;
	if( (fp = fopen(bin_file, "wb")) == NULL)fprintf(stderr,"CANNOT OPEN bin_file\n");

	int x = byte_num;
	int addbyte;
	unsigned long offset = 0;
	while(addbyte=((x<HPAGESIZE)?x:(HPAGESIZE)))
	{
		copy_dev2mem((uint8_t *)buffer,addbyte, in_buffer+offset);
		fwrite(buffer , sizeof(unsigned char), addbyte, fp);
		x -= addbyte;
		offset += addbyte;
	}
	printf("copy_mem2file offset=%d\n",offset);


	fclose(fp);

	free(buffer);

	return 0;
}



#endif

