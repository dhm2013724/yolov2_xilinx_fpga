#ifndef __MNV1_H
#define __MNV1_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <sys/ioctl.h>
#include <unistd.h>
#include <linux/fb.h>
#include <time.h>
#include <sys/time.h>
#include <stdint.h>
#include <math.h>

#define XFPGA_ACC_CTRL_BUS_ADDR_AP_CTRL            0x00
#define XFPGA_ACC_CTRL_BUS_ADDR_GIE                0x04
#define XFPGA_ACC_CTRL_BUS_ADDR_IER                0x08
#define XFPGA_ACC_CTRL_BUS_ADDR_ISR                0x0c
#define XFPGA_ACC_CTRL_BUS_ADDR_IFM_DATA           0x10
#define XFPGA_ACC_CTRL_BUS_BITS_IFM_DATA           32
#define XFPGA_ACC_CTRL_BUS_ADDR_OFM_DATA           0x18
#define XFPGA_ACC_CTRL_BUS_BITS_OFM_DATA           32
#define XFPGA_ACC_CTRL_BUS_ADDR_WEIGHT_DATA        0x20
#define XFPGA_ACC_CTRL_BUS_BITS_WEIGHT_DATA        32
#define XFPGA_ACC_CTRL_BUS_ADDR_BIAS_DATA          0x28
#define XFPGA_ACC_CTRL_BUS_BITS_BIAS_DATA          32
#define XFPGA_ACC_CTRL_BUS_ADDR_K_S_PAD_LTYPE_DATA 0x30
#define XFPGA_ACC_CTRL_BUS_BITS_K_S_PAD_LTYPE_DATA 32
#define XFPGA_ACC_CTRL_BUS_ADDR_IOFM_NUM_DATA      0x38
#define XFPGA_ACC_CTRL_BUS_BITS_IOFM_NUM_DATA      32
#define XFPGA_ACC_CTRL_BUS_ADDR_IFM_W_H_DATA       0x40
#define XFPGA_ACC_CTRL_BUS_BITS_IFM_W_H_DATA       32
#define XFPGA_ACC_CTRL_BUS_ADDR_OFM_W_H_DATA       0x48
#define XFPGA_ACC_CTRL_BUS_BITS_OFM_W_H_DATA       32
#define XFPGA_ACC_CTRL_BUS_ADDR_TRTC_DATA          0x50
#define XFPGA_ACC_CTRL_BUS_BITS_TRTC_DATA          32
#define XFPGA_ACC_CTRL_BUS_ADDR_TMTN_DATA          0x58
#define XFPGA_ACC_CTRL_BUS_BITS_TMTN_DATA          32
#define XFPGA_ACC_CTRL_BUS_ADDR_OFM_NUM_BOUND_DATA 0x60
#define XFPGA_ACC_CTRL_BUS_BITS_OFM_NUM_BOUND_DATA 32
#define XFPGA_ACC_CTRL_BUS_ADDR_MLOOPSXTM_DATA     0x68
#define XFPGA_ACC_CTRL_BUS_BITS_MLOOPSXTM_DATA     32
#define XFPGA_ACC_CTRL_BUS_ADDR_MLOOPS_A1XTM_DATA  0x70
#define XFPGA_ACC_CTRL_BUS_BITS_MLOOPS_A1XTM_DATA  32
#define XFPGA_ACC_CTRL_BUS_ADDR_PAD_VAL_DATA       0x78
#define XFPGA_ACC_CTRL_BUS_BITS_PAD_VAL_DATA       16
#define XFPGA_ACC_CTRL_BUS_ADDR_TROWTCOL_DATA      0x80
#define XFPGA_ACC_CTRL_BUS_BITS_TROWTCOL_DATA      32
#define XFPGA_ACC_CTRL_BUS_ADDR_IHW_DATA           0x88
#define XFPGA_ACC_CTRL_BUS_BITS_IHW_DATA           32
#define XFPGA_ACC_CTRL_BUS_ADDR_OHW_DATA           0x90
#define XFPGA_ACC_CTRL_BUS_BITS_OHW_DATA           32
#define XFPGA_ACC_CTRL_BUS_ADDR_KK_INUMXKK_DATA    0x98
#define XFPGA_ACC_CTRL_BUS_BITS_KK_INUMXKK_DATA    32
#define XFPGA_ACC_CTRL_BUS_ADDR_EN_BITS_DATA       0xa0
#define XFPGA_ACC_CTRL_BUS_BITS_EN_BITS_DATA       32
#define XFPGA_ACC_CTRL_BUS_ADDR_WEIGHTQ_DATA       0xa8
#define XFPGA_ACC_CTRL_BUS_BITS_WEIGHTQ_DATA       32
#define XFPGA_ACC_CTRL_BUS_ADDR_BETAQ_DATA         0xb0
#define XFPGA_ACC_CTRL_BUS_BITS_BETAQ_DATA         32
#define XFPGA_ACC_CTRL_BUS_ADDR_INPUTQ_DATA        0xb8
#define XFPGA_ACC_CTRL_BUS_BITS_INPUTQ_DATA        32
#define XFPGA_ACC_CTRL_BUS_ADDR_OUTPUTQ_DATA       0xc0
#define XFPGA_ACC_CTRL_BUS_BITS_OUTPUTQ_DATA       32

//#define ACC_BASEADDR     0x43c00000
//#define WEIGHT_BASEADDR  0x10000000//0x06129EC0 = 101883584 Bytes
//#define BETA_BASEADDR    0x16140000//0x00005414 =     21524 Bytes
//#define MEM_BASEADDR     0x16180000//

#define ACC_BASEADDR     0xA0000000
#define WEIGHT_BASEADDR  0x60000000//0x06129EC0 = 101883584 Bytes
#define BETA_BASEADDR    0x66140000//0x00005414 =     21524 Bytes
#define MEM_BASEADDR     0x66180000//

#define HW_S 2
#define K 3
#define Tn 2
#define Tm 60
#define Tr 26
#define Tc 26
#define MAX_BETA_LENGTH 1024
#define INTERWIDTH 19

#define WriteReg(BaseAddress, RegOffset, Data) *(volatile uint32_t*)((BaseAddress) + (RegOffset)) = (Data)
#define ReadReg(BaseAddress, RegOffset) *(volatile uint32_t*)((BaseAddress) + (RegOffset))

#define MIN_diy(x,y) ((x) < (y) ? (x) : (y))
#define MAX_diy(x,y) ((x) > (y) ? (x) : (y))

#define FALSE 0
#define TRUE  1

#define VALID 0
#define SAME  1

#define LT_CONV  0
#define LT_MAXPOOL 1

#define MIN_NEG (0x8001)

#define OnChipIB_Width  ((Tc-1)*HW_S+K)
#define OnChipIB_Height ((Tr-1)*HW_S+K)

#define HPAGESIZE (2*1024*1024)

void copy_mem2dev(uint8_t *orig,uint32_t byte_num, unsigned long in_buffer);

void copy_dev2mem(uint8_t *dst,uint32_t byte_num, unsigned long in_buffer);

int copy_file2mem(char *bin_file,uint32_t byte_num,unsigned long in_buffer);

int copy_mem2file(char *bin_file,uint32_t byte_num,unsigned long in_buffer);

int FPGA_Acc(uint64_t In_Address, uint64_t Out_Address, uint64_t Weight_offset, uint64_t Beta_offset, uint32_t k_s_pad_ltype, uint32_t iofm_num, uint32_t ifm_w_h, uint32_t ofm_w_h,
	uint32_t TRTC, uint32_t TMTN, int32_t OFM_num_bound, int32_t mLoopsxTM, int32_t mLoops_a1xTM, int16_t pad_val, uint32_t TRowTCol,
	uint32_t IHW, uint32_t OHW, uint32_t KK_INumxKK, uint32_t en_bits, int32_t WeightQ, int32_t BetaQ, int32_t InputQ, int32_t OutputQ);//enable_bits[2:0]={IsReLU, LoadBias, IsNotConv}

#endif
