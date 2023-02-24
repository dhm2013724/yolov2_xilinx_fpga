
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
#include <assert.h>
#include <math.h>

#define XFPGA_ACC_CTRL_BUS_ADDR_AP_CTRL            0x00
#define XFPGA_ACC_CTRL_BUS_ADDR_GIE                0x04
#define XFPGA_ACC_CTRL_BUS_ADDR_IER                0x08
#define XFPGA_ACC_CTRL_BUS_ADDR_ISR                0x0c
#define XFPGA_ACC_CTRL_BUS_ADDR_IFM_V_DATA         0x10
#define XFPGA_ACC_CTRL_BUS_BITS_IFM_V_DATA         32
#define XFPGA_ACC_CTRL_BUS_ADDR_OFM_V_DATA         0x18
#define XFPGA_ACC_CTRL_BUS_BITS_OFM_V_DATA         32
#define XFPGA_ACC_CTRL_BUS_ADDR_WEIGHT_V_DATA      0x20
#define XFPGA_ACC_CTRL_BUS_BITS_WEIGHT_V_DATA      32
#define XFPGA_ACC_CTRL_BUS_ADDR_BIAS_V_DATA        0x28
#define XFPGA_ACC_CTRL_BUS_BITS_BIAS_V_DATA        32
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
#define XFPGA_ACC_CTRL_BUS_ADDR_NTOY_DATA          0x60
#define XFPGA_ACC_CTRL_BUS_BITS_NTOY_DATA          32
#define XFPGA_ACC_CTRL_BUS_ADDR_NTOX_DATA          0x68
#define XFPGA_ACC_CTRL_BUS_BITS_NTOX_DATA          32
#define XFPGA_ACC_CTRL_BUS_ADDR_NTOF_DATA          0x70
#define XFPGA_ACC_CTRL_BUS_BITS_NTOF_DATA          32
#define XFPGA_ACC_CTRL_BUS_ADDR_NTCOMB_DATA        0x78
#define XFPGA_ACC_CTRL_BUS_BITS_NTCOMB_DATA        32
#define XFPGA_ACC_CTRL_BUS_ADDR_NTIF_DATA          0x80
#define XFPGA_ACC_CTRL_BUS_BITS_NTIF_DATA          32
#define XFPGA_ACC_CTRL_BUS_ADDR_LMODE_DATA         0x88
#define XFPGA_ACC_CTRL_BUS_BITS_LMODE_DATA         8
#define XFPGA_ACC_CTRL_BUS_ADDR_NTCOMB_L_DATA      0x90
#define XFPGA_ACC_CTRL_BUS_BITS_NTCOMB_L_DATA      32
#define XFPGA_ACC_CTRL_BUS_ADDR_PAD_VAL_DATA       0x98
#define XFPGA_ACC_CTRL_BUS_BITS_PAD_VAL_DATA       16
#define XFPGA_ACC_CTRL_BUS_ADDR_DIV_KK_DATA        0xa0
#define XFPGA_ACC_CTRL_BUS_BITS_DIV_KK_DATA        16
#define XFPGA_ACC_CTRL_BUS_ADDR_IHW_DATA           0xa8
#define XFPGA_ACC_CTRL_BUS_BITS_IHW_DATA           32
#define XFPGA_ACC_CTRL_BUS_ADDR_OHW_DATA           0xb0
#define XFPGA_ACC_CTRL_BUS_BITS_OHW_DATA           32
#define XFPGA_ACC_CTRL_BUS_ADDR_KK_INUMXKK_DATA    0xb8
#define XFPGA_ACC_CTRL_BUS_BITS_KK_INUMXKK_DATA    32
#define XFPGA_ACC_CTRL_BUS_ADDR_EN_BITS_DATA       0xc0
#define XFPGA_ACC_CTRL_BUS_BITS_EN_BITS_DATA       32
#define XFPGA_ACC_CTRL_BUS_ADDR_WEIGHTQ_DATA       0xc8
#define XFPGA_ACC_CTRL_BUS_BITS_WEIGHTQ_DATA       32
#define XFPGA_ACC_CTRL_BUS_ADDR_BIASQ_DATA         0xd0
#define XFPGA_ACC_CTRL_BUS_BITS_BIASQ_DATA         32
#define XFPGA_ACC_CTRL_BUS_ADDR_IFMQ_DATA          0xd8
#define XFPGA_ACC_CTRL_BUS_BITS_IFMQ_DATA          32
#define XFPGA_ACC_CTRL_BUS_ADDR_OFMQ_DATA          0xe0
#define XFPGA_ACC_CTRL_BUS_BITS_OFMQ_DATA          32
#define XFPGA_ACC_CTRL_BUS_ADDR_AVGQ_DATA          0xe8
#define XFPGA_ACC_CTRL_BUS_BITS_AVGQ_DATA          32
#define XFPGA_ACC_CTRL_BUS_ADDR_INTERQ_DATA        0xf0
#define XFPGA_ACC_CTRL_BUS_BITS_INTERQ_DATA        32

#define YOLO2_BASEADDR 0xA0000000
#define WEIGHT_BASE (0x60000000) //203779456 = C25 6D80
#define BETA_BASE (0x6C25F000) //43056 = 0xA830
#define MEM_BASE (0x6C400000) //416*416*32*4+208*208*32*4=173,056+43,264= 216,320*128 = 0x1A6_8000

#define MIN_diy(x,y) ((x) < (y) ? (x) : (y))
#define MAX_diy(x,y) ((x) > (y) ? (x) : (y))

#define FALSE 0
#define TRUE  1

#define VALID 0
#define SAME  1

#define LT_DCONV 0
#define LT_CONV  1
#define LT_AVGPOOL 2
#define LT_MAXPOOL 3

#define MIN_NEG (-1024*1024)
// #define MIN_NEG_INT16 (0x8000)
#define MIN_NEG_INT32 (0x80000000)

#define HW_S 2
#define K 3
#define Tn 8
#define Tm 24
#define Tr 26
#define Tc 26
#define MAX_BETA_LENGTH 1024
#define LANE_NUM 8
#define EXTRA_BIT 10
#define INTERQ_MAX 16
#define EXTRA_BIT 10

#define OnChipIB_Width  ((Tc-1)*HW_S+K)
#define OnChipIB_Height ((Tr-1)*HW_S+K)

#define PRAGMA_SUB(x) _Pragma (#x)
#define DO_PRAGMA(x) PRAGMA_SUB(x)

const uint32_t IB_W = OnChipIB_Width;
const uint32_t IB_H = OnChipIB_Height;
const uint32_t IB_HxW = IB_H*IB_W;
const uint32_t TnxIB_H = Tn*IB_H;	
const uint32_t TnxIB_HxIB_W = Tn*IB_H*IB_W;
const uint32_t TrxTc = Tr*Tc;
const uint32_t Tmax_dx = (Tm+LANE_NUM-1)/LANE_NUM;
const uint32_t Tnax_dx = (Tn+LANE_NUM-1)/LANE_NUM;

#define WriteReg(BaseAddress, RegOffset, Data) *(volatile unsigned int*)((BaseAddress) + (RegOffset)) = (Data)
#define ReadReg(BaseAddress, RegOffset) *(volatile unsigned int*)((BaseAddress) + (RegOffset))

#define HPAGESIZE (4*1024)

void copy_mem2dev(uint8_t *orig,uint32_t byte_num, unsigned long in_buffer);

void copy_dev2mem(uint8_t *dst,uint32_t byte_num, unsigned long in_buffer);

int copy_file2mem(char *bin_file,uint32_t byte_num,unsigned long in_buffer);

int copy_mem2file(char *bin_file,uint32_t byte_num,unsigned long in_buffer);

double what_time_is_it_now();
//typedef ap_uint<16*LANE_NUM> DT_IO;

int FPGA_Acc(unsigned int ifm_addr, unsigned int  ofm_addr, unsigned int weight_offset, unsigned int bias_offset, uint32_t k_s_pad_ltype, uint32_t iofm_num, uint32_t ifm_w_h, uint32_t ofm_w_h,
	uint32_t TRTC, uint32_t TMTN, int32_t NToy, int32_t NTox, int32_t NTof, int32_t NTcomb, int32_t NTif, uint8_t lmode, int32_t NTcomb_l, int16_t pad_val, int16_t div_kk,
	uint32_t IHW, uint32_t OHW, uint32_t KK_INumxKK, uint32_t en_bits, uint32_t weightQ, uint32_t biasQ, uint32_t ifmQ, uint32_t ofmQ, uint32_t avgQ, uint32_t interQ);//enable_bits[2:0]={IsReLU, LoadBias, IsNotConv}
