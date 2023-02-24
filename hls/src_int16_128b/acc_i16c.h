#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <ap_int.h>

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

typedef ap_uint<16*LANE_NUM> DT_IO;

void FPGA_Acc(DT_IO *ifm, DT_IO*ofm, DT_IO *weight, DT_IO* bias, uint32_t k_s_pad_ltype, uint32_t iofm_num, uint32_t ifm_w_h, uint32_t ofm_w_h,
	uint32_t TRTC, uint32_t TMTN, int32_t NToy, int32_t NTox, int32_t NTof, int32_t NTcomb, int32_t NTif, uint8_t lmode, int32_t NTcomb_l, int16_t pad_val, int16_t div_kk,
	uint32_t IHW, uint32_t OHW, uint32_t KK_INumxKK, uint32_t en_bits, uint32_t weightQ, uint32_t biasQ, uint32_t ifmQ, uint32_t ofmQ, uint32_t avgQ, uint32_t interQ);//enable_bits[2:0]={IsReLU, LoadBias, IsNotConv}
