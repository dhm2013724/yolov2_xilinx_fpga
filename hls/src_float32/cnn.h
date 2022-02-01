
#ifndef _CNN_CPP_
#define _CNN_CPP_

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <ap_int.h>
#include <assert.h>
#include <hls_stream.h>

#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))

//#define REORG_GEN
//#define REORG_TEST

#define PRAGMA_SUB(x) _Pragma (#x)
#define DO_PRAGMA(x) PRAGMA_SUB(x)

#define S 2
#define K 3
#define MAX_BETA_LENGTH 1024
#define Tn 4
#define Tm 28
#define Tr 26
#define Tc 32
#define OnChipIB_Width 65
#define OnChipIB_Height 53
#define TRow_max 53
#define TCol_max 65


const uint32_t IB_W = OnChipIB_Width;
const uint32_t IB_H = OnChipIB_Height;
const uint32_t IB_HxW = IB_H*IB_W;
//const uint32_t TnxIB_H = Tn*IB_H;
const uint32_t TnxIB_HxIB_W = Tn*IB_H*IB_W;
const uint32_t TrxTc = Tr*Tc;

void YOLO2_FPGA(float *Input,float *Output,float *Weight,float *Beta, int IFM_num, int OFM_num,
							   int Ksize, int Kstride, int Input_w, int Input_h, int Output_w, int Output_h, int Padding, bool IsNL,
							   int TM, int TN, int TR, int TC,
							   int32_t NToy, int32_t NTox, int32_t NTof, int32_t NTcomb, int32_t NTif, uint8_t lmode, int32_t NTcomb_l, int LayerType);

#endif
