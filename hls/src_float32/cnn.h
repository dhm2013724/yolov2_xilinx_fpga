
#ifndef _CNN_CPP_
#define _CNN_CPP_

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <ap_int.h>
#include <assert.h>

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


void YOLO2_FPGA(ap_uint<256> *Input,ap_uint<256> *Output,ap_uint<256> *Weight,ap_uint<256> *Beta, int IFM_num, int OFM_num,
							   int Ksize, int Kstride,
							   int Input_w, int Input_h, int Output_w, int Output_h, int Padding, bool IsNL, bool IsBN,
							   int TM, int TN, int TR, int TC,
							   int OFM_num_bound, int mLoopsxTM, int mLoops_a1xTM, int LayerType);

#endif
