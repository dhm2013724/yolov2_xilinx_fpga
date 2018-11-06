
#ifndef YOLOV2_HLS

#define YOLOV2_HLS

#include <stdio.h>
#include <string.h>
#include <ap_int.h>

//#define MAX(x,y) ((x)>(y)?(x):(y))
//#define MIN(x,y) ((x)<(y)?(x):(y))
//#define S 2
//#define K 3
//
//#define Tn 2
//#define Tm 32
//#define Tr 26
//#define Tc 26
//#define OnChipIB_Width  ((Tc-1)*S+K)
//#define OnChipIB_Height ((Tr-1)*S+K)
//#define MAX_BETA_LENGTH (1024)
//#define INTER_WIDTH (19)

void YOLO2_FPGA(int *Input,int *Output,int *Weight,int *Beta,const int InFM_num,const int OutFM_num,
							  const int Kernel_size,const int Kernel_stride,
							  const int Input_w,const int Input_h,const int Padding,const bool IsNL,const bool IsBN,
							  const int TM,const int TN,const int TR,const int TC,
							  const int mLoops,const int nLoops,const int rLoops,const int cLoops,const int LayerType,
							  const int InputQ,const int OutputQ,const int WeightQ,const int BetaQ);

#endif
