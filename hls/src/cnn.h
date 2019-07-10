#ifndef CNN_H

#define CNN_H


void YOLO2_FPGA(int *Input,int *Input1,int *Input2,int *Input3,int *Output,int *Output1,int *Weight,int *Beta,const int InFM_num,const int OutFM_num,
							  const int Kernel_size,const int Kernel_stride,
							  const int Input_w,const int Input_h,const int output_w,const int output_h,const int Padding,const bool IsNL,const bool IsBN,
							  const int TM,const int TN,const int TR,const int TC,
							  const int mLoops,const int nLoops,const int rLoops,const int cLoops,const int LayerType,
							  const int InputQ,const int OutputQ,const int WeightQ,const int BetaQ,int trow_loops);


#endif
