
#DEFS = -D REORG_TEST
#all: test

HW_S=2
K=3

Tn=8
Tm=24
Tr=26
Tc=26
MAX_BETA_LENGTH=1024

CFLAGS += -D HW_S=$(HW_S) -D K=$(K) -D Tn=$(Tn) -D Tm=$(Tm)
CFLAGS += -D Tr=$(Tr) -D Tc=$(Tc) -D MAX_BETA_LENGTH=$(MAX_BETA_LENGTH)
CFLAGS += -D LEAKY_RELU
CFLAGS += -D _DEF_IN_MAKEFILE_

LANE_NUM = 8
#INTERQ_MAX = {12....}
INTERQ_MAX = 16
EXTRA_BIT = 10

#for int16
# INTERWIDTH=19

#gen_f32_debug:
#	cp yolov2_acc_f32.h yolov2_acc_sim.h
#	g++ -g -w -o test_layers main.cpp -I . -I ../include -lm -D REORDER_GEN  $(CFLAGS)
diffg:
	diff yolo_last_layer_output.txt golden.txt
cp_bin:
	cp ../01_ExtractWeightAndBiasFromDarknet/darknet/weights.bin .
	cp ../01_ExtractWeightAndBiasFromDarknet/darknet/bias.bin .
cp2hls:
	cp bias_i16c* ../../hls/src_int16_128b
	cp kernel_w_* ../../hls/src_int16_128b	
	cp i16c_* ../../hls/src_int16_128b		
gen_i16c:
	cp acc_f32cti16.h acc.h
	cp yolov2_acc_i16c_gen.h yolov2_acc_sim.h
	g++ -O3 -w -o test_layers main.cpp -I . -lm -D REORDER_GEN  $(CFLAGS) -D LANE_NUM=$(LANE_NUM) -D INTERQ_MAX=$(INTERQ_MAX)
test_i16c:
	cp acc_i16c.h acc.h
	cp yolov2_acc_i16c_test.h yolov2_acc_sim.h
	g++ -O3 -w -o test_layers main.cpp -I . -lm -D REORDER_TEST $(CFLAGS) -D LANE_NUM=$(LANE_NUM) -D EXTRA_BIT=$(EXTRA_BIT)
clean:
	rm yolo_last_layer_output.txt
	rm *.png
	rm test_layers
	
# INTERQ_MAX = 16
# EXTRA_BIT = 0
# [248]:h=0.223988,w=0.224418,x=0.125542,y=0.501358,objectness=0.895340
# [256]:h=0.328778,w=0.076280,x=0.711160,y=0.526275,objectness=0.561455
# [266]:h=0.234967,w=0.109266,x=0.508800,y=0.636938,objectness=0.794704
# [267]:h=0.202159,w=0.121420,x=0.578518,y=0.619786,objectness=0.624294
# [268]:h=0.222028,w=0.085013,x=0.647592,y=0.606460,objectness=0.611614
# [272]:h=0.157134,w=0.125274,x=0.943035,y=0.586639,objectness=0.530480
# [297]:h=0.221379,w=0.096710,x=0.877632,y=0.839949,objectness=0.549882
# vase: 52%
# tvmonitor: 89%
# diningtable: 50%
# chair: 78%
# chair: 58%
# YOLOv2 TEST End

# INTERQ_MAX = 19
# EXTRA_BIT = 0
# [248]:h=0.213523,w=0.217513,x=0.125157,y=0.501641,objectness=0.832917
# [266]:h=0.216885,w=0.105594,x=0.508728,y=0.634842,objectness=0.655997
# [268]:h=0.202554,w=0.087370,x=0.647300,y=0.611225,objectness=0.500244
# tvmonitor: 83%
# chair: 64%
# YOLOv2 TEST End

# INTERQ_MAX = 20
# EXTRA_BIT = 7
# [248]:h=0.214987,w=0.222020,x=0.125367,y=0.499179,objectness=0.808067
# [266]:h=0.228183,w=0.108098,x=0.505244,y=0.628398,objectness=0.617859
# tvmonitor: 81%
# chair: 61%
# YOLOv2 TEST End

# i16_gen
# [248]:h=0.227738,w=0.231994,x=0.127181,y=0.502292,objectness=0.845049
# [266]:h=0.231777,w=0.103552,x=0.511694,y=0.636370,objectness=0.681515
# tvmonitor: 84%
# chair: 67%
# YOLOv2 TEST End

# i16_test
# [248]:h=0.223551,w=0.223980,x=0.125752,y=0.501132,objectness=0.895157
# [256]:h=0.327496,w=0.076354,x=0.710994,y=0.525914,objectness=0.560493
# [266]:h=0.235887,w=0.109372,x=0.508853,y=0.636246,objectness=0.795975
# [267]:h=0.200978,w=0.120357,x=0.578593,y=0.619532,objectness=0.617859
# [268]:h=0.221812,w=0.085847,x=0.647684,y=0.605883,objectness=0.610918
# [272]:h=0.157442,w=0.126256,x=0.942223,y=0.587404,objectness=0.540196
# [297]:h=0.220085,w=0.096521,x=0.877723,y=0.840422,objectness=0.546738
# vase: 51%
# tvmonitor: 89%
# chair: 78%
# chair: 58%
# YOLOv2 TEST End