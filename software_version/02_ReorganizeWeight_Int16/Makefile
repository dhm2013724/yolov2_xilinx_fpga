
#DEFS = -D REORG_TEST
#all: test

HW_S=2
K=3

Tn=2
Tm=60
Tr=26
Tc=26
MAX_BETA_LENGTH=1024

CFLAGS += -D HW_S=$(HW_S) -D K=$(K) -D Tn=$(Tn) -D Tm=$(Tm)
CFLAGS += -D Tr=$(Tr) -D Tc=$(Tc) -D MAX_BETA_LENGTH=$(MAX_BETA_LENGTH)

#for int16
INTERWIDTH=19

gen_i16:
	cp acc_f32.h acc.h
	cp yolov2_acc_i16_gen.h yolov2_acc_sim.h
	g++ -O3 -w -o test_layers main.cpp -I . -I ../include -lm -D REORDER_GEN  $(CFLAGS)
test_i16:
	cp acc_i16.h acc.h
	cp yolov2_acc_i16_test.h yolov2_acc_sim.h
	g++ -O3 -w -o test_layers main.cpp -I . -I ../include -lm $(CFLAGS) -D INTERWIDTH=$(INTERWIDTH)
clean:
	rm yolo_last_layer_output.txt
	rm *.png
	rm yolov2_acc_sim.h
	rm acc.h
	rm test_layers
	

