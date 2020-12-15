
#DEFS = -D REORG_TEST
#all: test


test:
	cp yolov2_acc_sim_test.h yolov2_acc_sim.h
	g++ -O3 -w -D REORG_TEST -o test main.cpp -I . -lm
gen:
	cp yolov2_acc_sim_gen.h yolov2_acc_sim.h
	g++ -O3 -w -D REORG_GEN -o test main.cpp -I . -lm
clean:
	rm test

