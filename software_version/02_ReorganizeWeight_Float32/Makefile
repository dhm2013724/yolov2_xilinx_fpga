


test:
	cp yolov2_acc_sim_test.h yolov2_acc_sim.h
	g++ -O3 -w -o test main.cpp -I . -D REORG_TEST -lm $(CFLAGS)
	cp bias.bin ../../hls/src_float32/bias.bin
	cp weights_reorg.bin ../../hls/src_float32/weights_reorg.bin
gen:
	cp yolov2_acc_sim_gen.h yolov2_acc_sim.h
	g++ -O3 -w -o test main.cpp -I . -D REORG_GEN -lm $(CFLAGS)
clean:
	rm *.png
	rm yolov2_acc_sim.h
	rm test

