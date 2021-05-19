

gen_i16:
	python3 hw_params_gen.py
	cp acc_f32.h acc.h
	cp yolov2_acc_i16_gen.h yolov2_acc_sim.h
	g++ -O3 -w -o test_layers main.cpp -I . -I ../include -lm -D REORDER_GEN
test_i16:
	cp acc_i16.h acc.h
	cp yolov2_acc_i16_test.h yolov2_acc_sim.h
	g++ -O3 -w -o test_layers main.cpp -I . -I ../include -lm
	cp bias_ap16* ../../hls/src_int16
	cp weight_reorg_ap16.bin ../../hls/src_int16
	cp weights_ap16_* ../../hls/src_int16
	cp iofm_Q.bin ../../hls/src_int16
clean:
	rm -rf yolo_last_layer_output.txt
	rm -rf *.png
	rm -rf yolov2_acc_sim.h
	rm -rf acc.h
	rm -rf test_layers
	

