


test:
	cp yolov2_acc_sim_test.h yolov2_acc_sim.h
	g++ -O3 -w -o test main.cpp -I . -D REORG_TEST -lm $(CFLAGS)
	cp bias.bin ../../hls/src_float32_mp/bias.bin
	cp weights_reorg.bin ../../hls/src_float32_mp/weights_reorg.bin
# cp bias.bin ../../hls/src_float32_fusion/bias.bin
# cp weights_reorg.bin ../../hls/src_float32_fusion/weights_reorg.bin	
gen:
	python3 hw_params_gen.py
	cp yolov2_acc_sim_gen.h yolov2_acc_sim.h
	g++ -O3 -w -o test main.cpp -I . -D REORG_GEN -lm $(CFLAGS)
cp2hls:
	cp bias.bin ../../hls/src_float32/bias.bin
	cp weights_reorg.bin ../../hls/src_float32/weights_reorg.bin
cp_bin:
	cp ../01_ExtractWeightAndBiasFromDarknet/darknet/weights.bin .
	cp ../01_ExtractWeightAndBiasFromDarknet/darknet/bias.bin .
clean_bin:
	rm -rf *.bin
clean:
	rm -rf *.png
	rm -rf yolov2_acc_sim.h
	rm -rf test

#//darknet origion result
# mask_scale: Using default '1.000000'
# Loading weights from yolov2.weights...Done!
# ../../test_imgs/000000000139.jpg: Predicted in 4.127025 seconds.
# 845
# [248]:h=0.227699,w=0.232511,x=0.127223,y=0.502272,objectness=0.846462
# [266]:h=0.231706,w=0.103453,x=0.511701,y=0.636352,objectness=0.683597
# tvmonitor: 85%
# chair: 67%

# make test result
# 845
# [248]:h=0.227699,w=0.232511,x=0.127223,y=0.502272,objectness=0.846463
# [266]:h=0.231706,w=0.103453,x=0.511701,y=0.636352,objectness=0.683597
# tvmonitor: 85%
# chair: 67%
# YOLOv2 TEST End