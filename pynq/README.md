This repo is about PYNQ.  
Weight and bias are available in [here](https://pan.baidu.com/s/1v1U78fdYJ0p8XWmWXA3P0Q)  
Continue with the step of creating yolov2.tcl and yolov2.bit
1. Copy yolov2.tcl and yolov2.bit to this repo
2. Open yolov2.ipynb in jupyter notebook
3. You can change the image path:  
```python
#image path
ORIG_IMG_PATH = 'dog.jpg'
```
4. You can also change the bitstream:(The name of .tcl must be same as the name of .bit)  
```python
overlay = Overlay("yolov2.bit")
```
5. Run all  
# Results as:  
![image1](https://github.com/dhm2013724/yolov2_xilinx_fpga/blob/master/pynq/result.jpg)

