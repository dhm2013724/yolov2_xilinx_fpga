# import numpy as np
import sys,os

S = 2
K = 3
MAX_BETA_LENGTH = 1024

Tn = 4
Tm = 28
Tr = 26
Tc = 32
OnChipIB_Width = (Tc-1)*S+K
OnChipIB_Height = (Tr-1)*S+K
TRow_max = OnChipIB_Height
TCol_max = OnChipIB_Width

head_define_str = ""
head_define_str += "#define S " + str(S) + "\n"
head_define_str += "#define K " + str(K) + "\n"
head_define_str += "#define MAX_BETA_LENGTH " + str(MAX_BETA_LENGTH) + "\n"
head_define_str += "#define Tn " + str(Tn) + "\n"
head_define_str += "#define Tm " + str(Tm) + "\n"
head_define_str += "#define Tr " + str(Tr) + "\n"
head_define_str += "#define Tc " + str(Tc) + "\n"
head_define_str += "#define OnChipIB_Width " + str(OnChipIB_Width) + "\n"
head_define_str += "#define OnChipIB_Height " + str(OnChipIB_Height) + "\n"
head_define_str += "#define TRow_max " + str(TRow_max) + "\n"
head_define_str += "#define TCol_max " + str(TCol_max) + "\n"

print(head_define_str)

f_head = open("yolov2_acc_gen_template.h","rb")
f_new = open("yolov2_acc_sim_gen.h","wb")

find_str = "#DEFINE_HEADER#".encode()
replace_str = head_define_str.encode()

for line in f_head:
    if find_str in line:
        line = line.replace(find_str,replace_str)
    f_new.write(line)
    
f_head.close()
f_new.close()

f_head = open("yolov2_acc_test_template.h","rb")
f_new = open("yolov2_acc_sim_test.h","wb")

find_str = "#DEFINE_HEADER#".encode()
replace_str = head_define_str.encode()

for line in f_head:
    if find_str in line:
        line = line.replace(find_str,replace_str)
    f_new.write(line)
    
f_head.close()
f_new.close()

f_head = open("../../hls/src_float32/cnn_template.h","rb")
f_new = open("../../hls/src_float32/cnn.h","wb")

find_str = "#DEFINE_HEADER#".encode()
replace_str = head_define_str.encode()

for line in f_head:
    if find_str in line:
        line = line.replace(find_str,replace_str)
    f_new.write(line)
    
f_head.close()
f_new.close()
