# import numpy as np
import sys,os

HW_S = 2
K = 3
Tn = 2
Tm = 60
Tr = 26
Tc = 26
MAX_BETA_LENGTH = 1024
#for int16
INTERWIDTH = 19
OnChipIB_Width  = (Tc-1)*HW_S+K
OnChipIB_Height = (Tr-1)*HW_S+K


head_define_str = ""
head_define_str += "#define HW_S " + str(HW_S) + "\n"
head_define_str += "#define K " + str(K) + "\n"
head_define_str += "#define MAX_BETA_LENGTH " + str(MAX_BETA_LENGTH) + "\n"
head_define_str += "#define Tn " + str(Tn) + "\n"
head_define_str += "#define Tm " + str(Tm) + "\n"
head_define_str += "#define Tr " + str(Tr) + "\n"
head_define_str += "#define Tc " + str(Tc) + "\n"
head_define_str += "#define OnChipIB_Width " + str(OnChipIB_Width) + "\n"
head_define_str += "#define OnChipIB_Height " + str(OnChipIB_Height) + "\n"
head_define_str += "#define INTERWIDTH " + str(INTERWIDTH) + "\n"

print(head_define_str)

f_head = open("acc_f32_t.h","rb")
f_new = open("acc_f32.h","wb")

find_str = "#DEFINE_HEADER#".encode()
replace_str = head_define_str.encode()

for line in f_head:
    if find_str in line:
        line = line.replace(find_str,replace_str)
    f_new.write(line)
    
f_head.close()
f_new.close()

f_head = open("acc_i16_t.h","rb")
f_new = open("acc_i16.h","wb")

find_str = "#DEFINE_HEADER#".encode()
replace_str = head_define_str.encode()

for line in f_head:
    if find_str in line:
        line = line.replace(find_str,replace_str)
    f_new.write(line)
    
f_head.close()
f_new.close()

f_head = open("../../hls/src_int16/cnn_t.h","rb")
f_new = open("../../hls/src_int16/cnn.h","wb")

find_str = "#DEFINE_HEADER#".encode()
replace_str = head_define_str.encode()

for line in f_head:
    if find_str in line:
        line = line.replace(find_str,replace_str)
    f_new.write(line)
    
f_head.close()
f_new.close()
