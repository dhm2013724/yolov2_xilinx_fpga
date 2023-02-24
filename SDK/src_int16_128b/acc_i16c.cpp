
#include "acc_i16c.h"

void copy_mem2dev(uint8_t *orig,uint32_t byte_num, unsigned long in_buffer)
{
	int fd = open("/dev/mem", O_RDWR);
	unsigned char *virtual_addr;
	uint32_t RequestByteNum;// must page
	if(byte_num%(HPAGESIZE)==0)
		RequestByteNum = byte_num;
	else
	{
		RequestByteNum = ceil(byte_num/(HPAGESIZE*1.0))*(HPAGESIZE);
	}
	virtual_addr = (unsigned char *)mmap(NULL, RequestByteNum, PROT_READ | PROT_WRITE, MAP_SHARED, fd, (off_t)in_buffer);
	if(virtual_addr == MAP_FAILED)
	{
		perror("Virtual_addr_in mappong for absolute memory access failed!\n");
		return;
	}
	memcpy(virtual_addr,orig,byte_num);

	munmap((void *)virtual_addr, byte_num);
	close(fd);
}

void copy_dev2mem(uint8_t *dst,uint32_t byte_num, unsigned long in_buffer)
{
	int fd = open("/dev/mem", O_RDWR);
	unsigned char *virtual_addr;
	uint32_t RequestByteNum;// must page
	if(byte_num%(HPAGESIZE)==0)
		RequestByteNum = byte_num;
	else
	{
		RequestByteNum = ceil(byte_num/(HPAGESIZE*1.0))*(HPAGESIZE);
	}
		virtual_addr = (unsigned char *)mmap(NULL, RequestByteNum, PROT_READ | PROT_WRITE, MAP_SHARED, fd, (off_t)in_buffer);
	if(virtual_addr == MAP_FAILED)
	{
		perror("Virtual_addr_in mappong for absolute memory access failed!\n");
		return;
	}
	printf("copy start-----byte_num=%d\n",byte_num);
	memcpy((uint8_t *)dst,virtual_addr,byte_num);
	printf("copy ok!\n");

	munmap((void *)virtual_addr, byte_num);
	close(fd);
}

int copy_file2mem(char *bin_file,uint32_t byte_num,unsigned long in_buffer)
{
	unsigned char *buffer = (unsigned char *)malloc(HPAGESIZE);
	if(buffer==NULL){
		printf("cannot malloc buffer %d byte\n", HPAGESIZE);
		return -1;
	}
	printf("Total Byte Num = %d\n Address 0x%X\n", byte_num, in_buffer);
	FILE *fp;
	if( (fp = fopen(bin_file, "rb")) == NULL)fprintf(stderr,"CANNOT OPEN bin_file\n");
	int rd_num;
	unsigned long offset = 0;
	while(rd_num = fread(buffer, sizeof(unsigned char), HPAGESIZE, fp))
	{
		if(rd_num < HPAGESIZE)
			rd_num = HPAGESIZE;
		copy_mem2dev(buffer,rd_num, in_buffer+offset);
//		printf("rd_num=%d, offset=%d\n", rd_num, offset);
		offset += rd_num;
	}
	printf("copy_file2mem offset=%d\n",offset);
	fclose(fp);

	free(buffer);


	return 0;
}

int copy_mem2file(char *bin_file,uint32_t byte_num,unsigned long in_buffer)
{
	void *buffer = malloc(HPAGESIZE);
	if(buffer==NULL){
		printf("cannot malloc buffer %d byte\n", HPAGESIZE);
		return -1;
	}

	FILE *fp;
	if( (fp = fopen(bin_file, "wb")) == NULL)fprintf(stderr,"CANNOT OPEN bin_file\n");

	int x = byte_num;
	int addbyte;
	unsigned long offset = 0;
	while(addbyte=((x<HPAGESIZE)?x:(HPAGESIZE)))
	{
		copy_dev2mem((uint8_t *)buffer,addbyte, in_buffer+offset);
		fwrite(buffer , sizeof(unsigned char), addbyte, fp);
		x -= addbyte;
		offset += addbyte;
	}
	printf("copy_mem2file offset=%d\n",offset);


	fclose(fp);

	free(buffer);

	return 0;
}

double what_time_is_it_now()
{
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

int FPGA_Acc(unsigned int ifm_addr, unsigned int  ofm_addr, unsigned int weight_offset, unsigned int bias_offset, uint32_t k_s_pad_ltype, uint32_t iofm_num, uint32_t ifm_w_h, uint32_t ofm_w_h,
	uint32_t TRTC, uint32_t TMTN, int32_t NToy, int32_t NTox, int32_t NTof, int32_t NTcomb, int32_t NTif, uint8_t lmode, int32_t NTcomb_l, int16_t pad_val, int16_t div_kk,
	uint32_t IHW, uint32_t OHW, uint32_t KK_INumxKK, uint32_t en_bits, uint32_t weightQ, uint32_t biasQ, uint32_t ifmQ, uint32_t ofmQ, uint32_t avgQ, uint32_t interQ)//enable_bits[2:0]={IsReLU, LoadBias, IsNotConv}
{
	unsigned int ap_idle;
	unsigned int ap_done;

	unsigned long int PhysicalAddress = YOLO2_BASEADDR;
	int map_len = 0x120;
	int fd = open("/dev/mem", O_RDWR);

	unsigned char *xbase_address;
	xbase_address = (unsigned char *)mmap(NULL, map_len, PROT_READ | PROT_WRITE, MAP_SHARED, fd, (off_t)PhysicalAddress);
	if(xbase_address == MAP_FAILED)
	{
		perror("1:Init Mapping memory for absolute memory access failed.\n");
		return -1;
	}

	while(1)
	{
		ap_idle = ((ReadReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_AP_CTRL) >> 2) && 0x1);
		if(ap_idle)
			break;
	}

	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_IFM_V_DATA,  ifm_addr);
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_OFM_V_DATA, ofm_addr);
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_WEIGHT_V_DATA, WEIGHT_BASE + weight_offset*2);
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_BIAS_V_DATA, BETA_BASE + bias_offset*2);

	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_K_S_PAD_LTYPE_DATA,  k_s_pad_ltype);
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_IOFM_NUM_DATA,  iofm_num);
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_IFM_W_H_DATA,  ifm_w_h);
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_OFM_W_H_DATA,  ofm_w_h);
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_TRTC_DATA,  TRTC);
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_TMTN_DATA,  TMTN);
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_NTOY_DATA,  NToy);
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_NTOX_DATA,  NTox);
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_NTOF_DATA,  NTof);
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_NTCOMB_DATA,  NTcomb);
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_NTIF_DATA,  NTif);
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_LMODE_DATA,  lmode);
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_NTCOMB_L_DATA,  NTcomb_l);
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_PAD_VAL_DATA,  pad_val);
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_DIV_KK_DATA,  div_kk);
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_IHW_DATA,  IHW);
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_OHW_DATA,  OHW);
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_KK_INUMXKK_DATA,  KK_INumxKK);
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_EN_BITS_DATA,  en_bits);
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_WEIGHTQ_DATA,  weightQ);
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_BIASQ_DATA,  biasQ);
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_IFMQ_DATA,  ifmQ);
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_OFMQ_DATA,  ofmQ);
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_AVGQ_DATA,  avgQ);
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_INTERQ_DATA,  interQ);

//	double time1,time2;
//	time1 = what_time_is_it_now();
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_GIE, 0x0);
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_AP_CTRL, 0x1);//Start
	while(1)
	{
		ap_done = ((ReadReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_AP_CTRL) >> 1) && 0x1);
		if(ap_done)
			break;
	}
//	time2 = what_time_is_it_now();
//	printf("START TO DONE in %f seconds.\n",time2 - time1);

	munmap((void *)xbase_address, map_len);
	close(fd);

	return 0;
}
