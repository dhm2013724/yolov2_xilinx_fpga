
#include "cnn.h"

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
	memcpy((uint8_t *)dst,virtual_addr,byte_num);

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
		printf("rd_num=%d, offset=%d\n", rd_num, offset);
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

int FPGA_Acc(uint64_t In_Address, uint64_t Out_Address, uint64_t Weight_offset, uint64_t Beta_offset, uint32_t k_s_pad_ltype, uint32_t iofm_num, uint32_t ifm_w_h, uint32_t ofm_w_h,
	uint32_t TRTC, uint32_t TMTN, int32_t OFM_num_bound, int32_t mLoopsxTM, int32_t mLoops_a1xTM, int16_t pad_val, uint32_t TRowTCol,
	uint32_t IHW, uint32_t OHW, uint32_t KK_INumxKK, uint32_t en_bits, int32_t WeightQ, int32_t BetaQ, int32_t InputQ, int32_t OutputQ)
{

	unsigned int ap_idle;
	unsigned int ap_done;

	uint64_t PhysicalAddress = ACC_BASEADDR;
	int map_len = 0x100;
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
		ap_idle = ((ReadReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_AP_CTRL) >> 2) & 0x1);
		if(ap_idle)
			break;
	}

	union {
		int32_t i32;
		float f32;
		uint32_t u32;
	} tmp_32;

	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_IFM_DATA,  In_Address);
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_OFM_DATA, Out_Address);
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_WEIGHT_DATA,   WEIGHT_BASEADDR + Weight_offset*4);
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_BIAS_DATA,     BETA_BASEADDR + Beta_offset*4);

	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_K_S_PAD_LTYPE_DATA, k_s_pad_ltype);
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_IOFM_NUM_DATA, iofm_num);
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_IFM_W_H_DATA, ifm_w_h);
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_OFM_W_H_DATA, ofm_w_h);
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_TRTC_DATA, TRTC);
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_TMTN_DATA, TMTN);

	tmp_32.i32 = OFM_num_bound;
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_OFM_NUM_BOUND_DATA, tmp_32.u32);
	tmp_32.i32 = mLoopsxTM;
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_MLOOPSXTM_DATA, tmp_32.u32);
	tmp_32.i32 = mLoops_a1xTM;
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_MLOOPS_A1XTM_DATA, tmp_32.u32);

	tmp_32.i32 = pad_val;
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_PAD_VAL_DATA, tmp_32.u32);

	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_TROWTCOL_DATA, TRowTCol);
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_IHW_DATA, IHW);
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_OHW_DATA, OHW);
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_KK_INUMXKK_DATA, KK_INumxKK);
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_EN_BITS_DATA, en_bits);

	tmp_32.i32 = WeightQ;
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_WEIGHTQ_DATA, tmp_32.u32);
	tmp_32.i32 = BetaQ;
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_BETAQ_DATA, tmp_32.u32);
	tmp_32.i32 = InputQ;
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_INPUTQ_DATA, tmp_32.u32);
	tmp_32.i32 = OutputQ;
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_OUTPUTQ_DATA, tmp_32.u32);

	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_GIE, 0x0);
	WriteReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_AP_CTRL, 0x1);//Start
	while(1)
	{
		ap_done = ((ReadReg(xbase_address, XFPGA_ACC_CTRL_BUS_ADDR_AP_CTRL) >> 1) & 0x1);
		if(ap_done)
			break;
	}


	munmap((void *)xbase_address, map_len);
	close(fd);

	return 0;


}
