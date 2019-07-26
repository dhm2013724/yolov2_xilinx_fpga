/*
 * hw_drivers.c
 *
 *  Created on: 2017Äê7ÔÂ26ÈÕ
 *      Author: shen_lee
 */

#include "hw_drivers.h"


int fd = 0;
/*
 *  hw_maamp_init
 *  input: void
 * output: fd of memeory
 */
int hw_mmap_init()
{
	fd = open("/dev/mem", O_RDWR);

	if(fd == -1) {
		printf("Open memory failed\r\n");
		return -1;
	}

	else return fd;
}
/*
 *  hw_mmap
 *  input:  physics address, map length
 * 	return: virtual address
 */
void * hw_mmap(uint32_t phy_addr, uint32_t map_len)
{

	void *virtual_addr = 0;
	if(map_len < 0x1000) map_len = 0x1000;
	virtual_addr = mmap(0, map_len, PROT_READ|PROT_WRITE, MAP_SHARED, fd, phy_addr);
	if(virtual_addr == NULL){
		printf("hw_mmap failed\r\n");
	}
	return virtual_addr;
}

/*
 *  hw_unmap
 *  input: start address, map_length
 */
void hw_unmap(void * start, uint32_t length)
{
	if(length < 0x1000) length = 0x1000;
	munmap(start, length);
}
