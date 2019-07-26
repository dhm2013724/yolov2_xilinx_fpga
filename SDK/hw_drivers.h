/*
 * hw_drivers.h
 *
 *  Created on: 2017Äê7ÔÂ26ÈÕ
 *      Author: shen_lee
 */
#ifndef __HW_DRIVERS_H_
#define __HW_DRIVERS_H_

#include <stdlib.h>
#include <stdint.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <stdio.h>

/*
 *  hw_maamp_init
 *  input: void
 * output: fd of memeory
 */
int hw_mmap_init();
/*
 *  hw_mmap
 *  input:  physics address, map length
 * 	return: virtual address
 */
void *hw_mmap(uint32_t phy_addr, uint32_t map_len);

/*
 *  hw_unmap
 *   input: start address, map_length
 *  output: none
 */
void hw_unmap(void * start, uint32_t length);

#endif
