/*
 * utils.c
 *
 *  Created on: Nov 16, 2009
 *      Author: Nikolay Balov
 */

#include "utils.h"

extern size_t g_memcounter;

void * CATNET_MALLOC(size_t nsize) { 
	if(nsize <= 0)
		return 0;
	return malloc(nsize);
	/*g_memcounter += nsize;
	char str[128];
	sprintf(str, "+%d    %d\n", (int)nsize, (int)g_memcounter);
	printf(str);
	void *pMem = malloc(sizeof(int) + nsize);
	if(pMem) {
		*(int*)pMem = nsize;
		pMem = (void*)((int*)pMem + 1);
	}
	return pMem;*/
}

void CATNET_FREE(void *pMem) {
	if(!pMem)	
		return;
	free(pMem);
	return;
	/*pMem = (void*)((int*)pMem-1);
	size_t nsize = *((int*)pMem);
	g_memcounter -= nsize;
	char str[128];
	sprintf(str, "-%d    %d\n", (int)nsize, (int)g_memcounter);
	printf(str);
	free(pMem);*/
}

