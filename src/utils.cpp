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
	void *pMem;
	pMem = malloc(nsize);
	if(!pMem) {
		error("Insufficient memory");
		return 0;
	}
	return pMem;
	// memmonitor
	g_memcounter += nsize;
	pMem = malloc(sizeof(int) + nsize);
	if(!pMem) {
		error("Insufficient memory");
		return 0;
	}
	*(int*)pMem = nsize;
	pMem = (void*)((int*)pMem + 1);
	//char str[128];
	//sprintf(str, "+%d    %d        %p\n", (int)nsize, (int)g_memcounter, pMem);
	//fprintf(g_hf,str);
	//printf(str);
	return pMem;
}

void CATNET_FREE(void *pMem) {
	if(!pMem)	
		return;
	free(pMem);
	return;
	// memmonitor
	pMem = (void*)((int*)pMem-1);
	size_t nsize = *((int*)pMem);
	g_memcounter -= nsize;
	//char str[128];
	//sprintf(str, "-%d    %d\n", (int)nsize, (int)g_memcounter);
	//printf(str);
	free(pMem);
}

void CATNET_MEM_ERR() {
	// generate R-errorx
	error("Insufficient memory");
}
