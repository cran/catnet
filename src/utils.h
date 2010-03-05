/*
 * utils.h
 *
 *  Created on: Nov 16, 2009
 *      Author: Nikolay Balov
 */

#ifndef UTILS_H_
#define UTILS_H_

#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <float.h>
#include <time.h>

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

#ifndef CATNET_PI
#define CATNET_PI	(double)3.14159265358979323846264338327950288
#endif
#ifndef CATNET_PI2
#define CATNET_PI2	(2*(double)CATNET_PI)
#endif

void * CATNET_MALLOC(size_t nsize);
void CATNET_FREE(void *pMem);

template<class t_elem>
void _quick_sort(t_elem *plist, int nlist){
	t_elem pivot;
	int j, nless, ngreater;
	if(nlist <= 1)
		return;
	t_elem *paux = (t_elem*)malloc(nlist*sizeof(t_elem));
	nless = 0;
	ngreater = nlist-1;
	pivot = plist[0];
	for(j = 1; j < nlist; j++) {
		if(plist[j] <= pivot) {
			paux[nless] = plist[j];
			nless++;
		}
		else {
			paux[ngreater] = plist[j];
			ngreater--;
		}
	}
	_quick_sort<t_elem>(paux, nless);
	_quick_sort<t_elem>(paux + ngreater + 1, nlist - ngreater - 1);
	paux[nless] = pivot;
	memcpy(plist, paux, nlist*sizeof(t_elem));
	free(paux);
	return;
}


template<class t_elem>
t_elem _gen_std_normal_var() {
	/* ISO C pseudo random generator */
	/* include stdlib.h and math.h */
	t_elem u, v;
	u = (t_elem)rand() / (t_elem)RAND_MAX;
	v = (t_elem)rand() / (t_elem)RAND_MAX;
	return(sqrt(-2*log(u)) * cos(CATNET_PI2*v));
}

#endif /* UTILS_H_ */

