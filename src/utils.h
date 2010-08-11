/*
 *  catnet : categorical Bayesian network inference
 *  Copyright (C) 2009--2010  Nikolay Balov
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.gnu.org/licenses/gpl-2.0.html
 */

/*
 * utils.h
 *
 *  Created on: Nov 16, 2009
 *      Author: nbalov
 */

#ifndef UTILS_H_
#define UTILS_H_

#include "thread.h"

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

template<class t_elem>
int _gen_permutation(t_elem *psample, int nsample) {
	int i, j;
	if(nsample < 1 || !psample)
		return -1;
	double *paux = (double*)malloc(nsample*sizeof(double));
	for(i = 0; i < nsample; i++)
		paux[i] = (double)rand() / (double)RAND_MAX;
	for(j = 0; j < nsample; j++) {
		psample[j] = 0;
		for(i = 0; i < nsample; i++) {
			if(paux[j] >= paux[i])
				psample[j]++;
		}
	}
	free(paux);
//printf("permutation ");
//for(i=0;i<nsample;i++)
//printf("%d ", psample[i]);
//printf("\n");
	return 0;
}

template<class t_elem>
int _gen_binomial(int size, t_elem prob) {
	int i, r;
	if(size < 1 || prob <= 0)
		return 0;
	r = 0;
	for(i = 0; i < size; i++) {
		if((t_elem)rand() / (t_elem)RAND_MAX <= prob)
			r++;
	}
	return r;
}


template<class t_prob>
t_prob _gamma_upper_bound(t_prob p, int m, t_prob delta) {
	t_prob fr1, fr2, fp, fres;
	if(p <= 0 || p >= 1 || m < 1)
		return 0;
	if(delta <= 0)
		return 1;
	fp = 1 + log((double)p);
	if(fp < 0) fp = -fp;
	fp /= delta;

	fr1 = fp + sqrt((double)(fp*fp + 2/(p*delta)));
	fr2 = 1/p + fp + sqrt((double)((1/p-fp)*(1/p-fp) + 2/(p*delta)));
	fres = fr1*fr1 + fr2*fr2;

	return(p*(1-p)*fres/(8*m));
}

#endif /* UTILS_H_ */

