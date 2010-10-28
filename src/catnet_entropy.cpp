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
 * rcatnet_kl.cpp
 *
 *  Created on: Sep 21, 2009
 *      Author: nbalov
 */

#include "utils.h"
#include "rcatnet.h"

extern "C" {

extern size_t g_memcounter;

SEXP catnetEntropyPairwise(SEXP rSamples, SEXP rPerturbations) {

	int *pSamples, *pPerturbations;
	int *psubSamples, numsubsamples, *psamples1, numsamples1;
	int *pNodeNumCats, **pNodeCats, mincat, maxcat, maxCategories, *pprobs;
	int numsamples, numnodes, i, j, k, d, nnode1, nnode2;
	double floglik, faux, fsum, *pvec, *klmat;
	SEXP dim, rvec = R_NilValue;

	if(!isMatrix(rSamples))
		error("Data should be a matrix");
	if(!isNull(rPerturbations) && !isMatrix(rPerturbations))
		error("Perturbations should be a matrix");

	PROTECT(rSamples = AS_INTEGER(rSamples));
	pSamples = INTEGER(rSamples);

	dim = GET_DIM(rSamples);
	numnodes = INTEGER(dim)[0];
	numsamples = INTEGER(dim)[1];

	// pSamples are assumed positive indices
	for(j = 0; j < numnodes*numsamples; j++) {
		pSamples[j]--;
	}

	// categoies
	pNodeNumCats = (int*)CATNET_MALLOC(numnodes*sizeof(int));
	pNodeCats = (int**)CATNET_MALLOC(numnodes*sizeof(int*));
	memset(pNodeCats, 0, numnodes*sizeof(int*));
	memset(pNodeNumCats, 0, numnodes*sizeof(int));

	maxCategories = 1;
	for(i = 0; i < numnodes; i++) {
		mincat = INT_MAX;
		maxcat = -INT_MAX;
		for(j = 0; j < numsamples; j++) {
			if(pSamples[j*numnodes + i] < mincat)
				mincat = pSamples[j*numnodes + i];
			if(pSamples[j*numnodes + i] > maxcat)
				maxcat = pSamples[j*numnodes + i];
		}
		pNodeNumCats[i] = maxcat - mincat + 1;
		pNodeCats[i] = (int*)CATNET_MALLOC(pNodeNumCats[i]*sizeof(int));
		for(j = 0; j < pNodeNumCats[i]; j++)
			pNodeCats[i][j] = mincat + j;
	}
	for(i = 0; i < numnodes; i++) {
		/* order pNodeNumCats[i] */
		for(j = 0; j < pNodeNumCats[i]; j++) {
			for(k = j + 1; k < pNodeNumCats[i]; k++) {
				if(pNodeCats[i][j] > pNodeCats[i][k]) {
					d = pNodeCats[i][j]; 
					pNodeCats[i][j] = pNodeCats[i][k];
					pNodeCats[i][k] = d;
				}
			}
		} 
		for(j = 0; j < numsamples; j++) {
			for(d = 0; d < pNodeNumCats[i]; d++)
				if(pNodeCats[i][d] == pSamples[j*numnodes + i])
					break;
			pSamples[j*numnodes + i] = d;
		}
		if(maxCategories < pNodeNumCats[i])
			maxCategories = pNodeNumCats[i];
	}
	//printf("maxCategories = %d\n", maxCategories);
	pprobs = (int*)CATNET_MALLOC(maxCategories*maxCategories*sizeof(int));

	klmat = (double*)CATNET_MALLOC(numnodes*numnodes*sizeof(double));
	memset(klmat, 0, numnodes*numnodes*sizeof(double));

	psubSamples = 0;
	pPerturbations = 0;
	if(!isNull(rPerturbations)) {
		PROTECT(rPerturbations = AS_INTEGER(rPerturbations));
		pPerturbations = INTEGER(rPerturbations);
		psubSamples = (int*)CATNET_MALLOC(numnodes*numsamples*sizeof(int));
	}

	for(nnode1 = 0; nnode1 < numnodes; nnode1++) {
		psamples1 = pSamples;
		numsamples1 = numsamples;
		if(pPerturbations) {
			numsubsamples = 0;
			for(j = 0; j < numsamples; j++) {
				if(!pPerturbations[j * numnodes + nnode1]) {
					memcpy(psubSamples + numsubsamples*numnodes, pSamples + j*numnodes, numnodes*sizeof(int));
					numsubsamples++;
				}
				psamples1 = psubSamples;
				numsamples1 = numsubsamples;
			}
		}
		//printf("nnode = %d (%d)\n", nnode1, numsamples1);
		for(nnode2 = 0; nnode2 < numnodes; nnode2++) {
	
			memset(pprobs, 0, maxCategories*maxCategories*sizeof(int));
		
			if(nnode2 == nnode1) {
				for(j = 0; j < numsamples1; j++) 
					pprobs[psamples1[j*numnodes + nnode1]]++;
				floglik = 0;
				fsum  = 0;
				for(j = 0; j < pNodeNumCats[nnode1]; j++) {
					fsum += pprobs[j];
					if(pprobs[j] > 0)
						floglik += pprobs[j]*(double)log((double)pprobs[j]);
				}
				if(fsum > 0) {
					floglik -= fsum*(double)log((double)fsum);
					floglik /= fsum;
					//floglik *= numsamples;
				}
				klmat[nnode2*numnodes + nnode1] = -floglik;
				continue;
			}

			// estimate logP(nnode1|nnode2)
			for(j = 0; j < numsamples1; j++) 
				pprobs[maxCategories*psamples1[j*numnodes + nnode2] + psamples1[j*numnodes + nnode1]]++;

			floglik = 0;
			fsum  = 0;
			for(i = 0; i < pNodeNumCats[nnode2]; i++) {
				faux = 0;
				for(j = 0; j < pNodeNumCats[nnode1]; j++) {
					faux += pprobs[maxCategories*i+j];
					if(pprobs[maxCategories*i+j] > 0)
						floglik += pprobs[maxCategories*i+j]*(double)log((double)pprobs[maxCategories*i+j]);
				}
				fsum += faux;
				if(faux > 0) {
					floglik -= faux*(double)log((double)faux);
				}
			}
			if(fsum > 0) {
				floglik /= fsum;
				//floglik *= numsamples;
			}
			klmat[nnode2*numnodes + nnode1] = -floglik;
		}
	}

	UNPROTECT(1); //rSamples
	if(!isNull(rPerturbations))
		UNPROTECT(1); //rPerturbations

	if(psubSamples)
		CATNET_FREE(psubSamples);

	if(pprobs)
		CATNET_FREE(pprobs);

	if(pNodeCats) {
		for(i = 0; i < numnodes; i++) 
			if(pNodeCats[i])
				CATNET_FREE(pNodeCats[i]);
		CATNET_FREE(pNodeCats);
	}

	if(pNodeNumCats) 
		CATNET_FREE(pNodeNumCats);

	if(klmat) {
		PROTECT(rvec = NEW_NUMERIC(numnodes*numnodes));
		pvec = NUMERIC_POINTER(rvec);
		memcpy(pvec, klmat, numnodes*numnodes*sizeof(double));
		UNPROTECT(1);
	}

	CATNET_FREE(klmat);

	//char str[128];
	//sprintf(str, "Mem Balance  %d\n", (int)g_memcounter);
	//printf(str);

	return rvec;

}


SEXP catnetEntropyOrder(SEXP rSamples, SEXP rPerturbations) {

	int *pSamples, *pPerturbations;
	int *psubSamples, numsubsamples, *psamples1, numsamples1;
	int *pNodeNumCats, **pNodeCats, mincat, maxcat, maxCategories, *pprobs;
	int numsamples, numnodes, i, j, k, d, nnode1, nnode2;
	int *porder, *pvec;
	double floglik, faux, fsum, *klmat, *pnoderanks;
	SEXP dim, rvec = R_NilValue;

	if(!isMatrix(rSamples))
		error("Data should be a matrix");
	if(!isNull(rPerturbations) && !isMatrix(rPerturbations))
		error("Perturbations should be a matrix");

	PROTECT(rSamples = AS_INTEGER(rSamples));
	pSamples = INTEGER(rSamples);

	dim = GET_DIM(rSamples);
	numnodes = INTEGER(dim)[0];
	numsamples = INTEGER(dim)[1];

	// pSamples are assumed positive indices
	for(j = 0; j < numnodes*numsamples; j++) {
		pSamples[j]--;
	}

	porder = (int*)CATNET_MALLOC(numnodes*sizeof(int));
	pnoderanks = (double*)CATNET_MALLOC(2*numnodes*sizeof(double));

	// categoies
	pNodeNumCats = (int*)CATNET_MALLOC(numnodes*sizeof(int));
	pNodeCats = (int**)CATNET_MALLOC(numnodes*sizeof(int*));
	memset(pNodeCats, 0, numnodes*sizeof(int*));
	memset(pNodeNumCats, 0, numnodes*sizeof(int));

	maxCategories = 1;
	for(i = 0; i < numnodes; i++) {
		mincat = INT_MAX;
		maxcat = -INT_MAX;
		for(j = 0; j < numsamples; j++) {
			if(pSamples[j*numnodes + i] < mincat)
				mincat = pSamples[j*numnodes + i];
			if(pSamples[j*numnodes + i] > maxcat)
				maxcat = pSamples[j*numnodes + i];
		}
		pNodeNumCats[i] = maxcat - mincat + 1;
		pNodeCats[i] = (int*)CATNET_MALLOC(pNodeNumCats[i]*sizeof(int));
		for(j = 0; j < pNodeNumCats[i]; j++)
			pNodeCats[i][j] = mincat + j;
	}
	for(i = 0; i < numnodes; i++) {
		/* order pNodeNumCats[i] */
		for(j = 0; j < pNodeNumCats[i]; j++) {
			for(k = j + 1; k < pNodeNumCats[i]; k++) {
				if(pNodeCats[i][j] > pNodeCats[i][k]) {
					d = pNodeCats[i][j]; 
					pNodeCats[i][j] = pNodeCats[i][k];
					pNodeCats[i][k] = d;
				}
			}
		} 
		for(j = 0; j < numsamples; j++) {
			for(d = 0; d < pNodeNumCats[i]; d++)
				if(pNodeCats[i][d] == pSamples[j*numnodes + i])
					break;
			pSamples[j*numnodes + i] = d;
		}
		if(maxCategories < pNodeNumCats[i])
			maxCategories = pNodeNumCats[i];
	}
	//printf("maxCategories = %d\n", maxCategories);
	pprobs = (int*)CATNET_MALLOC(maxCategories*maxCategories*sizeof(int));

	klmat = (double*)CATNET_MALLOC(numnodes*numnodes*sizeof(double));
	memset(klmat, 0, numnodes*numnodes*sizeof(double));

	psubSamples = 0;
	pPerturbations = 0;
	if(!isNull(rPerturbations)) {
		PROTECT(rPerturbations = AS_INTEGER(rPerturbations));
		pPerturbations = INTEGER(rPerturbations);
		psubSamples = (int*)CATNET_MALLOC(numnodes*numsamples*sizeof(int));
	}

	for(nnode1 = 0; nnode1 < numnodes; nnode1++) {
		psamples1 = pSamples;
		numsamples1 = numsamples;
		if(pPerturbations) {
			numsubsamples = 0;
			for(j = 0; j < numsamples; j++) {
				if(!pPerturbations[j * numnodes + nnode1]) {
					memcpy(psubSamples + numsubsamples*numnodes, pSamples + j*numnodes, numnodes*sizeof(int));
					numsubsamples++;
				}
				psamples1 = psubSamples;
				numsamples1 = numsubsamples;
			}
		}
		//printf("nnode = %d (%d)\n", nnode1, numsamples1);
		for(nnode2 = 0; nnode2 < numnodes; nnode2++) {
	
			memset(pprobs, 0, maxCategories*maxCategories*sizeof(int));
		
			if(nnode2 == nnode1) {
				for(j = 0; j < numsamples1; j++) 
					pprobs[psamples1[j*numnodes + nnode1]]++;
				floglik = 0;
				fsum  = 0;
				for(j = 0; j < pNodeNumCats[nnode1]; j++) {
					fsum += pprobs[j];
					if(pprobs[j] > 0)
						floglik += pprobs[j]*(double)log((double)pprobs[j]);
				}
				if(fsum > 0) {
					floglik -= fsum*(double)log((double)fsum);
					floglik /= fsum;
					//floglik *= numsamples;
				}
				klmat[nnode2*numnodes + nnode1] = -floglik;
				continue;
			}

			// estimate logP(nnode1|nnode2)
			for(j = 0; j < numsamples1; j++) 
				pprobs[maxCategories*psamples1[j*numnodes + nnode2] + psamples1[j*numnodes + nnode1]]++;

			floglik = 0;
			fsum  = 0;
			for(i = 0; i < pNodeNumCats[nnode2]; i++) {
				faux = 0;
				for(j = 0; j < pNodeNumCats[nnode1]; j++) {
					faux += pprobs[maxCategories*i+j];
					if(pprobs[maxCategories*i+j] > 0)
						floglik += pprobs[maxCategories*i+j]*(double)log((double)pprobs[maxCategories*i+j]);
				}
				fsum += faux;
				if(faux > 0) {
					floglik -= faux*(double)log((double)faux);
				}
			}
			if(fsum > 0) {
				floglik /= fsum;
				//floglik *= numsamples;
			}
			klmat[nnode2*numnodes + nnode1] = -floglik;
		}
	}

	UNPROTECT(1); //rSamples
	if(!isNull(rPerturbations))
		UNPROTECT(1); //rPerturbations

	if(psubSamples)
		CATNET_FREE(psubSamples);

	if(pprobs)
		CATNET_FREE(pprobs);

	if(pNodeCats) {
		for(i = 0; i < numnodes; i++) 
			if(pNodeCats[i])
				CATNET_FREE(pNodeCats[i]);
		CATNET_FREE(pNodeCats);
	}

	if(pNodeNumCats) 
		CATNET_FREE(pNodeNumCats);

	for(nnode1 = 0; nnode1 < numnodes; nnode1++) {
		faux = klmat[nnode1*numnodes + nnode1];
		floglik = 0;
		for(nnode2 = 0; nnode2 < numnodes; nnode2++) {
			klmat[nnode2*numnodes + nnode1] = faux - klmat[nnode2*numnodes + nnode1];
			floglik += klmat[nnode2*numnodes + nnode1];
		}
		pnoderanks[nnode1] = floglik;
	}

//printf("pnoderanks = ");
//for(i = 0; i < numnodes; i++)
//printf(" %f", pnoderanks[i]);
//printf("\n");
	_order<double>(pnoderanks, numnodes, porder, 1);
//printf("porder = ");
//for(i = 0; i < numnodes; i++)
//printf(" %d", porder[i]+1);
//printf("\n");

	/*fsum = 0;
	pnoderanks[porder[0]] = 0;
	for(nnode1 = 1; nnode1 < numnodes; nnode1++) {
		faux = 0;
		for(i = 0; i < nnode1; i++)
			if(faux < klmat[porder[i]*numnodes + porder[nnode1]])
				faux = klmat[porder[i]*numnodes + porder[nnode1]];
		pnoderanks[porder[nnode1]] = faux;
printf("pnoderanks[%d] = %f\n", porder[nnode1], pnoderanks[porder[nnode1]]);
		fsum += faux;
	}

	for(nnode1 = 0; nnode1 < numnodes; nnode1++) {
		for(nnode2 = nnode1 + 1; nnode2 < numnodes; nnode2++) {
			memcpy(pnoderanks + numnodes, pnoderanks, numnodes*sizeof(double));

			for(d = nnode1; d <= nnode2; d++) {
				faux = 0;
				for(i = 0; i < d; i++) {
					if(i == nnode1 && faux < klmat[porder[nnode2]*numnodes + porder[d]])
						faux = klmat[porder[nnode2]*numnodes + porder[d]];
					else if(i == nnode2 && faux < klmat[porder[nnode1]*numnodes + porder[d]])
						faux = klmat[porder[nnode1]*numnodes + porder[d]];
					else if(faux < klmat[porder[i]*numnodes + porder[d]])
						faux = klmat[porder[i]*numnodes + porder[d]];
				}
				if(d == nnode1)
					pnoderanks[numnodes+porder[nnode2]] = faux;
				else if(d == nnode2)
					pnoderanks[numnodes+porder[nnode1]] = faux;
				else
					pnoderanks[numnodes+porder[d]] = faux;
			}
			faux = 0;
			for(i = 0; i < numnodes; i++)
				faux += pnoderanks[numnodes+i];

			if(faux > fsum) {
printf("change %d <-> %d    %f (%f)\n", porder[nnode1]+1, porder[nnode2]+1, fsum, faux);
				fsum = faux;
				i = porder[nnode1];
				porder[nnode1] = porder[nnode2];
				porder[nnode2] = i;
				memcpy(pnoderanks, pnoderanks + numnodes, numnodes*sizeof(double));
			}
		}
	}*/

	CATNET_FREE(klmat);

	for(i = 0; i < numnodes; i++)
		porder[i]++;

	PROTECT(rvec = NEW_INTEGER(numnodes));
	pvec = INTEGER_POINTER(rvec);
	memcpy(pvec, porder, numnodes*sizeof(int));
	UNPROTECT(1);

	CATNET_FREE(porder);
	CATNET_FREE(pnoderanks);

	//char str[128];
	//sprintf(str, "Mem Balance  %d\n", (int)g_memcounter);
	//printf(str);

	return rvec;

}

SEXP catnetKLpairwise(SEXP rSamples, SEXP rPerturbations) {

	int *pSamples, *pPerturbations;
	int *pSamplesPert, numsamplesPert;
	int *pNodeNumCats, **pNodeCats, mincat, maxcat, maxCategories;
	double *pprobs1, *pprobs2;
	int numsamples, numnodes, i, j, k, d, nnode1, nnode2;
	double floglik, faux, fsum, *pvec, *klmat;
	SEXP dim, rvec = R_NilValue;

	if(!isMatrix(rSamples))
		error("Data should be a matrix");
	if(!isNull(rPerturbations) && !isMatrix(rPerturbations))
		error("Perturbations should be a matrix");

	PROTECT(rSamples = AS_INTEGER(rSamples));
	pSamples = INTEGER(rSamples);

	dim = GET_DIM(rSamples);
	numnodes = INTEGER(dim)[0];
	numsamples = INTEGER(dim)[1];	

	if(isNull(rPerturbations)) {
		UNPROTECT(1); //rSamples
		PROTECT(rvec = NEW_NUMERIC(numnodes*numnodes));
		pvec = NUMERIC_POINTER(rvec);
		memset(pvec, 0, numnodes*numnodes*sizeof(double));
		UNPROTECT(1);
		return rvec;
	}

	// pSamples are assumed positive indices
	for(j = 0; j < numnodes*numsamples; j++) {
		pSamples[j]--;
	}

	// categoies
	pNodeNumCats = (int*)CATNET_MALLOC(numnodes*sizeof(int));
	pNodeCats = (int**)CATNET_MALLOC(numnodes*sizeof(int*));
	memset(pNodeCats, 0, numnodes*sizeof(int*));
	memset(pNodeNumCats, 0, numnodes*sizeof(int));

	maxCategories = 1;
	for(i = 0; i < numnodes; i++) {
		mincat = INT_MAX;
		maxcat = -INT_MAX;
		for(j = 0; j < numsamples; j++) {
			if(pSamples[j*numnodes + i] < mincat)
				mincat = pSamples[j*numnodes + i];
			if(pSamples[j*numnodes + i] > maxcat)
				maxcat = pSamples[j*numnodes + i];
		}
		pNodeNumCats[i] = maxcat - mincat + 1;
		pNodeCats[i] = (int*)CATNET_MALLOC(pNodeNumCats[i]*sizeof(int));
		for(j = 0; j < pNodeNumCats[i]; j++)
			pNodeCats[i][j] = mincat + j;
	}
	for(i = 0; i < numnodes; i++) {
		/* order pNodeNumCats[i] */
		for(j = 0; j < pNodeNumCats[i]; j++) {
			for(k = j + 1; k < pNodeNumCats[i]; k++) {
				if(pNodeCats[i][j] > pNodeCats[i][k]) {
					d = pNodeCats[i][j]; 
					pNodeCats[i][j] = pNodeCats[i][k];
					pNodeCats[i][k] = d;
				}
			}
		} 
		for(j = 0; j < numsamples; j++) {
			for(d = 0; d < pNodeNumCats[i]; d++)
				if(pNodeCats[i][d] == pSamples[j*numnodes + i])
					break;
			pSamples[j*numnodes + i] = d;
		}
		if(maxCategories < pNodeNumCats[i])
			maxCategories = pNodeNumCats[i];
	}

	pprobs1 = (double*)CATNET_MALLOC(maxCategories*maxCategories*sizeof(double));
	pprobs2 = (double*)CATNET_MALLOC(maxCategories*maxCategories*sizeof(double));

	pSamplesPert = 0;
	pPerturbations = 0;
	if(!isNull(rPerturbations)) {
		PROTECT(rPerturbations = AS_INTEGER(rPerturbations));
		pPerturbations = INTEGER(rPerturbations);
		pSamplesPert = (int*)CATNET_MALLOC(numnodes*numsamples*sizeof(int));
	}

	klmat = (double*)CATNET_MALLOC(numnodes*numnodes*sizeof(double));
	memset(klmat, 0, numnodes*numnodes*sizeof(double));

	for(nnode1 = 0; nnode1 < numnodes; nnode1++) {
		numsamplesPert = 0;
		if(pPerturbations) {
			for(j = 0; j < numsamples; j++) {
				if(!pPerturbations[j * numnodes + nnode1]) {
					memcpy(pSamplesPert + numsamplesPert*numnodes, pSamples + j*numnodes, numnodes*sizeof(int));
					numsamplesPert++;
				}
			}
		}
		//printf("\nnnode = %d (%d)\n", nnode1, numsamplesPert);
		for(nnode2 = 0; nnode2 < numnodes; nnode2++) {
	
			if(nnode1 == nnode2)
				continue;

			memset(pprobs2, 0, maxCategories*maxCategories*sizeof(double));
			// estimate logP(nnode1|nnode2) for the whole sample
			for(j = 0; j < numsamples; j++) 
				pprobs2[maxCategories*pSamples[j*numnodes + nnode2] + pSamples[j*numnodes + nnode1]]+=1;

			memset(pprobs1, 0, maxCategories*maxCategories*sizeof(double));
			// estimate logP(nnode1|nnode2) for the perturbed sub-sample only
			for(j = 0; j < numsamplesPert; j++) 
				pprobs1[maxCategories*pSamplesPert[j*numnodes + nnode2] + pSamplesPert[j*numnodes + nnode1]]+=1;

			for(i = 0; i < pNodeNumCats[nnode2]; i++) {
				fsum = 0;
				for(j = 0; j < pNodeNumCats[nnode1]; j++)
					fsum += pprobs1[maxCategories*i+j];
				if(fsum <= 0)
					continue;
				faux = 1 / fsum;
				for(j = 0; j < pNodeNumCats[nnode1]; j++)
					pprobs1[maxCategories*i+j] *= faux;
			}
			for(i = 0; i < pNodeNumCats[nnode2]; i++) {
				fsum = 0;
				for(j = 0; j < pNodeNumCats[nnode1]; j++)
					fsum += pprobs2[maxCategories*i+j];
				if(fsum <= 0)
					continue;
				faux = 1 / fsum;
				for(j = 0; j < pNodeNumCats[nnode1]; j++)
					pprobs2[maxCategories*i+j] *= faux;
			}

			floglik = 0;
			for(i = 0; i < pNodeNumCats[nnode2]; i++) {
				for(j = 0; j < pNodeNumCats[nnode1]; j++) {
					if(pprobs1[maxCategories*i+j] > 0 && pprobs2[maxCategories*i+j] > 0)
						floglik += pprobs1[maxCategories*i+j]*
						(double)log((double)pprobs1[maxCategories*i+j] / (double)pprobs2[maxCategories*i+j]);
					else if(pprobs1[maxCategories*i+j] != 0 && pprobs2[maxCategories*i+j] == 0)
						floglik = FLT_MAX;
				}
			}
			klmat[nnode2*numnodes + nnode1] += floglik;

			/*memset(pprobs2, 0, maxCategories*maxCategories*sizeof(double));
			// estimate logP(nnode2|nnode1) for the whole sample
			for(j = 0; j < numsamples; j++) 
				pprobs2[maxCategories*pSamples[j*numnodes + nnode1] + pSamples[j*numnodes + nnode2]]++;

			memset(pprobs1, 0, maxCategories*maxCategories*sizeof(double));
			// estimate logP(nnode2|nnode1) for the perturbed sub-sample only
			for(j = 0; j < numsamplesPert; j++) 
				pprobs1[maxCategories*pSamplesPert[j*numnodes + nnode1] + pSamplesPert[j*numnodes + nnode2]]++;

			for(i = 0; i < pNodeNumCats[nnode1]; i++) {
				fsum = 0;
				for(j = 0; j < pNodeNumCats[nnode2]; j++)
					fsum += pprobs1[maxCategories*i+j];
				if(fsum <= 0)
					continue;
				faux = 1 / fsum;
				for(j = 0; j < pNodeNumCats[nnode2]; j++)
					pprobs1[maxCategories*i+j] *= faux;
			}
			for(i = 0; i < pNodeNumCats[nnode1]; i++) {
				fsum = 0;
				for(j = 0; j < pNodeNumCats[nnode2]; j++)
					fsum += pprobs2[maxCategories*i+j];
				if(fsum <= 0)
					continue;
				faux = 1 / fsum;
				for(j = 0; j < pNodeNumCats[nnode2]; j++)
					pprobs2[maxCategories*i+j] *= faux;
			}

			floglik = 0;
			for(i = 0; i < pNodeNumCats[nnode1]; i++) {
				for(j = 0; j < pNodeNumCats[nnode2]; j++) {
					if(pprobs1[maxCategories*i+j] > 0 && pprobs2[maxCategories*i+j] > 0)
						floglik += pprobs1[maxCategories*i+j]*
						(double)log((double)pprobs1[maxCategories*i+j] / (double)pprobs2[maxCategories*i+j]);
					else if(pprobs1[maxCategories*i+j] != 0 && pprobs2[maxCategories*i+j] == 0)
						floglik = FLT_MAX;
				}
			}
			klmat[nnode2*numnodes + nnode1] += floglik;*/

			/*printf("nnode2 = %d: ", nnode2);
			for(i = 0; i < pNodeNumCats[nnode2]; i++) {
				printf(" [");
				for(j = 0; j < pNodeNumCats[nnode1]; j++) {
					printf("%f ", pprobs1[maxCategories*i+j]);
				}
				printf("] ");
			}
			printf("\n");
			for(i = 0; i < pNodeNumCats[nnode2]; i++) {
				printf(" [");
				for(j = 0; j < pNodeNumCats[nnode1]; j++) {
					printf("%f ", pprobs2[maxCategories*i+j]);
				}
				printf("] ");
			}
			printf(" KL = %f\n", floglik);*/
		}
	}

	UNPROTECT(1); //rSamples
	if(!isNull(rPerturbations))
		UNPROTECT(1); //rPerturbations

	if(pSamplesPert)
		CATNET_FREE(pSamplesPert);

	if(pprobs1)
		CATNET_FREE(pprobs1);
	if(pprobs2)
		CATNET_FREE(pprobs2);

	if(pNodeCats) {
		for(i = 0; i < numnodes; i++) 
			if(pNodeCats[i])
				CATNET_FREE(pNodeCats[i]);
		CATNET_FREE(pNodeCats);
	}

	if(pNodeNumCats) 
		CATNET_FREE(pNodeNumCats);

	PROTECT(rvec = NEW_NUMERIC(numnodes*numnodes));
	pvec = NUMERIC_POINTER(rvec);
	memcpy(pvec, klmat, numnodes*numnodes*sizeof(double));
	UNPROTECT(1);

	CATNET_FREE(klmat);

	//char str[128];
	//sprintf(str, "Mem Balance  %d\n", (int)g_memcounter);
	//printf(str);

	return rvec;

}

SEXP catnetPearsonPairwise(SEXP rSamples, SEXP rPerturbations) {

	int *pSamples, *pPerturbations;
	int *pSamplesPert, numsamplesPert;
	int *pNodeNumCats, **pNodeCats, mincat, maxcat, maxCategories;
	double *pprobs1, *pprobs2;
	int numsamples, numnodes, i, j, k, d, nnode1, nnode2;
	double floglik, faux, fsum, *pvec, *klmat;
	SEXP dim, rvec = R_NilValue;

	if(!isMatrix(rSamples))
		error("Data should be a matrix");
	if(!isNull(rPerturbations) && !isMatrix(rPerturbations))
		error("Perturbations should be a matrix");

	PROTECT(rSamples = AS_INTEGER(rSamples));
	pSamples = INTEGER(rSamples);

	dim = GET_DIM(rSamples);
	numnodes = INTEGER(dim)[0];
	numsamples = INTEGER(dim)[1];	

	if(isNull(rPerturbations)) {
		UNPROTECT(1); //rSamples
		PROTECT(rvec = NEW_NUMERIC(numnodes*numnodes));
		pvec = NUMERIC_POINTER(rvec);
		memset(pvec, 0, numnodes*numnodes*sizeof(double));
		UNPROTECT(1);
		return rvec;
	}

	// pSamples are assumed positive indices
	for(j = 0; j < numnodes*numsamples; j++) {
		pSamples[j]--;
	}

	// categoies
	pNodeNumCats = (int*)CATNET_MALLOC(numnodes*sizeof(int));
	pNodeCats = (int**)CATNET_MALLOC(numnodes*sizeof(int*));
	memset(pNodeCats, 0, numnodes*sizeof(int*));
	memset(pNodeNumCats, 0, numnodes*sizeof(int));

	maxCategories = 1;
	for(i = 0; i < numnodes; i++) {
		mincat = INT_MAX;
		maxcat = -INT_MAX;
		for(j = 0; j < numsamples; j++) {
			if(pSamples[j*numnodes + i] < mincat)
				mincat = pSamples[j*numnodes + i];
			if(pSamples[j*numnodes + i] > maxcat)
				maxcat = pSamples[j*numnodes + i];
		}
		pNodeNumCats[i] = maxcat - mincat + 1;
		pNodeCats[i] = (int*)CATNET_MALLOC(pNodeNumCats[i]*sizeof(int));
		for(j = 0; j < pNodeNumCats[i]; j++)
			pNodeCats[i][j] = mincat + j;
	}
	for(i = 0; i < numnodes; i++) {
		/* order pNodeNumCats[i] */
		for(j = 0; j < pNodeNumCats[i]; j++) {
			for(k = j + 1; k < pNodeNumCats[i]; k++) {
				if(pNodeCats[i][j] > pNodeCats[i][k]) {
					d = pNodeCats[i][j]; 
					pNodeCats[i][j] = pNodeCats[i][k];
					pNodeCats[i][k] = d;
				}
			}
		} 
		for(j = 0; j < numsamples; j++) {
			for(d = 0; d < pNodeNumCats[i]; d++)
				if(pNodeCats[i][d] == pSamples[j*numnodes + i])
					break;
			pSamples[j*numnodes + i] = d;
		}
		if(maxCategories < pNodeNumCats[i])
			maxCategories = pNodeNumCats[i];
	}

	pprobs1 = (double*)CATNET_MALLOC(maxCategories*maxCategories*sizeof(double));
	pprobs2 = (double*)CATNET_MALLOC(maxCategories*maxCategories*sizeof(double));

	pSamplesPert = 0;
	pPerturbations = 0;
	if(!isNull(rPerturbations)) {
		PROTECT(rPerturbations = AS_INTEGER(rPerturbations));
		pPerturbations = INTEGER(rPerturbations);
		pSamplesPert = (int*)CATNET_MALLOC(numnodes*numsamples*sizeof(int));
	}

	klmat = (double*)CATNET_MALLOC(numnodes*numnodes*sizeof(double));
	memset(klmat, 0, numnodes*numnodes*sizeof(double));

	for(nnode1 = 0; nnode1 < numnodes; nnode1++) {
		numsamplesPert = 0;
		if(pPerturbations) {
			for(j = 0; j < numsamples; j++) {
				if(!pPerturbations[j * numnodes + nnode1]) {
					memcpy(pSamplesPert + numsamplesPert*numnodes, pSamples + j*numnodes, numnodes*sizeof(int));
					numsamplesPert++;
				}
			}
		}
		//printf("\nnnode = %d (%d)\n", nnode1, numsamplesPert);
		for(nnode2 = 0; nnode2 < numnodes; nnode2++) {
	
			if(nnode1 == nnode2)
				continue;

			memset(pprobs2, 0, maxCategories*maxCategories*sizeof(double));
			// estimate logP(nnode1|nnode2) for the whole sample
			for(j = 0; j < numsamples; j++) 
				pprobs2[maxCategories*pSamples[j*numnodes + nnode2] + pSamples[j*numnodes + nnode1]]+=1;

			memset(pprobs1, 0, maxCategories*maxCategories*sizeof(double));
			// estimate logP(nnode1|nnode2) for the perturbed sub-sample only
			for(j = 0; j < numsamplesPert; j++) 
				pprobs1[maxCategories*pSamplesPert[j*numnodes + nnode2] + pSamplesPert[j*numnodes + nnode1]]+=1;

			for(i = 0; i < pNodeNumCats[nnode2]; i++) {
				fsum = 0;
				for(j = 0; j < pNodeNumCats[nnode1]; j++)
					fsum += pprobs2[maxCategories*i+j];
				if(fsum <= 0)
					continue;
				faux = 1 / fsum;
				for(j = 0; j < pNodeNumCats[nnode1]; j++)
					pprobs2[maxCategories*i+j] *= faux;
			}

			floglik = 0;
			for(i = 0; i < pNodeNumCats[nnode2]; i++) {
				fsum = 0;
				for(j = 0; j < pNodeNumCats[nnode1]; j++)
					fsum += pprobs1[maxCategories*i+j];
				if(fsum <= 0)
					continue;
				for(j = 0; j < pNodeNumCats[nnode1]; j++) {
					faux = pprobs1[maxCategories*i+j] - fsum*pprobs2[maxCategories*i+j];
					if(pprobs2[maxCategories*i+j] > 0)
						floglik += (double)(faux*faux) / (double)(fsum*pprobs2[maxCategories*i+j]);
					else if(faux != 0 && pprobs2[maxCategories*i+j] == 0)
						floglik = FLT_MAX;
				}
			}
			klmat[nnode2*numnodes + nnode1] += floglik;

			/*printf("nnode2 = %d: ", nnode2);
			for(i = 0; i < pNodeNumCats[nnode2]; i++) {
				printf(" [");
				for(j = 0; j < pNodeNumCats[nnode1]; j++) {
					printf("%f ", pprobs1[maxCategories*i+j]);
				}
				printf("] ");
			}
			printf("\n");
			for(i = 0; i < pNodeNumCats[nnode2]; i++) {
				printf(" [");
				for(j = 0; j < pNodeNumCats[nnode1]; j++) {
					printf("%f ", pprobs2[maxCategories*i+j]);
				}
				printf("] ");
			}
			printf(" KL = %f\n", floglik);*/
		}
	}

	UNPROTECT(1); //rSamples
	if(!isNull(rPerturbations))
		UNPROTECT(1); //rPerturbations

	if(pSamplesPert)
		CATNET_FREE(pSamplesPert);

	if(pprobs1)
		CATNET_FREE(pprobs1);
	if(pprobs2)
		CATNET_FREE(pprobs2);

	if(pNodeCats) {
		for(i = 0; i < numnodes; i++) 
			if(pNodeCats[i])
				CATNET_FREE(pNodeCats[i]);
		CATNET_FREE(pNodeCats);
	}

	if(pNodeNumCats) 
		CATNET_FREE(pNodeNumCats);

	PROTECT(rvec = NEW_NUMERIC(numnodes*numnodes));
	pvec = NUMERIC_POINTER(rvec);
	memcpy(pvec, klmat, numnodes*numnodes*sizeof(double));
	UNPROTECT(1);

	CATNET_FREE(klmat);

	//char str[128];
	//sprintf(str, "Mem Balance  %d\n", (int)g_memcounter);
	//printf(str);

	return rvec;

}

} // extern "C"
