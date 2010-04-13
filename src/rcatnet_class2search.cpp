/*
 * rcatnet.cpp
 *
 *  Created on: Sep 21, 2009
 *      Author: nbalov
 */

#include "utils.h"
#include "catnet_class.h"
#include "rcatnet_class2search.h"
#include "rcatnet.h"

RCatnetClass2Search::RCatnetClass2Search() {
	m_pRorder = 0;
	m_pRorderInverse = 0;
	m_parBuff1 = 0;
	m_parBuff2 = 0;
	m_bUseCache = 1;
}

RCatnetClass2Search::~RCatnetClass2Search() {
	if(m_pRorder)
		CATNET_FREE(m_pRorder);
	m_pRorder = 0;
	if(m_pRorderInverse)
		CATNET_FREE(m_pRorderInverse);
	m_pRorderInverse = 0;
	if(m_parBuff1)
		CATNET_FREE(m_parBuff1);
	m_parBuff1 = 0;
	if(m_parBuff2)
		CATNET_FREE(m_parBuff2);
	m_parBuff2 = 0;
}

SEXP RCatnetClass2Search::estimateCatnets(SEXP rSamples1, SEXP rPerturbations1,
			SEXP rSamples2, SEXP rPerturbations2, 
			SEXP rMaxParents, SEXP rMaxComplexity, SEXP rOrder,
			SEXP rParentsPool, SEXP rFixedParentsPool, SEXP rUseCache, SEXP rEcho) {

	int i, j, k, len, maxComplexity, numnets, inet, echo;
 	int *pRsamples1, *pRsamples2, *pRperturbations1, *pRperturbations2, *pSamples1, *pPerturbations1, *pSamples2, *pPerturbations2, **parentsPool, **fixedParentsPool, *pPool;

	RCatnet rcatnet;
	SEXP dim, rparpool, cnetlist, cnetnode;

	if(!isMatrix(rSamples1) || !isMatrix(rSamples2))
		error("Data is not a matrix");

	PROTECT(rSamples1 = AS_INTEGER(rSamples1));
	PROTECT(rSamples2 = AS_INTEGER(rSamples2));
	PROTECT(rMaxParents = AS_INTEGER(rMaxParents));
	PROTECT(rMaxComplexity = AS_INTEGER(rMaxComplexity));
	PROTECT(rOrder = AS_INTEGER(rOrder));

	PROTECT(rUseCache = AS_LOGICAL(rUseCache));
	m_bUseCache = LOGICAL(rUseCache)[0];
	//printf("bUseCache = %d\n", m_bUseCache);
	UNPROTECT(1);

	PROTECT(rEcho = AS_LOGICAL(rEcho));
	echo = LOGICAL(rEcho)[0];
	UNPROTECT(1);

	dim = GET_DIM(rSamples1);
	m_numNodes = INTEGER(dim)[0];
	m_numSamples1 = INTEGER(dim)[1];
	dim = GET_DIM(rSamples2);
	if(m_numNodes != INTEGER(dim)[0]) {
		UNPROTECT(5);
		error("Samples should have equal number of nodes");
	}
	m_numSamples2 = INTEGER(dim)[1];

        m_maxParentSet = INTEGER_POINTER(rMaxParents)[0];
	if(m_parBuff1)
		CATNET_FREE(m_parBuff1);
	m_parBuff1 = (int*)CATNET_MALLOC(m_numNodes*sizeof(int));
	if(m_parBuff2)
		CATNET_FREE(m_parBuff2);
	m_parBuff2 = (int*)CATNET_MALLOC(m_numNodes*sizeof(int));

	maxComplexity = INTEGER_POINTER(rMaxComplexity)[0];

	pSamples1 = (int*)CATNET_MALLOC(m_numNodes*m_numSamples1*sizeof(int));
	pRsamples1 = INTEGER(rSamples1);
	pSamples2 = (int*)CATNET_MALLOC(m_numNodes*m_numSamples2*sizeof(int));
	pRsamples2 = INTEGER(rSamples2);

	if(m_pRorder)
		CATNET_FREE(m_pRorder);
	m_pRorder = (int*)CATNET_MALLOC(m_numNodes*sizeof(int));
	if(m_pRorderInverse)
		CATNET_FREE(m_pRorderInverse);
	m_pRorderInverse = (int*)CATNET_MALLOC(m_numNodes*sizeof(int));

	if(length(rOrder) < m_numNodes) {
		warning("Invalid nodeOrder parameter - reset to default node order.");
		for(i = 0; i < m_numNodes; i++)
			m_pRorder[i] = i + 1;
	}
	else {
		memcpy(m_pRorder, INTEGER(rOrder), m_numNodes*sizeof(int));
	}
	for(i = 0; i < m_numNodes; i++)
		m_pRorderInverse[m_pRorder[i]-1] = i + 1;

	for(j = 0; j < m_numSamples1; j++) {
		for(i = 0; i < m_numNodes; i++) {
			pSamples1[j*m_numNodes + i] = pRsamples1[j*m_numNodes + m_pRorder[i] - 1];
		}
	}

	for(j = 0; j < m_numSamples2; j++) {
		for(i = 0; i < m_numNodes; i++) {
			pSamples2[j*m_numNodes + i] = pRsamples2[j*m_numNodes + m_pRorder[i] - 1];
		}
	}

	pPerturbations1 = 0;
	if(!isNull(rPerturbations1)) {
		PROTECT(rPerturbations1 = AS_INTEGER(rPerturbations1));
		pPerturbations1 = (int*)CATNET_MALLOC(m_numNodes*m_numSamples1*sizeof(int));
		pRperturbations1 = INTEGER(rPerturbations1);
		for(j = 0; j < m_numSamples1; j++) {
			for(i = 0; i < m_numNodes; i++) {
				pPerturbations1[j*m_numNodes + i] = pRperturbations1[j*m_numNodes + m_pRorder[i] - 1];
			}
		}
		UNPROTECT(1);
	}

	pPerturbations2 = 0;
	if(!isNull(rPerturbations2)) {
		PROTECT(rPerturbations2 = AS_INTEGER(rPerturbations2));
		pPerturbations2 = (int*)CATNET_MALLOC(m_numNodes*m_numSamples2*sizeof(int));
		pRperturbations2 = INTEGER(rPerturbations2);
		for(j = 0; j < m_numSamples2; j++) {
			for(i = 0; i < m_numNodes; i++) {
				pPerturbations2[j*m_numNodes + i] = pRperturbations2[j*m_numNodes + m_pRorder[i] - 1];
			}
		}
		UNPROTECT(1);
	}

	parentsPool = 0;
	if(!isNull(rParentsPool)) {
		PROTECT(rParentsPool = AS_LIST(rParentsPool));

		parentsPool = (int**)CATNET_MALLOC(m_numNodes*sizeof(int*));
		memset(parentsPool, 0, m_numNodes*sizeof(int*));
		for(i = 0; i < m_numNodes; i++) {
			rparpool = AS_INTEGER(VECTOR_ELT(rParentsPool, (int)(m_pRorder[i] - 1)));
			len = length(rparpool);
			if(isVector(rparpool) && len > 0 && len <= m_numNodes) {
				pPool = INTEGER(rparpool);
				parentsPool[i] = (int*)CATNET_MALLOC(m_numNodes*sizeof(int));
				for(j = 0; j < len; j++) {
					if(pPool[j] > 0 && pPool[j] <= m_numNodes) {
						for(k = 0; k < m_numNodes; k++)
							if(pPool[j] == m_pRorder[k])
								break;
						if(k < m_numNodes)
							parentsPool[i][j] = k;
						else
							parentsPool[i][j] = -1;
					}
				}
				for(; j < m_numNodes; j++)
					parentsPool[i][j] = -1;
			}
		}
		UNPROTECT(1);
	}

	fixedParentsPool = 0;
	if(!isNull(rFixedParentsPool)) {
		PROTECT(rFixedParentsPool = AS_LIST(rFixedParentsPool));

		fixedParentsPool = (int**)CATNET_MALLOC(m_numNodes*sizeof(int*));
		memset(fixedParentsPool, 0, m_numNodes*sizeof(int*));
		for(i = 0; i < m_numNodes; i++) {
			rparpool = AS_INTEGER(VECTOR_ELT(rFixedParentsPool, (int)(m_pRorder[i] - 1)));
			len = length(rparpool);
			if(isVector(rparpool) && len > 0 && len <= m_numNodes) {
			 	if(m_maxParentSet < len)
			    		m_maxParentSet = len;
				pPool = INTEGER(rparpool);
				fixedParentsPool[i] = (int*)CATNET_MALLOC(m_numNodes*sizeof(int));
				for(j = 0; j < len; j++) {
					if(pPool[j] > 0 && pPool[j] <= m_numNodes) {
						for(k = 0; k < m_numNodes; k++)
							if(pPool[j] == m_pRorder[k])
								break;
						if(k < m_numNodes)
							fixedParentsPool[i][j] = k;
						else
							fixedParentsPool[i][j] = -1;
					}
				}
				for(; j < m_numNodes; j++)
					fixedParentsPool[i][j] = -1;
			}
		}
		UNPROTECT(1);
	}

	UNPROTECT(5);

printf("estimate %d, %d\n", m_numSamples1, m_numSamples2);

	estimate(m_numNodes, m_numSamples1, pSamples1, pPerturbations1, m_numSamples2, pSamples2, pPerturbations2, m_maxParentSet, maxComplexity, parentsPool, fixedParentsPool, echo);

	if(pSamples1)
		CATNET_FREE(pSamples1);
	if(pPerturbations1)
		CATNET_FREE(pPerturbations1);
	if(pSamples2)
		CATNET_FREE(pSamples2);
	if(pPerturbations2)
		CATNET_FREE(pPerturbations2);
	if(parentsPool) {
		for(i = 0; i < m_numNodes; i++)
			if(parentsPool[i])
				CATNET_FREE(parentsPool[i]);
		CATNET_FREE(parentsPool);
	}
	if(fixedParentsPool) {
		for(i = 0; i < m_numNodes; i++)
			if(fixedParentsPool[i])
				CATNET_FREE(fixedParentsPool[i]);
		CATNET_FREE(fixedParentsPool);
	}

	if(!m_nCatnets || !m_pCatnets1)
		return R_NilValue;

	// create a R-list of catNetworks
	numnets = 0;
	for(i = 0; i < m_nCatnets; i++) {
		if(m_pCatnets1[i]) {
			m_pCatnets1[i]->setNodesOrder(m_pRorder);
			numnets++;
		}
	}

	PROTECT(cnetlist = allocVector(VECSXP, numnets));

	inet = 0;
	for(i = 0; i < m_nCatnets; i++) {
		if(!m_pCatnets1[i])
			continue;

		rcatnet = *m_pCatnets1[i];

		rcatnet.setCategoryIndices(m_pNodeNumCats, m_pNodeCats);

		PROTECT(cnetnode = rcatnet.genRcatnet("catNetwork"));

		SET_VECTOR_ELT(cnetlist, inet, cnetnode);
		UNPROTECT(1);
		inet++;
	}

	UNPROTECT(1);

	return cnetlist;
}

