/*
 * rcatnet.cpp
 *
 *  Created on: Sep 21, 2009
 *      Author: nbalov
 */

#include "utils.h"
#include "catnet_class.h"
#include "rcatnet_search.h"
#include "rcatnet.h"

#define CACHE_MAX_LEN_POOL 50
int PRIMES_NUM = 180;
unsigned PRIMES_1000[] = {
      2,      3,      5,      7,     11,     13,     17,     19,     23,     29, 
     31,     37,     41,     43,     47,     53,     59,     61,     67,     71, 
     73,     79,     83,     89,     97,    101,    103,    107,    109,    113, 
    127,    131,    137,    139,    149,    151,    157,    163,    167,    173, 
    179,    181,    191,    193,    197,    199,    211,    223,    227,    229, 
    233,    239,    241,    251,    257,    263,    269,    271,    277,    281, 
    283,    293,    307,    311,    313,    317,    331,    337,    347,    349, 
    353,    359,    367,    373,    379,    383,    389,    397,    401,    409, 
    419,    421,    431,    433,    439,    443,    449,    457,    461,    463, 
    467,    479,    487,    491,    499,    503,    509,    521,    523,    541, 
    547,    557,    563,    569,    571,    577,    587,    593,    599,    601, 
    607,    613,    617,    619,    631,    641,    643,    647,    653,    659, 
    661,    673,    677,    683,    691,    701,    709,    719,    727,    733, 
    739,    743,    751,    757,    761,    769,    773,    787,    797,    809, 
    811,    821,    823,    827,    829,    839,    853,    857,    859,    863, 
    877,    881,    883,    887,    907,    911,    919,    929,    937,    941, 
    947,    953,    967,    971,    977,    983,    991,    997,   1009,   1013, 
   1019,   1021,   1031,   1033,   1039,   1049,   1051,   1061,   1063,   1069 
};

//#define DEBUG_INFO

CATNET_CACHE_EL<double>		**g_pcache = 0;
unsigned int			g_ncache = 0;
unsigned int			g_nCacheBits = 0;

void ReleaseCache() {
	unsigned i;
	//printf("\nRELEASE CACHE\n");
	if(g_pcache && g_ncache > 0) {
		for(i = 0; i < g_ncache; i++) {
			if(g_pcache[i])
				delete g_pcache[i];
			g_pcache[i] = 0;
		}
	}
	if(g_pcache)
		CATNET_FREE(g_pcache);
	g_pcache = 0;
	g_ncache = 0;
	g_nCacheBits = 0;
}
	
void InitializeCache(int nslots, int ncachebits) {
	ReleaseCache();
	//printf("\nINITIALIZE CACHE\n");
	if(nslots < 1)
		nslots = 1;
	g_ncache = nslots;
	g_pcache = (CATNET_CACHE_EL<double>**)CATNET_MALLOC(g_ncache*sizeof(CATNET_CACHE_EL<double>*));
	memset(g_pcache, 0, g_ncache*sizeof(CATNET_CACHE_EL<double>*));
	g_nCacheBits = ncachebits;
}


RCatnetSearch::RCatnetSearch() {
	m_pRorder = 0;
	m_pRorderInverse = 0;
	m_parBuff1 = 0;
	m_parBuff2 = 0;
	m_bUseCache = 1;
}

RCatnetSearch::~RCatnetSearch() {
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

int RCatnetSearch::getCachedProb(int *ppool, int poolsize, int node, int *parset, int parsize, 
				PROB_LIST<double> *probNode, double *pflik) {

		int i, j;
		unsigned nlookup;
		CATNET_CACHE_EL<double> *pCacheEl = 0;

		if(!m_bUseCache)
			return 0;
		if(!g_pcache)
			return 0;
#ifdef DEBUG_INFO
char str[256];
sprintf(str,"getCachedProb node=%d, pool = ", node);printf(str);
for(i = 0; i < poolsize; i++) {
	sprintf(str,"%d  ", ppool[i]);printf(str);
}printf("\n");
#endif
		node = m_pRorder[node];
		for(i = 0; i < poolsize; i++)
			m_parBuff1[i] = m_pRorder[ppool[i]];
		_quick_sort<int>(m_parBuff1, poolsize);

#ifdef DEBUG_INFO
sprintf(str,"    reordered node=%d, pool = ", node);printf(str);
for(i = 0; i < poolsize; i++) {
	sprintf(str,"%d  ", m_parBuff1[i]); printf(str);
}printf("\n");
#endif
		nlookup = 1;
		for(i = 0; i < poolsize; i++) {
			j = m_parBuff1[i] - 1;
			while(j >= PRIMES_NUM)
				j -= PRIMES_NUM;
			nlookup *= PRIMES_1000[PRIMES_NUM - j - 1];
			while(nlookup >= g_ncache)
				nlookup -= g_ncache;		
		}
		nlookup = (nlookup << g_nCacheBits) + node + m_numNodes*parsize;
		while(nlookup >= g_ncache)
			nlookup -= g_ncache;

#ifdef DEBUG_INFO
sprintf(str,"nlookup=%d\n", nlookup);printf(str);
#endif

		pCacheEl = g_pcache[nlookup];
		if(!pCacheEl) 
			return 0;
		if(pCacheEl->nnode != node)
			return 0;
		if(pCacheEl->npars != parsize)
			return 0;
		if(pCacheEl->nPool != poolsize)
			return 0;
		for(i = 0; i < poolsize; i++) {
			if(pCacheEl->pPool[i] != m_parBuff1[i])
				return 0;
		}

		for(i = 0; i < parsize; i++) {
			parset[i] = m_pRorderInverse[pCacheEl->pPars[i] - 1] - 1;
		}
		*probNode = *pCacheEl->pNodeProb;
		*pflik = pCacheEl->fLogLik;

 //clock_t lt = clock();
 //printf("%ld\n", lt);;

#ifdef DEBUG_INFO
printf("\n    HIT, parset = ");
for(i = 0; i < parsize; i++) {
	sprintf(str,"%d  ", parset[i]);printf(str);
}printf("\n");
#endif
		return 1;
	}

int RCatnetSearch::setCachedProb(int *ppool, int poolsize, int node, int *parset, int parsize, 
				PROB_LIST<double> *probNode, double flik) {
		int i, j, prime;
		unsigned nlookup;

		if(!m_bUseCache)
			return 0;
//char str[256];

		if(!g_pcache) {
			nlookup = m_numNodes;
			for(i = 0; i < m_maxParentSet; i++)
				nlookup *= m_numNodes;
			if(nlookup < PRIMES_1000[PRIMES_NUM-1]) {
				for(i = 0; i < PRIMES_NUM; i++)
					if(nlookup >= PRIMES_1000[i]) {
						nlookup = PRIMES_1000[i];
						break;
					}
			}
			else {
				prime = PRIMES_1000[PRIMES_NUM-1];
				for(i = 0; i < PRIMES_NUM; i++) {
					if(PRIMES_1000[i] >= (unsigned)m_numNodes) {
						prime = PRIMES_1000[i];
						break;
					}
				}
				nlookup = prime;
				for(i = 0; i < m_maxParentSet; i++)
					nlookup *= prime;
			}
			if(nlookup < PRIMES_1000[10])
				nlookup = PRIMES_1000[10];
			
			g_nCacheBits = 1;
			i = 1;
			while(i < m_numNodes*m_maxParentSet) {
				g_nCacheBits++;
				i <<= 1;
			}
#ifdef DEBUG_INFO
sprintf(str,"nCacheBits = %d  ", g_nCacheBits); printf(str);
#endif
			InitializeCache(nlookup, g_nCacheBits);
		}

#ifdef DEBUG_INFO
sprintf(str,"setCachedProb node=%d, poolsize = %d, parsize = %d, pool = ", node, poolsize, parsize); printf(str);
for(i = 0; i < poolsize; i++) {
	sprintf(str,"%d  ", ppool[i]);printf(str);
}printf("\n");
#endif

		node = m_pRorder[node];
		for(i = 0; i < poolsize; i++)
			m_parBuff1[i] = m_pRorder[ppool[i]];
		_quick_sort<int>(m_parBuff1, poolsize);

		for(i = 0; i < parsize; i++)
			m_parBuff2[i] = m_pRorder[parset[i]];

#ifdef DEBUG_INFO
sprintf(str,"    reordered node=%d, pool = ", node);printf(str);
for(i = 0; i < poolsize; i++) {
	sprintf(str,"%d  ", m_parBuff1[i]);printf(str);
}printf("\n par = ");
for(i = 0; i < parsize; i++) {
	sprintf(str,"%d  ", m_parBuff2[i]);printf(str);
}printf("\n");
#endif
		CATNET_CACHE_EL<double> *pCacheEl = new CATNET_CACHE_EL<double>(
			m_parBuff1, poolsize, node, m_parBuff2, parsize, probNode, flik);

		nlookup = 1;
		for(i = 0; i < poolsize; i++) {
			j = m_parBuff1[i] - 1;
			while(j >= PRIMES_NUM)
				j -= PRIMES_NUM;
			nlookup *= PRIMES_1000[PRIMES_NUM - j - 1];
			while(nlookup >= g_ncache)
				nlookup -= g_ncache;		
		}
		nlookup = (nlookup << g_nCacheBits) + node + m_numNodes*parsize;
		while(nlookup >= g_ncache)
			nlookup -= g_ncache;	

#ifdef DEBUG_INFO
sprintf(str,"\nnlookup=%d\n", nlookup);printf(str);
#endif
		if(g_pcache[nlookup]) {

/*			j = 0;
			for(i = 0; i < g_ncache; i++)
				if(g_pcache[i]) j++;
printf("fill %d/%d\n", j, g_ncache);
char str[256];
sprintf(str,"  Overwrite node=%d, numpars=%d, pool = ", node, parsize);printf(str);
for(i = 0; i < poolsize; i++) {
	sprintf(str,"%d  ", m_parBuff1[i]);printf(str);
}printf("\n");
sprintf(str,"  Old node=%d, numpars=%d, pool = ", g_pcache[nlookup]->nnode, g_pcache[nlookup]->npars);printf(str);
for(i = 0; i < g_pcache[nlookup]->nPool; i++) {
	sprintf(str,"%d  ", g_pcache[nlookup]->pPool[i]);printf(str);
}printf("\n");*/
			delete g_pcache[nlookup];
		}
		g_pcache[nlookup] = pCacheEl;

  //clock_t lt = clock();
  //printf("%ld\n", lt);

		return 1;
	}

SEXP RCatnetSearch::estimateCatnets(SEXP rSamples, SEXP rPerturbations, 
                       SEXP rMaxParents, SEXP rMaxComplexity, SEXP rOrder,
                       SEXP rParentsPool, SEXP rFixedParentsPool, SEXP rUseCache, SEXP rEcho) {

	int i, j, k, len, maxComplexity, numnets, inet, echo;
 	int *pRsamples, *pRperturbations, *pSamples, *pPerturbations, **parentsPool, **fixedParentsPool, *pPool;

	RCatnet rcatnet;
	SEXP dim, rparpool, cnetlist, cnetnode;

	if(!isMatrix(rSamples))
		error("Data is not a matrix");

	PROTECT(rSamples = AS_INTEGER(rSamples));
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

	dim = GET_DIM(rSamples);
	m_numNodes = INTEGER(dim)[0];
	m_numSamples = INTEGER(dim)[1];
 
        m_maxParentSet = INTEGER_POINTER(rMaxParents)[0];
	if(m_parBuff1)
		CATNET_FREE(m_parBuff1);
	m_parBuff1 = (int*)CATNET_MALLOC(m_numNodes*sizeof(int));
	if(m_parBuff2)
		CATNET_FREE(m_parBuff2);
	m_parBuff2 = (int*)CATNET_MALLOC(m_numNodes*sizeof(int));

	maxComplexity = INTEGER_POINTER(rMaxComplexity)[0];

	pSamples = (int*)CATNET_MALLOC(m_numNodes*m_numSamples*sizeof(int));
	pRsamples = INTEGER(rSamples);

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

	for(j = 0; j < m_numSamples; j++) {
		for(i = 0; i < m_numNodes; i++) {
			pSamples[j*m_numNodes + i] = pRsamples[j*m_numNodes + m_pRorder[i] - 1];
		}
	}

#ifdef DEBUG_INFO
char str[256];
printf("rOrder: \n");
for(i = 0; i < m_numNodes; i++) {
	sprintf(str,"%d  ",m_pRorder[i]);printf(str);
}printf("\n");
printf("rOrderInverse: \n");
for(i = 0; i < m_numNodes; i++) {
	sprintf(str,"%d  ",m_pRorderInverse[i]);printf(str);
}printf("\n");
#endif

	pPerturbations = 0;
	if(!isNull(rPerturbations)) {
		PROTECT(rPerturbations = AS_INTEGER(rPerturbations));
		pPerturbations = (int*)CATNET_MALLOC(m_numNodes*m_numSamples*sizeof(int));
		pRperturbations = INTEGER(rPerturbations);
		for(j = 0; j < m_numSamples; j++) {
			for(i = 0; i < m_numNodes; i++) {
				pPerturbations[j*m_numNodes + i] = pRperturbations[j*m_numNodes + m_pRorder[i] - 1];
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

	UNPROTECT(4);

	estimate(m_numNodes, m_numSamples, pSamples, pPerturbations, m_maxParentSet, maxComplexity, parentsPool, fixedParentsPool, echo);

	if(pSamples)
		CATNET_FREE(pSamples);
	if(pPerturbations)
		CATNET_FREE(pPerturbations);
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

	if(!m_nCatnets || !m_pCatnets)
		return R_NilValue;

	// create a R-list of catNetworks
	numnets = 0;
	for(i = 0; i < m_nCatnets; i++) {
		if(m_pCatnets[i]) {
			m_pCatnets[i]->setNodesOrder(m_pRorder);
			numnets++;
		}
	}

	PROTECT(cnetlist = allocVector(VECSXP, numnets));

	inet = 0;
	for(i = 0; i < m_nCatnets; i++) {
		if(!m_pCatnets[i])
			continue;

		rcatnet = *m_pCatnets[i];

		rcatnet.setCategoryIndices(m_pNodeNumCats, m_pNodeCats);

		PROTECT(cnetnode = rcatnet.genRcatnet("catNetwork"));

		SET_VECTOR_ELT(cnetlist, inet, cnetnode);
		UNPROTECT(1);
		inet++;
	}

	UNPROTECT(1);

	return cnetlist;
}

