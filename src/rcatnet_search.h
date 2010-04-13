/*
 * rcatnet.h
 *
 *  Created on: Sep 21, 2009
 *      Author: nbalov
 */

#ifndef RCATNET_SEARCH_H
#define RCATNET_SEARCH_H

#include "catnet_search.h"

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

#define MAX_NODE_NAME	16

extern CATNET_CACHE_EL<double>		**g_pcache;
extern unsigned int			g_ncache;

void ReleaseCache();
void InitializeCache(int nslots, int ncachebits);

class RCatnetSearch : public CATNET_SEARCH<char, MAX_NODE_NAME, double> {
protected:
	int m_maxParentSet, *m_pRorder, *m_pRorderInverse, *m_parBuff1, *m_parBuff2;
	int m_bUseCache, m_nCacheBits;

public:
	RCatnetSearch();
	~RCatnetSearch();

	int getCachedProb(int *ppool, int poolsize, int node, int *parset, int parsize, 
						PROB_LIST<double> *probNode, double *pflik);
	int setCachedProb(int *ppool, int poolsize, int node, int *parset, int parsize, 
						PROB_LIST<double> *probNode, double flik);

	SEXP estimateCatnets(SEXP rSamples, SEXP rPerturbations,
                       SEXP rMaxParents, SEXP rMaxComplexity, SEXP rOrder,
                       SEXP rParentsPool, SEXP rFixedParentsPool, SEXP rUseCache, SEXP rEcho);

};

#endif /* RCATNET_SEARCH_H */
