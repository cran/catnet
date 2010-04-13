/*
 * rcatnet_class2search.h
 *
 *  Created on: Sep 21, 2009
 *      Author: nbalov
 */

#ifndef RCATNET_CLASS2SEARCH_H
#define RCATNET_CLASS2SEARCH_H

#include "catnet_class2search.h"

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

#define MAX_NODE_NAME	16

class RCatnetClass2Search : public CATNET_CLASS2SEARCH<char, MAX_NODE_NAME, double> {
protected:
	int m_maxParentSet, *m_pRorder, *m_pRorderInverse, *m_parBuff1, *m_parBuff2, m_bUseCache;

public:
	RCatnetClass2Search();
	~RCatnetClass2Search();

	SEXP estimateCatnets(SEXP rSamples1, SEXP rPerturbations1,
			SEXP rSamples2, SEXP rPerturbations2,
			SEXP rMaxParents, SEXP rMaxComplexity, SEXP rOrder,
			SEXP rParentsPool, SEXP rFixedParentsPool, SEXP rUseCache, SEXP rEcho);

};

#endif /* RCATNET_CLASS2SEARCH_H */
