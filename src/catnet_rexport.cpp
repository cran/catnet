/*
 * rcatnet.cpp
 *
 *  Created on: Sep 21, 2009
 *      Author: nbalov
 */

#include "utils.h"
#include "catnet_class.h"
#include "rcatnet.h"
#include "rcatnet_search.h"
#include "rcatnet_class2search.h"

extern "C" {

extern size_t g_memcounter;

SEXP catnetReleaseCache()
{
	ReleaseCache();
	return R_NilValue;
}

SEXP createRCatnet(SEXP cnet)
{
	SEXP pcnet = R_NilValue;
	PROTECT(cnet);
	RCatnet *rnet = new RCatnet(cnet);
	UNPROTECT(1);
	pcnet = rnet->genRcatnet((const char*)"catNetwork");
	delete rnet;
	return pcnet;
}

SEXP catnetMarginalProb(SEXP cnet, SEXP rnode)
{
	int node, i, ncats;
	SEXP rvec = R_NilValue;
	double *pvec;

	PROTECT(rnode = AS_INTEGER(rnode));
	node = INTEGER_VALUE(rnode);
	UNPROTECT(1);

	PROTECT(cnet);
	RCatnet *rnet = new RCatnet(cnet);
	UNPROTECT(1);

	if(!rnet)
		return rvec;

	if(node < 1 || node > rnet->numNodes())
		return rvec;
	node--;

	double *pmarg = rnet->marginal_prob(node);
	if(!pmarg)
		return rvec;

	ncats = rnet->numCategories(node);
	PROTECT(rvec = NEW_NUMERIC(ncats));
	pvec = NUMERIC_POINTER(rvec);
	for(i = 0; i < ncats; i++) {
		pvec[i] = pmarg[i];
	}
	UNPROTECT(1);

	CATNET_FREE(pmarg);
	delete rnet;

	//char str[128];
	//sprintf(str, "Mem Balance  %d\n", (int)g_memcounter);
	//printf(str);

	return rvec;
}

SEXP catnetJointProb(SEXP cnet, SEXP rnode)
{
	int node, jointprobsize;
	SEXP rvec = R_NilValue;
	double *pvec;

	PROTECT(rnode = AS_INTEGER(rnode));
	node = INTEGER_VALUE(rnode);
	UNPROTECT(1);

	PROTECT(cnet);
	RCatnet *rnet = new RCatnet(cnet);
	UNPROTECT(1);

	if(!rnet)
		return rvec;
	if(node < 1 || node > rnet->numNodes())
		return rvec;
	node--;

	jointprobsize = 0;
	double *pjoint = rnet->findJointProb(node, jointprobsize);
	if(!pjoint) {
		delete rnet;
		return rvec;
	}

	PROTECT(rvec = NEW_NUMERIC(jointprobsize));
	pvec = NUMERIC_POINTER(rvec);
	memcpy(pvec, pjoint, jointprobsize*sizeof(double));
	UNPROTECT(1);

	CATNET_FREE(pjoint);

	delete rnet;
	return rvec;
}

SEXP catnetFindParentPool(SEXP cnet, SEXP rnode)
{
	int node, i, poolsize;
	SEXP rvec = R_NilValue;
	int *ppool, *pvec;

	PROTECT(rnode = AS_INTEGER(rnode));
	node = INTEGER_VALUE(rnode);
	UNPROTECT(1);

	PROTECT(cnet);
	RCatnet *rnet = new RCatnet(cnet);
	UNPROTECT(1);

	if(!rnet)
		return rvec;

	if(node < 1 || node > rnet->numNodes())
		return rvec;
	node--;

	ppool = rnet->findParentPool(node, poolsize);
	if (!ppool) {
		delete rnet;
		return rvec;
	}

	PROTECT(rvec = NEW_INTEGER(poolsize));
	pvec = INTEGER_POINTER(rvec);
	for(i = 0; i < poolsize; i++) {
		pvec[i] = ppool[i]+1;
	}
	UNPROTECT(1);

	CATNET_FREE(ppool);

	delete rnet;
	return rvec;
}

SEXP show_catnet(SEXP rnodes, SEXP rparents, SEXP rcatlist, SEXP rproblist)
{
	int i, m_numNodes, nnode;
	char str[1024];
	SEXP pf, pstr;

	PROTECT(rnodes = AS_LIST(rnodes));
	PROTECT(rparents = AS_LIST(rparents));
	PROTECT(rcatlist = AS_LIST(rcatlist));
	PROTECT(rproblist = AS_LIST(rproblist));

	//printf("Call show_object.\n");

	PROTECT(pstr = allocVector(STRSXP, 3));

	m_numNodes = length(rnodes);
	sprintf(str, "Nodes = %d: ", m_numNodes);
	//printf(str);

	for(nnode = 0; nnode < m_numNodes; nnode++) {
		PROTECT(pf = VECTOR_ELT(rnodes, nnode));
		if(IS_VECTOR(pf)) {
			sprintf(str, "%s%s, ", str, CHAR(STRING_ELT(pf, 0)));
			//printf(str);
		}
		else {
			//printf("");
		}
		UNPROTECT(1);
	}
	sprintf(str, "%s\n", str);
	SET_STRING_ELT(pstr, 0, mkChar(str));

	//printf("\n\nParents:\n");
	sprintf(str, "Parents:\n");

	for(nnode = 0; nnode < m_numNodes; nnode++) {
		PROTECT(pf = VECTOR_ELT(rparents, nnode));
		sprintf(str, "%s[%d] ", str, nnode);
		if(IS_VECTOR(pf)) {
			for(i = 0; i < length(pf); i++)
				sprintf(str, "%s%d, ", str, INTEGER_POINTER(pf)[i]-1);
			sprintf(str, "%s\n", str);
			//printf(str);
		}
		else
			sprintf(str, "%s\n", str);
		UNPROTECT(1);
	}
	SET_STRING_ELT(pstr, 1, mkChar(str));

	//printf("\nCategories:\n");
	sprintf(str, "Categories:\n");

	for(nnode = 0; nnode < m_numNodes; nnode++) {
		PROTECT(pf = VECTOR_ELT(rcatlist, nnode));
		if(IS_VECTOR(pf)) {
			for(i = 0; i < length(pf); i++)
				sprintf(str, "%s%s, ", str, CHAR(STRING_ELT(pf,i)));
			sprintf(str, "%s\n", str);
			//printf(str);
		}
		else {
			//printf("null\n");
		}
		UNPROTECT(1);
	}
	SET_STRING_ELT(pstr, 2, mkChar(str));

	UNPROTECT(5);

	return pstr;
}

SEXP showCatnet(SEXP cnet)
{
	SEXP rnodes, rparents, rcatlist, rproblist;

	PROTECT(cnet);
	//if(isS4(cnet))
	//	printf("cnet is an object\n");

	rnodes = GET_SLOT(cnet, install("nodes"));
	//if(IS_VECTOR(m_nodeNames))
	//	printf("m_nodeNames is a vector\n");
	rparents = GET_SLOT(cnet, install("parents"));
	rcatlist = GET_SLOT(cnet, install("categories"));
	rproblist = GET_SLOT(cnet, install("probabilities"));

	if(rnodes == R_NilValue || rparents == R_NilValue || rcatlist == R_NilValue || rproblist == R_NilValue) {
		UNPROTECT(1);
		return R_NilValue;
	}

	SEXP res = show_catnet(rnodes, rparents, rcatlist, rproblist);

	UNPROTECT(1);

	return res;
}

SEXP catnetOptimalNetsForOrder(SEXP rSamples, SEXP rPerturbations, 
                              SEXP rMaxParents, SEXP rMaxComplexity, SEXP rOrder, 
                              SEXP rParentsPool, SEXP rFixedParentsPool, SEXP rUseCache, SEXP rEcho) {

  //clock_t lt = clock();
  //printf("Call estimateCatnets %ld\n", lt);

	//if(!isMatrix(rSamples))
	//	error("Data is not a matrix");
	if(!isNull(rPerturbations) && !isMatrix(rPerturbations))
		error("Perturbations should be a matrix");
	if(!isInteger(AS_INTEGER(rMaxParents)))
		error("maxParents should be an integer");
	if(!isInteger(AS_INTEGER(rMaxComplexity)))
		error("maxComplexity should be an integer");
	if(!isVector(rOrder))
		error("Order should be a vector");
	if(!isNull(rParentsPool) && !isVector(rParentsPool))
		error("ParentsPool should be a list");
	if(!isNull(rFixedParentsPool) && !isVector(rFixedParentsPool))
		error("FixedParentsPool should be a list");
	if(!isNull(rUseCache) && !isLogical(rUseCache))
		error("UseCache should be logical");
	if(!isNull(rEcho) && !isLogical(rEcho))
		error("Echo should be logical");

	RCatnetSearch * pengine = new RCatnetSearch;
	SEXP res = pengine -> estimateCatnets(rSamples, rPerturbations, 
						rMaxParents, rMaxComplexity, rOrder, 
	                                        rParentsPool, rFixedParentsPool, rUseCache, rEcho);
	delete pengine;

	//lt = clock();
	//printf("Exit estimateCatnets %ld\n", lt);
	//char str[128];
	//sprintf(str, "Mem Balance  %d\n", (int)g_memcounter);
	//printf(str);

	return res;
}

SEXP catnetOptimalClass2NetsForOrder(SEXP rSamples1, SEXP rPerturbations1,
				SEXP rSamples2, SEXP rPerturbations2,
				SEXP rMaxParents, SEXP rMaxComplexity, SEXP rOrder,
				SEXP rParentsPool, SEXP rFixedParentsPool, SEXP rUseCache, SEXP rEcho) {

	//if(!isMatrix(rSamples))
	//	error("Data is not a matrix");
	if(!isNull(rPerturbations1) && !isMatrix(rPerturbations1))
		error("Perturbations should be a matrix");
	if(!isNull(rPerturbations2) && !isMatrix(rPerturbations2))
		error("Perturbations should be a matrix");
	if(!isInteger(AS_INTEGER(rMaxParents)))
		error("maxParents should be an integer");
	if(!isInteger(AS_INTEGER(rMaxComplexity)))
		error("maxComplexity should be an integer");
	if(!isVector(rOrder))
		error("Order should be a vector");
	if(!isNull(rParentsPool) && !isVector(rParentsPool))
		error("ParentsPool should be a list");
	if(!isNull(rFixedParentsPool) && !isVector(rFixedParentsPool))
		error("FixedParentsPool should be a list");
	if(!isNull(rUseCache) && !isLogical(rUseCache))
		error("UseCache should be logical");
	if(!isNull(rEcho) && !isLogical(rEcho))
		error("Echo should be logical");

	RCatnetClass2Search * pengine = new RCatnetClass2Search;
	SEXP res = pengine -> estimateCatnets(rSamples1, rPerturbations1, 
						rSamples2, rPerturbations2, 
						rMaxParents, rMaxComplexity, rOrder, 
	                                        rParentsPool, rFixedParentsPool, rUseCache, rEcho);
	delete pengine;

	//char str[128];
	//sprintf(str, "Mem Balance  %d\n", (int)g_memcounter);
	//printf(str);

	return res;
}

SEXP catnetSetProb(SEXP cnet, SEXP rSamples, SEXP rPerturbations) {

	int *pSamples, *pSubsamples, *pPerturbations;
	int nnode, numsamples, numnodes, numsubsamples, j;
	SEXP dim;
	SEXP pcnet = R_NilValue;

	if(!isMatrix(rSamples))
		error("Data is not a matrix");
	if(!isNull(rPerturbations) && !isMatrix(rPerturbations))
		error("Perturbations is not a vector");

	PROTECT(cnet);
	RCatnet *rnet = new RCatnet(cnet);
	UNPROTECT(1);

	PROTECT(rSamples = AS_INTEGER(rSamples));
	pSamples = INTEGER(rSamples);

	dim = GET_DIM(rSamples);
	numnodes = INTEGER(dim)[0];
	numsamples = INTEGER(dim)[1];

	////////////////////////////////////////
	// Danger Ahead
	// We don's check that sample nodes actually correspond to the cnet's nodes
	// Missmatch of categories possible

	// pSamples are assumed positive indices
	for(j = 0; j < numnodes*numsamples; j++) {
		pSamples[j]--;
	}

	pSubsamples = 0;
	pPerturbations = 0;
	if(!isNull(rPerturbations)) {
		PROTECT(rPerturbations = AS_INTEGER(rPerturbations));
		pPerturbations = INTEGER(rPerturbations);
		pSubsamples = (int*)CATNET_MALLOC(numnodes*numsamples*sizeof(int));
	}

	for(nnode = 0; nnode < numnodes; nnode++) {
		if(pPerturbations) {
			numsubsamples = 0;
			for(j = 0; j < numsamples; j++) {		
				if(!pPerturbations[j * numnodes + nnode]) {
					memcpy(pSubsamples + numsubsamples*numnodes, 
						pSamples + j*numnodes, numnodes*sizeof(int));
					numsubsamples++;
				}
			}
			rnet->setNodeSampleProb(nnode, pSubsamples, numsubsamples, 1);
		}
		else {
			rnet->setNodeSampleProb(nnode, pSamples, numsamples, 1);
		}
	}

	UNPROTECT(1);
	if(pPerturbations) {
		UNPROTECT(1);
	}

	if(pSubsamples)
		CATNET_FREE(pSubsamples);
	pcnet = rnet->genRcatnet((const char*)"catNetwork");
	delete rnet;

	//char str[128];
	//sprintf(str, "Mem Balance  %d\n", (int)g_memcounter);
	//printf(str);

	return pcnet;

}

SEXP catnetLoglik(SEXP cnet, SEXP rSamples, SEXP rPerturbations) {

	int *pSamples, *pPerturbations;
	int numsamples, numnodes, j;
	double *floglik, *pvec;
	SEXP dim, rvec = R_NilValue;

	if(!isMatrix(rSamples))
		error("Data is not a matrix");
	if(!isNull(rPerturbations) && !isMatrix(rPerturbations))
		error("Perturbations is not a vector");

	////////////////////////////////////////
	// Danger Ahead
	// We don's check that sample nodes actually correspond to the cnet's nodes
	// Missmatch of categories possible

	PROTECT(cnet);
	RCatnet *rnet = new RCatnet(cnet);
	UNPROTECT(1);

	PROTECT(rSamples = AS_INTEGER(rSamples));
	pSamples = INTEGER(rSamples);

	dim = GET_DIM(rSamples);
	numnodes = INTEGER(dim)[0];
	numsamples = INTEGER(dim)[1];

	// pSamples are assumed positive indices
	for(j = 0; j < numnodes*numsamples; j++) {
		pSamples[j]--;
	}

	pPerturbations = 0;
	if(!isNull(rPerturbations)) {
		PROTECT(rPerturbations = AS_INTEGER(rPerturbations));
		pPerturbations = INTEGER(rPerturbations);
		UNPROTECT(1);
	}

	floglik = rnet->sampleLoglikVector(pSamples, numsamples);

	UNPROTECT(1);

	delete rnet;

	
	if(floglik) {
		PROTECT(rvec = NEW_NUMERIC(numsamples));
		pvec = NUMERIC_POINTER(rvec);
		memcpy(pvec, floglik, numsamples*sizeof(double));
		UNPROTECT(1);
		CATNET_FREE(floglik);
	}

	//char str[128];
	//sprintf(str, "Mem Balance  %d\n", (int)g_memcounter);
	//printf(str);

	return rvec;

}

} // extern "C"
