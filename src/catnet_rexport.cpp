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

extern "C" {

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

	delete rnet;
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

char *gen_prob_string(int node, SEXP parlist, int paridx, SEXP catlist, SEXP problist, char *str) {
	int j, npar;
	SEXP parprobs, pcats;
	char *newstr, *aux, *aux2, *aux3;

	//char ss[128];

	if(!str) {
		str = (char*)CATNET_MALLOC(1);
		str[0] = 0;
	}

	if(paridx >= length(parlist)) {
		pcats = VECTOR_ELT(catlist, node);
		newstr = (char*)CATNET_MALLOC(((strlen(str)+1+32)*length(pcats))*sizeof(char));

		newstr[0] = 0;
		for(j = 0; j < length(pcats); j++) {
			sprintf(newstr, "%s%s%s %f\n", newstr, str, CHAR(STRING_ELT(pcats, j)), NUMERIC_POINTER(problist)[j]);
		}

		CATNET_FREE(str);
		str = newstr;
		return str;
	}

	npar = INTEGER_POINTER(parlist)[paridx] - 1;
	pcats = VECTOR_ELT(catlist, npar);

	newstr = (char*)CATNET_MALLOC(sizeof(char));
	newstr[0] = 0;
	for(j = 0; j < length(pcats); j++) {
		parprobs = VECTOR_ELT(problist, j);

		aux = (char*)CATNET_MALLOC((strlen(str)+1+8)*sizeof(char));

		//sprintf(aux, "%d, %d, %d,  %s\n", node, npar, length(parprobs), CHAR(STRING_ELT(pcats, j)));
		//printf(aux);

		sprintf(aux, "%s%s", str, CHAR(STRING_ELT(pcats, j)));
		aux2 = gen_prob_string(node, parlist, paridx + 1, catlist, parprobs, aux);

		aux3 = (char*)CATNET_MALLOC((strlen(newstr)+strlen(aux2)+2)*sizeof(char));
		sprintf(aux3, "%s%s", newstr, aux2);
		//sprintf(ss, "5 CATNET_FREE: %p, %d\n", newstr, strlen(newstr)); printf(ss);
		CATNET_FREE(newstr);
		newstr = aux3;

		CATNET_FREE(aux2);
	}
	CATNET_FREE(str);
	str = newstr;

	return str;
}

SEXP prob_string(SEXP rnodes, SEXP rparents, SEXP rcatlist, SEXP rproblist) {

	int node;
	SEXP nodepars, nodeproblist, pstr;
	char *str = NULL, *newstr, *aux;

	PROTECT(rnodes = AS_LIST(rnodes));
	PROTECT(rparents = AS_LIST(rparents));
	PROTECT(rcatlist = AS_LIST(rcatlist));
	PROTECT(rproblist = AS_LIST(rproblist));

	for(node = 0; node < length(rnodes); node++) {
		nodepars = VECTOR_ELT(rparents, node);
		nodeproblist = VECTOR_ELT(rproblist, node);
		newstr = gen_prob_string(node, nodepars, 0, rcatlist, nodeproblist, NULL);
		if(str) {
			aux = (char*)CATNET_MALLOC((strlen(str)+strlen(newstr)+1+16)*sizeof(char));
			sprintf(aux, "%sNode [%d]:\n%s", str, node, newstr);
			CATNET_FREE(str);
			CATNET_FREE(newstr);
			str = aux;
		}
		else
			str = newstr;
	}

	PROTECT(pstr = allocVector(STRSXP, 1));
	SET_STRING_ELT(pstr, 0, mkChar(str));

	UNPROTECT(5);
	return pstr;
}


void gen_prob_vector(int node, SEXP parlist, int paridx, SEXP catlist, SEXP problist, double *&pvec, int &nvec) {
	int j, npar;
	SEXP parprobs, pcats;
	double *newvec;

	if(!pvec) {
		pvec = (double*)CATNET_MALLOC(sizeof(double));
		nvec = 0;
	}

	if(paridx >= length(parlist)) {
		pcats = VECTOR_ELT(catlist, node);
		if (length(problist) != length(pcats)) {
			printf("gen_prob_vector: length(problist) != length(pcats))\n");
			return;
		}
		newvec = (double*)CATNET_MALLOC((nvec + length(pcats))*sizeof(double));
		memcpy(newvec, pvec, nvec*sizeof(double));
		for(j = 0; j < length(pcats); j++) {
			newvec[nvec+j] = NUMERIC_POINTER(problist)[j];
		}
		CATNET_FREE(pvec);
		pvec = newvec;
		nvec += length(pcats);
		return;
	}

	npar = INTEGER_POINTER(parlist)[paridx] - 1;
	pcats = VECTOR_ELT(catlist, npar);
	if (length(problist) != length(pcats)) {
		printf("gen_prob_vector: length(problist) != length(pcats))\n");
		return;
	}
	for(j = 0; j < length(pcats); j++) {
		parprobs = VECTOR_ELT(problist, j);
		gen_prob_vector(node, parlist, paridx + 1, catlist, parprobs, pvec, nvec);
	}
}

SEXP prob_vector(SEXP rnodes, SEXP rparents, SEXP rcatlist, SEXP rproblist) {

	int node, nvec;
	SEXP nodepars, nodeproblist, pstr;
	SEXP rvec;
	double *pvec, *prvec;

	PROTECT(rnodes = AS_LIST(rnodes));
	PROTECT(rparents = AS_LIST(rparents));
	PROTECT(rcatlist = AS_LIST(rcatlist));
	PROTECT(rproblist = AS_LIST(rproblist));

	PROTECT(pstr = allocVector(VECSXP, length(rnodes)));

	for(node = 0; node < length(rnodes); node++) {
		nodepars = VECTOR_ELT(rparents, node);
		nodeproblist = VECTOR_ELT(rproblist, node);
		pvec = 0;
		nvec = 0;
		gen_prob_vector(node, nodepars, 0, rcatlist, nodeproblist, pvec, nvec);
		PROTECT(rvec = NEW_NUMERIC(nvec));
		prvec = NUMERIC_POINTER(rvec);
		memcpy(prvec, pvec, nvec*sizeof(double));
		CATNET_FREE(pvec);
		SET_VECTOR_ELT(pstr, node, rvec);
		UNPROTECT(1);
	}

	UNPROTECT(5);
	return pstr;
}

SEXP catnetOptimalNetsForOrder(SEXP rSamples, SEXP rPerturbations, 
                              SEXP rMaxParents, SEXP rMaxComplexity, SEXP rOrder, 
                              SEXP rParentsPool, SEXP rFixedParentsPool, SEXP rUseCache) {

  //clock_t lt = clock();
  //printf("Call estimateCatnets %ld\n", lt);

	//if(!isMatrix(rSamples))
	//	error("Data is not a matrix");
	if(!isNull(rPerturbations) && !isMatrix(rPerturbations))
		error("Perturbations is not a vector");
	if(!isInteger(AS_INTEGER(rMaxParents)))
		error("MaxParents is not an integer");
	if(!isInteger(AS_INTEGER(rMaxComplexity)))
		error("rMaxComplexity is not an integer");
	if(!isVector(rOrder))
		error("rOrder is not a vector");
	if(!isNull(rParentsPool) && !isVector(rParentsPool))
		error("ParentsPool is not a list");
	if(!isNull(rFixedParentsPool) && !isVector(rFixedParentsPool))
		error("FixedParentsPool is not a list");
	if(!isNull(rUseCache) && !isLogical(rUseCache))
		error("rUseCache is not a boolean");

	RCatnetSearch * pengine = new RCatnetSearch;
	SEXP res = pengine -> estimateCatnets(rSamples, rPerturbations, 
						rMaxParents, rMaxComplexity, rOrder, 
	                                        rParentsPool, rFixedParentsPool, rUseCache);
	delete pengine;

  //lt = clock();
  //printf("Exit estimateCatnets %ld\n", lt);

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
	return pcnet;

}

SEXP catnetLoglik(SEXP cnet, SEXP rSamples, SEXP rPerturbations) {

	int *pSamples, *pPerturbations;
	int numsamples, numnodes, j;
	double floglik, *pvec;
	SEXP dim, rvec;

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

	floglik = rnet->sampleLoglik(pSamples, numsamples);


	UNPROTECT(1);

	delete rnet;

	PROTECT(rvec = NEW_NUMERIC(1));
	pvec = NUMERIC_POINTER(rvec);
	pvec[0] = floglik;
	UNPROTECT(1);

	return rvec;

}

} // extern "C"
