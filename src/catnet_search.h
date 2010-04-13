/*
 * catnet_search.h
 *
 *  Created on: Sep 25, 2009
 *      Author: nbalov
 */

#include "utils.h"
#include "catnet_class.h"

#ifndef CATNET_SEARCH_H
#define CATNET_SEARCH_H

template<class t_prob>
struct CATNET_CACHE_EL {
	int		nnode;
	PROB_LIST<t_prob> *pNodeProb;
	int		npars, *pPars;
	int		nPool, *pPool;
	t_prob		fLogLik;

	CATNET_CACHE_EL(int *ppool, int npool, int node, int *parset, int parsetsize, PROB_LIST<t_prob> *probNode, t_prob flik) {
		nPool = npool;
		pPool = (int*)CATNET_MALLOC(npool*sizeof(int));
		memcpy(pPool, ppool, npool*sizeof(int));
		nnode = node;
		npars = parsetsize;
		pPars = (int*)CATNET_MALLOC(npars*sizeof(int));
		memcpy(pPars, parset, npars*sizeof(int));
		pNodeProb = new PROB_LIST<t_prob>;
		*pNodeProb = *probNode;
		fLogLik = flik;
	}

	~CATNET_CACHE_EL() {
		if(pPool)
			CATNET_FREE(pPool);
		if(pPars)
			CATNET_FREE(pPars);
		if(pNodeProb)
			delete pNodeProb;
	}
};

template<class t_node, int t_node_size, class t_prob>
class CATNET_SEARCH {

protected:
	int m_nCatnets;
	CATNET<t_node, t_node_size, t_prob> **m_pCatnets;

	int m_numNodes, *m_pNodeNumCats, **m_pNodeCats, m_numSamples;

protected:
	/* returns increasing subsets of 'parset' of size 'parsize' */
	void combinationSets(int **&plist, int &nlist, int *curset, int *parset, int nparset, int parid, int parsize) {
		int i, ancestor;
 
		if(parid < 0 || parid >= parsize)
			return;
		ancestor = -1;
		if(parid > 0)
			ancestor = curset[parid-1];

		if(parid == parsize - 1) {
			for(i = 0; i < nparset; i++) {
				if(parset[i] <= ancestor)
					continue;
				int **pnewlist = (int**)CATNET_MALLOC((nlist+1)*sizeof(int*));
				if(nlist > 0)
					memcpy(pnewlist, plist, nlist*sizeof(int*));
				pnewlist[nlist] = (int*)CATNET_MALLOC(parsize*sizeof(int));
				if(curset) {
					memcpy(pnewlist[nlist], curset, (parsize-1)*sizeof(int));
				}
				pnewlist[nlist][parsize-1] = parset[i];

				CATNET_FREE(plist);
				plist = pnewlist;
				nlist++;
			}
			if(curset) {
				CATNET_FREE(curset);
				curset = 0;
			}
			return;
		}
		for(i = 0; i < nparset; i++) {
			if(parset[i] <= ancestor)
				continue;
			int *pnewset = (int*)CATNET_MALLOC((parid+1)*sizeof(int));
			if(curset && parid > 0)
				memcpy(pnewset, curset, parid*sizeof(int));
			pnewset[parid] = parset[i];
			combinationSets(plist, nlist, pnewset, parset, nparset, parid+1, parsize);
		}
		if(curset) {
			CATNET_FREE(curset);
			curset = 0;
		}
	}

public:
	CATNET_SEARCH() {
		m_nCatnets = 0;
		m_pCatnets = 0;
		m_numNodes = 0;
		m_pNodeNumCats = 0;
		m_pNodeCats = 0;
	}

	~CATNET_SEARCH() {
		int i;
		if(m_pCatnets) {
			for(i = 0; i < m_nCatnets; i++)
				if(m_pCatnets[i]) {
					delete m_pCatnets[i];
					m_pCatnets[i] = 0;
				}
			CATNET_FREE(m_pCatnets);
		}
		m_pCatnets = 0;
		m_nCatnets = 0;
		if(m_pNodeCats) {
			for(i = 0; i < m_numNodes; i++) 
				if(m_pNodeCats[i])
					CATNET_FREE(m_pNodeCats[i]);
			CATNET_FREE(m_pNodeCats);
		}
		if(m_pNodeNumCats) 
			CATNET_FREE(m_pNodeNumCats);
	}

	virtual int getCachedProb(int *ppool, int poolsize, int node, int *parset, int parsize, 
							PROB_LIST<t_prob> *probNode, t_prob *pflik) {
		return 0;
	}

	virtual int setCachedProb(int *ppool, int poolsize, int node, int *parset, int parsize, 
							PROB_LIST<t_prob> *probNode, t_prob flik) {
		return 0;
	}

	/* psamples and perturbations are sample=columns and node=rows. */
	/* Each parentsPool[i] is numnodes long ! */
	int estimate(int numnodes, int numsamples, int *psamples, int *perturbations, 
               int maxParentSet, int maxComplexity,    
               int **parentsPool, int **fixedParentsPool, int becho) {

		int i, j, k, d, ncomb, ncombMaxLogLik, nnode;
		int maxCategories, numsubsamples, complx;
		int *parset, parsetsize, *idparset, *fixparset, fixparsetsize;
		int *paux, *psubsamples, **pcomblist, ncomblist, maxpars, ballow, bfixallow;
		CATNET<t_node, t_node_size, t_prob> baseCatnet, *pNewNet, **pCurCatnetList;

		t_prob fLogLik, fMaxLogLik, loglik;
		PROB_LIST<t_prob> probMaxNode, *pProbNode;

		if(numnodes < 1 || numsamples < 1 || !psamples)
			return 0;
		if(maxComplexity < numnodes)
			maxComplexity = numnodes;

		m_numNodes = numnodes;
		m_numSamples = numsamples;

		maxCategories = 0;

		m_pNodeNumCats = (int*)CATNET_MALLOC(numnodes*sizeof(int));
		m_pNodeCats = (int**)CATNET_MALLOC(numnodes*sizeof(int*));
		memset(m_pNodeCats, 0, numnodes*sizeof(int*));
		memset(m_pNodeNumCats, 0, numnodes*sizeof(int));

		for(j = 0; j < numsamples; j++) {
			for(i = 0; i < numnodes; i++) {
				/* convert sample categories to indices */
				d = -1;
				if(m_pNodeCats[i]) {
					for(k = 0; k < m_pNodeNumCats[i]; k++)
						if(m_pNodeCats[i][k] == psamples[j*numnodes + i])
							d = k;
				}

				if(d == -1) {
					paux = (int*)CATNET_MALLOC((m_pNodeNumCats[i]+1)*sizeof(int));
					if(m_pNodeCats[i]) {
						memcpy(paux, m_pNodeCats[i], m_pNodeNumCats[i]*sizeof(int));
						CATNET_FREE(m_pNodeCats[i]);
						m_pNodeCats[i] = paux;
					}
					else
						m_pNodeCats[i] = paux;
					d = m_pNodeNumCats[i];
					m_pNodeCats[i][d] = psamples[j*numnodes + i];
					m_pNodeNumCats[i]++;
				}
			}
		}
		for(i = 0; i < numnodes; i++) {
			/* order m_pNodeNumCats[i] */
			for(j = 0; j < m_pNodeNumCats[i]; j++) {
				for(k = j + 1; k < m_pNodeNumCats[i]; k++) {
					if(m_pNodeCats[i][j] > m_pNodeCats[i][k]) {
						d = m_pNodeCats[i][j]; 
						m_pNodeCats[i][j] = m_pNodeCats[i][k];
						m_pNodeCats[i][k] = d;
					}
				}
			} 
			for(j = 0; j < numsamples; j++) {
				for(d = 0; d < m_pNodeNumCats[i]; d++)
					if(m_pNodeCats[i][d] == psamples[j*numnodes + i])
						break;
				psamples[j*numnodes + i] = d;
			}
			if(maxCategories < m_pNodeNumCats[i])
				maxCategories = m_pNodeNumCats[i];
		}


		parset = (int*)CATNET_MALLOC(numnodes*sizeof(int));
		idparset = (int*)CATNET_MALLOC(numnodes*sizeof(int));
		fixparset = (int*)CATNET_MALLOC(numnodes*sizeof(int));

		m_nCatnets = maxComplexity + 1;
		m_pCatnets = (CATNET<t_node, t_node_size, t_prob>**)CATNET_MALLOC(m_nCatnets*sizeof(CATNET<t_node, t_node_size, t_prob>*));
		memset(m_pCatnets, 0, m_nCatnets*sizeof(CATNET<t_node, t_node_size, t_prob>*));

		pCurCatnetList = (CATNET<t_node, t_node_size, t_prob>**)CATNET_MALLOC(m_nCatnets*sizeof(CATNET<t_node, t_node_size, t_prob>*));

		psubsamples = 0;
		if(perturbations) {
			psubsamples = (int*)CATNET_MALLOC(numnodes*numsamples*sizeof(int));
		}

		/* create a network without edges*/
		pNewNet = new CATNET<t_node, t_node_size, t_prob>
				(numnodes, 0/*maxParentSet*/, maxCategories, 0, 0, 0, m_pNodeNumCats);

		/* set parents */
		for(nnode = 0; nnode < numnodes; nnode++) {
			if(fixedParentsPool && fixedParentsPool[nnode]) {
				fixparsetsize = 0;
				for(j = 0; j < nnode; j++) {
					ballow = 1;
					if(parentsPool && parentsPool[nnode]) {
						ballow = 0;
						for(k = 0; k < numnodes; k++) {
							if(j == parentsPool[nnode][k])
								ballow = 1;
						}
					}
					if(parentsPool && !parentsPool[nnode])
						ballow = 0;
					bfixallow = 0;
					if(fixedParentsPool && fixedParentsPool[nnode]) {
						for(k = 0; k < numnodes; k++) {
							if(j == fixedParentsPool[nnode][k])
								bfixallow = 1;
						}
					}
					if(!ballow)
						continue;
					if(bfixallow) {
					  fixparset[fixparsetsize] = j;
					  fixparsetsize++;
					}
				}
				if(fixparsetsize > 0)
					pNewNet -> setParents(nnode, fixparset, fixparsetsize);
			}
		}

		baseCatnet.init(numnodes, maxParentSet, maxCategories, 0, 0, 0, m_pNodeNumCats);

		// set sample probabilities and calculate log-likelihood
		for(nnode = 0; nnode < numnodes; nnode++) {
			if(perturbations) {
				numsubsamples = 0;
				for(j = 0; j < numsamples; j++) {
					if(!perturbations[j * numnodes + nnode]) {
						memcpy(psubsamples + numsubsamples*numnodes, psamples + j*numnodes, numnodes*sizeof(int));
						numsubsamples++;
					}
				}
				pNewNet->setNodeSampleProb(nnode, psubsamples, numsubsamples);
			}
			else {
				pNewNet->setNodeSampleProb(nnode, psamples, numsamples);
			}
		}
		loglik = pNewNet->loglik();
		complx = pNewNet->complexity();
		m_pCatnets[complx] = pNewNet;

		/* main loop of consequential non-empty-parenthood-node additions */
		for(nnode = 1; nnode < numnodes; nnode++) {

			if(becho) {
				printf("processing node %d\n", nnode+1);
				printf("    [#parents][#combinations] = ");
			}

			fixparsetsize = 0;
			parsetsize = 0;
			for(j = 0; j < nnode; j++) {
				ballow = 1;
				if(parentsPool && parentsPool[nnode]) {
					ballow = 0;
					for(k = 0; k < numnodes; k++) {
						if(j == parentsPool[nnode][k])
							ballow = 1;
					}
				}
				if(parentsPool && !parentsPool[nnode])
					ballow = 0;
				bfixallow = 0;
				if(fixedParentsPool && fixedParentsPool[nnode]) {
					for(k = 0; k < numnodes; k++) {
						if(j == fixedParentsPool[nnode][k])
							bfixallow = 1;
					}
				}
				if(!ballow)
					continue;
				if(bfixallow) {
				  fixparset[fixparsetsize] = j;
				  fixparsetsize++;
				}
				else {
				  parset[parsetsize] = j;
				  parsetsize++;
				}
			}
			/* extend the content before sending to cache; parsetsize + fixparsetsize < numnodes */
			memcpy(parset + parsetsize, fixparset, fixparsetsize*sizeof(int));

			maxpars = maxParentSet;
			if(maxpars > parsetsize + fixparsetsize)
				maxpars = parsetsize + fixparsetsize;

			memset(pCurCatnetList, 0, m_nCatnets*sizeof(CATNET<t_node, t_node_size, t_prob>*));

			for(d = fixparsetsize + 1; d <= maxpars; d++) {
				if(!getCachedProb(parset, parsetsize + fixparsetsize, nnode, idparset, d, &probMaxNode, &fMaxLogLik)) { 

					pcomblist = 0;
					ncomblist = 0;
					combinationSets(pcomblist, ncomblist, 0, parset, parsetsize, 0, d - fixparsetsize);

				        if(fixparsetsize > 0) {
				        	if(!pcomblist || ncomblist < 1) {
				        	    	pcomblist = (int**)CATNET_MALLOC(1*sizeof(int*));
				            		pcomblist[0] = 0;	
				            		ncomblist = 1;
				          	}
				        	for(k = 0; k < ncomblist; k++) {
				            		paux = (int*)CATNET_MALLOC(d*sizeof(int));
							for(j = 0; j < fixparsetsize; j++)
				            			paux[j] = fixparset[j];
					        	if(pcomblist[k] && d > fixparsetsize) {
				            			memcpy(paux + fixparsetsize, pcomblist[k], (d-fixparsetsize)*sizeof(int));
							}
				            		if(pcomblist[k])
				            			CATNET_FREE(pcomblist[k]); 
				           		pcomblist[k] = paux;
						}
					}
			
					if(becho)
						printf("[%d]%d  ", d, ncomblist);

					fMaxLogLik = -FLT_MAX;
					ncombMaxLogLik = -1;
					probMaxNode.reset();

					for(ncomb = 0; ncomb < ncomblist; ncomb++) {
						// add pcomplist[j] parent set to nnode
						baseCatnet.setParents(nnode, pcomblist[ncomb], d);
			     
						// add perturbation
						if(perturbations) {
							numsubsamples = 0;
							for(j = 0; j < numsamples; j++) {
								if(!perturbations[j * numnodes + nnode]) {
									memcpy(psubsamples + numsubsamples*numnodes, psamples + j*numnodes, numnodes*sizeof(int));
									numsubsamples++;
								}
							}
							fLogLik = baseCatnet.setNodeSampleProb(nnode, psubsamples, numsubsamples);
						}
						else {
							fLogLik = baseCatnet.setNodeSampleProb(nnode, psamples, numsamples);
						}

						if(fMaxLogLik < fLogLik) {
							pProbNode = baseCatnet.getNodeProb(nnode); 
							fMaxLogLik = fLogLik;
							ncombMaxLogLik = ncomb;
							if(pProbNode)
								probMaxNode = *pProbNode;
						}
					} /* for ncomb */

					if(ncombMaxLogLik < 0){
						/* retain the same m_pCatnets list */
						continue;
					}

					memcpy(idparset, pcomblist[ncombMaxLogLik], d*sizeof(int));
					setCachedProb(parset, parsetsize + fixparsetsize, 
								nnode, idparset, d, &probMaxNode, fMaxLogLik);

					/* release combination set */
        				for(ncomb = 0; ncomb < ncomblist; ncomb++) {
        	  				if(pcomblist[ncomb])
        	    					CATNET_FREE(pcomblist[ncomb]);
        	  				pcomblist[ncomb] = NULL;
					} /* for ncomb */
        				CATNET_FREE(pcomblist);
        				pcomblist = 0;
					ncomblist = 0;
	
				} /* if(!getCachedProb) */

				for(k = 0; k < m_nCatnets; k++) {
					if(!m_pCatnets[k])
						continue;
					baseCatnet = *m_pCatnets[k];
					baseCatnet.setParents(nnode, idparset, d);
					baseCatnet.setNodeProb(nnode, &probMaxNode);
					complx = baseCatnet.complexity();
					if(complx > maxComplexity) {
						continue;
					}
					loglik = baseCatnet.loglik();
					if(pCurCatnetList[complx]) {
						if(pCurCatnetList[complx]->getLoglik() < loglik) {
							*pCurCatnetList[complx] = baseCatnet;
						}
					}
					else {
						pCurCatnetList[complx] = new CATNET<t_node, t_node_size, t_prob>;
						*pCurCatnetList[complx] = baseCatnet;
					}
				} /* for k */

			} /* for d */

			if(becho)
				printf("\n");

			for(j = 0; j < m_nCatnets; j++) {
				if(m_pCatnets[j]) {
					if(pCurCatnetList[j]) {
						if(m_pCatnets[j]->getLoglik() < pCurCatnetList[j]->getLoglik()) {
							delete m_pCatnets[j];
							m_pCatnets[j] = pCurCatnetList[j];
							pCurCatnetList[j] = 0;
						}
						else {
							delete pCurCatnetList[j];
							pCurCatnetList[j] = 0;
						}
					}
				}
				else {
					m_pCatnets[j] = pCurCatnetList[j];
					pCurCatnetList[j] = 0;
				}
			}
		} // for(nnode = 0; nnode < numnodes; nnode++)

		CATNET_FREE(pCurCatnetList);
		CATNET_FREE(parset);
		CATNET_FREE(fixparset);
		CATNET_FREE(idparset);

		if(psubsamples)
			CATNET_FREE(psubsamples);

		if(m_pNodeCats) {
			for(i = 0; i < m_numNodes; i++) 
				if(m_pNodeCats[i])
					CATNET_FREE(m_pNodeCats[i]);
			CATNET_FREE(m_pNodeCats);
		}
		m_pNodeCats = 0;
		if(m_pNodeNumCats) 
			CATNET_FREE(m_pNodeNumCats);
		m_pNodeNumCats = 0;

		for(j = 0; j < m_nCatnets; j++) {
			if(m_pCatnets[j]) {
				m_pCatnets[j] -> normalizeProbabilities();
			}
		} 
		return 0;
	}

};

#endif /* CATNET_SEARCH_H */
