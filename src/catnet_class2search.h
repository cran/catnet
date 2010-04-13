/*
 * catnet_class2search.h
 *
 *  Created on: Aprt 6, 2010
 *      Author: nbalov
 */

#include "utils.h"
#include "catnet_class.h"

#ifndef CATNET_CLASS2SEARCH_H
#define CATNET_CLASS2SEARCH_H

template<class t_node, int t_node_size, class t_prob>
class CATNET_CLASS2SEARCH {

protected:
	int m_nCatnets;
	t_prob *m_pCatnetScores;
	CATNET<t_node, t_node_size, t_prob> **m_pCatnets1, **m_pCatnets2;

	int m_numNodes, *m_pNodeNumCats, **m_pNodeCats, m_numSamples1, m_numSamples2;

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
	CATNET_CLASS2SEARCH() {
		m_nCatnets = 0;
		m_pCatnets1 = 0;
		m_pCatnets2 = 0;
		m_pCatnetScores = 0;
		m_numNodes = 0;
		m_pNodeNumCats = 0;
		m_pNodeCats = 0;
	}

	~CATNET_CLASS2SEARCH() {
		int i;
		if(m_pCatnets1) {
			for(i = 0; i < m_nCatnets; i++)
				if(m_pCatnets1[i]) {
					delete m_pCatnets1[i];
					m_pCatnets1[i] = 0;
				}
			CATNET_FREE(m_pCatnets1);
		}
		m_pCatnets1 = 0;
		if(m_pCatnets2) {
			for(i = 0; i < m_nCatnets; i++)
				if(m_pCatnets2[i]) {
					delete m_pCatnets2[i];
					m_pCatnets2[i] = 0;
				}
			CATNET_FREE(m_pCatnets2);
		}
		m_pCatnets2 = 0;
		m_nCatnets = 0;
		if(m_pCatnetScores)
			CATNET_FREE(m_pCatnetScores);
		m_pCatnetScores = 0;
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
							PROB_LIST<t_prob> *probNode1, PROB_LIST<t_prob> *probNode2, t_prob *pflik) {
		return 0;
	}

	virtual int setCachedProb(int *ppool, int poolsize, int node, int *parset, int parsize, 
							PROB_LIST<t_prob> *probNode1, PROB_LIST<t_prob> *probNode2, t_prob flik) {
		return 0;
	}

	/* psamples and perturbations are sample=columns and node=rows. */
	/* Each parentsPool[i] is numnodes long ! */
	int estimate(int numnodes, 
		int numsamples1, int *psamples1, int *perturbations1, 
		int numsamples2, int *psamples2, int *perturbations2, 
		int maxParentSet, int maxComplexity,    
		int **parentsPool, int **fixedParentsPool, int becho) {

		int i, j, k, d, ncomb, ncombMaxLogLik, nnode;
		int maxCategories, complx;
		int *parset, parsetsize, *idparset, *fixparset, fixparsetsize;
		int *paux, **pcomblist, ncomblist, maxpars, ballow, bfixallow;
		CATNET<t_node, t_node_size, t_prob> baseCatnet1, baseCatnet2, *pNewNet1, *pNewNet2, **pCurCatnetList1, **pCurCatnetList2;

		int *psubsamples1, numsubsamples1, *psubsamples2, numsubsamples2; 

		t_prob fdiff, fLogLik1, fLogLik2, loglik1, loglik2, fMaxLogLik, *pCurCatnetScores;
		PROB_LIST<t_prob> probMaxNode1, probMaxNode2, *pProbNode;

		if(numnodes < 1 || numsamples1 < 1 || numsamples2 < 1 || !psamples1 || !psamples2)
			return 0;
		if(maxComplexity < numnodes)
			maxComplexity = numnodes;

		m_numNodes = numnodes;
		m_numSamples1 = numsamples1;
		m_numSamples2 = numsamples2;

		maxCategories = 0;

		if(m_pNodeNumCats)
			CATNET_FREE(m_pNodeNumCats);
		m_pNodeNumCats = (int*)CATNET_MALLOC(numnodes*sizeof(int));
		if(m_pNodeCats)
			CATNET_FREE(m_pNodeCats);
		m_pNodeCats = (int**)CATNET_MALLOC(numnodes*sizeof(int*));
		memset(m_pNodeCats, 0, numnodes*sizeof(int*));
		memset(m_pNodeNumCats, 0, numnodes*sizeof(int));

		for(j = 0; j < numsamples1; j++) {
			for(i = 0; i < numnodes; i++) {
				/* convert sample categories to indices */
				d = -1;
				if(m_pNodeCats[i]) {
					for(k = 0; k < m_pNodeNumCats[i]; k++)
						if(m_pNodeCats[i][k] == psamples1[j*numnodes + i])
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
					m_pNodeCats[i][d] = psamples1[j*numnodes + i];
					m_pNodeNumCats[i]++;
				}
			}
		}
		for(j = 0; j < numsamples2; j++) {
			for(i = 0; i < numnodes; i++) {
				/* convert sample categories to indices */
				d = -1;
				if(m_pNodeCats[i]) {
					for(k = 0; k < m_pNodeNumCats[i]; k++)
						if(m_pNodeCats[i][k] == psamples2[j*numnodes + i])
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
					m_pNodeCats[i][d] = psamples2[j*numnodes + i];
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
			for(j = 0; j < numsamples1; j++) {
				for(d = 0; d < m_pNodeNumCats[i]; d++)
					if(m_pNodeCats[i][d] == psamples1[j*numnodes + i])
						break;
				psamples1[j*numnodes + i] = d;
			}
			for(j = 0; j < numsamples2; j++) {
				for(d = 0; d < m_pNodeNumCats[i]; d++)
					if(m_pNodeCats[i][d] == psamples2[j*numnodes + i])
						break;
				psamples2[j*numnodes + i] = d;
			}
			if(maxCategories < m_pNodeNumCats[i])
				maxCategories = m_pNodeNumCats[i];
		}


		parset = (int*)CATNET_MALLOC(numnodes*sizeof(int));
		idparset = (int*)CATNET_MALLOC(numnodes*sizeof(int));
		fixparset = (int*)CATNET_MALLOC(numnodes*sizeof(int));

		if(m_pCatnets1) {
			for(i = 0; i < m_nCatnets; i++)
				if(m_pCatnets1[i]) {
					delete m_pCatnets1[i];
					m_pCatnets1[i] = 0;
				}
			CATNET_FREE(m_pCatnets1);
		}
		m_pCatnets1 = 0;
		if(m_pCatnets2) {
			for(i = 0; i < m_nCatnets; i++)
				if(m_pCatnets2[i]) {
					delete m_pCatnets2[i];
					m_pCatnets2[i] = 0;
				}
			CATNET_FREE(m_pCatnets2);
		}
		m_pCatnets2 = 0;
		m_nCatnets = maxComplexity + 1;
		m_pCatnets1 = (CATNET<t_node, t_node_size, t_prob>**)CATNET_MALLOC(m_nCatnets*sizeof(CATNET<t_node, t_node_size, t_prob>*));
		memset(m_pCatnets1, 0, m_nCatnets*sizeof(CATNET<t_node, t_node_size, t_prob>*));
		m_pCatnets2 = (CATNET<t_node, t_node_size, t_prob>**)CATNET_MALLOC(m_nCatnets*sizeof(CATNET<t_node, t_node_size, t_prob>*));
		memset(m_pCatnets2, 0, m_nCatnets*sizeof(CATNET<t_node, t_node_size, t_prob>*));

		pCurCatnetList1 = (CATNET<t_node, t_node_size, t_prob>**)CATNET_MALLOC(m_nCatnets*sizeof(CATNET<t_node, t_node_size, t_prob>*));
		pCurCatnetList2 = (CATNET<t_node, t_node_size, t_prob>**)CATNET_MALLOC(m_nCatnets*sizeof(CATNET<t_node, t_node_size, t_prob>*));

		if(m_pCatnetScores)
			CATNET_FREE(m_pCatnetScores);
		m_pCatnetScores = (t_prob*)CATNET_MALLOC(m_nCatnets*sizeof(t_prob));
		for(j = 0; j < m_nCatnets; j++)
			m_pCatnetScores[j] = -FLT_MAX;
		pCurCatnetScores = (t_prob*)CATNET_MALLOC(m_nCatnets*sizeof(t_prob));

		psubsamples1 = 0;
		if(perturbations1) {
			psubsamples1 = (int*)CATNET_MALLOC(numnodes*numsamples1*sizeof(int));
		}

		psubsamples2 = 0;
		if(perturbations2) {
			psubsamples2 = (int*)CATNET_MALLOC(numnodes*numsamples2*sizeof(int));
		}

		/* create a network without edges*/
		pNewNet1 = new CATNET<t_node, t_node_size, t_prob>
				(numnodes, 0/*maxParentSet*/, maxCategories, 0, 0, 0, m_pNodeNumCats);
		pNewNet2 = new CATNET<t_node, t_node_size, t_prob>
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
				if(fixparsetsize > 0) {
					pNewNet1 -> setParents(nnode, fixparset, fixparsetsize);
					pNewNet2 -> setParents(nnode, fixparset, fixparsetsize);
				}
			}
		}

		baseCatnet1.init(numnodes, maxParentSet, maxCategories, 0, 0, 0, m_pNodeNumCats);
		baseCatnet2.init(numnodes, maxParentSet, maxCategories, 0, 0, 0, m_pNodeNumCats);

		// set sample probabilities and calculate log-likelihood
		for(nnode = 0; nnode < numnodes; nnode++) {
			if(perturbations2) {
				numsubsamples2 = 0;
				for(j = 0; j < numsamples1; j++) {
					if(!perturbations2[j * numnodes + nnode]) {
						memcpy(psubsamples2 + numsubsamples2*numnodes, psamples1 + j*numnodes, numnodes*sizeof(int));
						numsubsamples2++;
					}
				}
				pNewNet2->setNodeSampleProb(nnode, psubsamples2, numsubsamples2);
			}
			else {
				pNewNet2->setNodeSampleProb(nnode, psamples2, numsamples2);
			}
		}
		for(nnode = 0; nnode < numnodes; nnode++) {
			if(perturbations1) {
				numsubsamples1 = 0;
				for(j = 0; j < numsamples1; j++) {
					if(!perturbations1[j * numnodes + nnode]) {
						memcpy(psubsamples1 + numsubsamples1*numnodes, psamples1 + j*numnodes, numnodes*sizeof(int));
						numsubsamples1++;
					}
				}
				pNewNet1->setNodeSampleProb(nnode, psubsamples1, numsubsamples1);
			}
			else {
				pNewNet1->setNodeSampleProb(nnode, psamples1, numsamples1);
			}
		}
		loglik1 = pNewNet1->loglik();
		loglik2 = pNewNet2->loglik();
		complx = pNewNet1->complexity();
		m_pCatnets1[complx] = pNewNet1;
		m_pCatnets2[complx] = pNewNet2;
		m_pCatnetScores[complx] = loglik1 - loglik2;

		//printf("pCurCatnetScores[%d] = %f\n", complx, pCurCatnetScores[complx]);

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

			memset(pCurCatnetList1, 0, m_nCatnets*sizeof(CATNET<t_node, t_node_size, t_prob>*));
			memset(pCurCatnetList2, 0, m_nCatnets*sizeof(CATNET<t_node, t_node_size, t_prob>*));
			for(j = 0; j < m_nCatnets; j++)
				pCurCatnetScores[j] = -FLT_MAX;

			for(d = fixparsetsize + 1; d <= maxpars; d++) {
				//if(!getCachedProb(parset, parsetsize + fixparsetsize, nnode, idparset, d, &probMaxNode, &fMaxLogLik)) { 

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
					probMaxNode1.reset();
					probMaxNode2.reset();

					for(ncomb = 0; ncomb < ncomblist; ncomb++) {
						// add pcomplist[j] parent set to nnode
						baseCatnet1.setParents(nnode, pcomblist[ncomb], d);
						baseCatnet2.setParents(nnode, pcomblist[ncomb], d);

						// add perturbation
						if(perturbations2) {
							numsubsamples2 = 0;
							for(j = 0; j < numsamples2; j++) {
								if(!perturbations2[j * numnodes + nnode]) {
									memcpy(psubsamples2 + numsubsamples2*numnodes, psamples2 + j*numnodes, numnodes*sizeof(int));
									numsubsamples2++;
								}
							}
							fLogLik2 = baseCatnet2.setNodeSampleProb(nnode, psubsamples2, numsubsamples2);
						}
						else {
							fLogLik2 = baseCatnet2.setNodeSampleProb(nnode, psamples2, numsamples2);
						}
						if(perturbations1) {
							numsubsamples1 = 0;
							for(j = 0; j < numsamples1; j++) {
								if(!perturbations1[j * numnodes + nnode]) {
									memcpy(psubsamples1 + numsubsamples1*numnodes, psamples1 + j*numnodes, 	numnodes*sizeof(int));
									numsubsamples1++;
								}
							}
							fLogLik1 = baseCatnet1.setNodeSampleProb(nnode, psubsamples1, numsubsamples1);
						}
						else {
							fLogLik1 = baseCatnet1.setNodeSampleProb(nnode, psamples1, numsamples1);
						}

						if(fMaxLogLik < fLogLik1 - fLogLik2) {
							fMaxLogLik = fLogLik1 - fLogLik2;
							ncombMaxLogLik = ncomb;
							pProbNode = baseCatnet1.getNodeProb(nnode); 
							if(pProbNode)
								probMaxNode1 = *pProbNode;
							pProbNode = baseCatnet2.getNodeProb(nnode); 
							if(pProbNode)
								probMaxNode2 = *pProbNode;
						}
					} /* for ncomb */

					if(ncombMaxLogLik < 0){
						/* retain the same m_pCatnets list */
						continue;
					}

					memcpy(idparset, pcomblist[ncombMaxLogLik], d*sizeof(int));
					//setCachedProb(parset, parsetsize + fixparsetsize, 
					//			nnode, idparset, d, &probMaxNode1, fMaxLogLik);

					/* release combination set */
        				for(ncomb = 0; ncomb < ncomblist; ncomb++) {
        	  				if(pcomblist[ncomb])
        	    					CATNET_FREE(pcomblist[ncomb]);
        	  				pcomblist[ncomb] = NULL;
					} /* for ncomb */
        				CATNET_FREE(pcomblist);
        				pcomblist = 0;
					ncomblist = 0;
	
				//} /* if(!getCachedProb) */

				for(k = 0; k < m_nCatnets; k++) {
					if(!m_pCatnets1[k] || !m_pCatnets2[k])
						continue;
					baseCatnet1 = *m_pCatnets1[k];
					baseCatnet1.setParents(nnode, idparset, d);
					baseCatnet1.setNodeProb(nnode, &probMaxNode1);
					complx = baseCatnet1.complexity();
					if(complx > maxComplexity) {
						continue;
					}
					loglik1 = baseCatnet1.loglik();
					baseCatnet2 = *m_pCatnets2[k];
					baseCatnet2.setParents(nnode, idparset, d);
					baseCatnet2.setNodeProb(nnode, &probMaxNode1);
					loglik2 = baseCatnet2.loglik();
					fdiff = loglik1/m_numSamples1 - loglik2/m_numSamples2;
					if(pCurCatnetScores[complx] < fdiff) {
						pCurCatnetScores[complx] = fdiff;
						//printf("pCurCatnetScores[%d] = %f\n", complx, pCurCatnetScores[complx]);
						if(pCurCatnetList1[complx])
							*pCurCatnetList1[complx] = baseCatnet1;
						else {
							pCurCatnetList1[complx] = new CATNET<t_node, t_node_size, t_prob>;
							*pCurCatnetList1[complx] = baseCatnet1;
						}
						if(pCurCatnetList2[complx])
							*pCurCatnetList2[complx] = baseCatnet2;
						else {
							pCurCatnetList2[complx] = new CATNET<t_node, t_node_size, t_prob>;
							*pCurCatnetList2[complx] = baseCatnet2;
						}
					}
				} /* for k */

			} /* for d */

			if(becho)
				printf("\n");

			for(j = 0; j < m_nCatnets; j++) {
				if(m_pCatnetScores[j] < pCurCatnetScores[j]) {
					//assert(pCurCatnetList1[j]);assert(pCurCatnetList2[j]);
					m_pCatnetScores[j] = pCurCatnetScores[j];
					if(m_pCatnets1[j]) 
						delete m_pCatnets1[j];
					m_pCatnets1[j] = pCurCatnetList1[j];
					pCurCatnetList1[j] = 0;
					if(m_pCatnets2[j]) 
						delete m_pCatnets2[j];
					m_pCatnets2[j] = pCurCatnetList2[j];
					pCurCatnetList2[j] = 0;
				}
				else {
					if(pCurCatnetList1[j]) {
						delete pCurCatnetList1[j];
						pCurCatnetList1[j] = 0;
					}
					if(pCurCatnetList2[j]) {
						delete pCurCatnetList2[j];
						pCurCatnetList2[j] = 0;
					}
				}
			}
		} // for(nnode = 0; nnode < numnodes; nnode++)

		CATNET_FREE(pCurCatnetList1);
		CATNET_FREE(pCurCatnetList2);
		CATNET_FREE(parset);
		CATNET_FREE(fixparset);
		CATNET_FREE(idparset);
		CATNET_FREE(pCurCatnetScores);

		if(psubsamples1)
			CATNET_FREE(psubsamples1);
		if(psubsamples2)
			CATNET_FREE(psubsamples2);

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
			if(m_pCatnets1[j]) {
				m_pCatnets1[j] -> normalizeProbabilities();
			}
			if(m_pCatnets2[j]) {
				delete m_pCatnets2[j];
				m_pCatnets2[j] = 0;
			}
		}

		return 0;
	}

};

#endif /* CATNET_CLASS2SEARCH_H */
