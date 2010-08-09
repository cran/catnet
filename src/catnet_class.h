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
 * catnet_class.h
 *
 *  Created on: Sep 18, 2009
 *      Author: nbalov
 */

#include "utils.h"
#include "problist.h"

#ifndef CATNET_CLASS_H_
#define CATNET_CLASS_H_

//#define DEBUG_ON	1 

template<class t_node, int t_node_size, class t_prob>
class CATNET {
protected:
	/* nodes are assumed ordered */
	int m_numNodes;
	t_node **m_nodeNames;
	int m_maxParents;
	int *m_numParents;
	int **m_parents;
	int m_maxCategories;
	int *m_numCategories;
	int **m_catIndices;
	int m_complexity;
	t_prob m_loglik;
	PROB_LIST<t_prob> **m_pProbLists;

public:
	CATNET() {
		_reset();
	}

	CATNET(int nnodes, int maxpars, int maxcats = 2, const t_node **nodes = 0,
			const int * pnumpars = 0, const int **ppars = 0, const int *pcats = 0) {
		_reset();
		init(nnodes, maxpars, maxcats, nodes, pnumpars, ppars, pcats);
	}

	virtual ~CATNET() {
		_release();
	}

	CATNET<t_node, t_node_size, t_prob>& operator =(const CATNET<t_node,
			t_node_size, t_prob> &cnet) {
		init(cnet.m_numNodes, cnet.m_maxParents, cnet.m_maxCategories,
				(const t_node**)cnet.m_nodeNames, (const int*)cnet.m_numParents, (const int**)cnet.m_parents,
				(const int*)cnet.m_numCategories, (const PROB_LIST<t_prob> **)cnet.m_pProbLists);
		m_complexity = cnet.m_complexity;
		m_loglik = cnet.m_loglik;
		return *this;
	}

	virtual void _release() {
		int i;
		for (i = 0; i < m_numNodes; i++) {
			if (m_pProbLists && m_pProbLists[i]) {
				delete m_pProbLists[i];
				m_pProbLists[i] = 0;
			}
			if (m_parents && m_parents[i]) {
				CATNET_FREE(m_parents[i]);
				m_parents[i] = 0;
			}
			if (m_nodeNames && m_nodeNames[i]) {
				CATNET_FREE(m_nodeNames[i]);
				m_nodeNames[i] = 0;
			}
			if(m_catIndices && m_catIndices[i]) {
				CATNET_FREE(m_catIndices[i]);
				m_catIndices[i] = 0;
			}
		}
		if (m_numParents)
			CATNET_FREE(m_numParents);
		if (m_parents)
			CATNET_FREE(m_parents);
		if (m_numCategories)
			CATNET_FREE(m_numCategories);
		if (m_nodeNames)
			CATNET_FREE(m_nodeNames);
		if(m_catIndices)
			CATNET_FREE(m_catIndices);
		if (m_pProbLists)
			CATNET_FREE(m_pProbLists);

		_reset();
	}

	virtual void _reset() {
		m_numNodes = 0;
		m_maxParents = 0;
		m_maxCategories = 0;
		m_nodeNames = 0;
		m_numParents = 0;
		m_parents = 0;
		m_numCategories = 0;
		m_catIndices = 0;
		m_pProbLists = 0;
		m_complexity = 0;
		m_loglik = 0;
	}

	void init(int nnodes, int maxpars, int maxcats = 2, const t_node **nodes = 0,
			const int *pnumpars = 0, const int **ppars = 0, const int *pcats = 0, 
			const PROB_LIST<t_prob> **pprobs = 0) {

		_release();

		int i, j, *nodeparcats, nodenamelen;

		if (nnodes < 1 || maxpars < 0)
			return;
		m_numNodes = nnodes;
		m_maxParents = maxpars;
		m_maxCategories = maxcats;

		m_numParents = (int*) CATNET_MALLOC(m_numNodes * sizeof(int));
		m_parents = (int**) CATNET_MALLOC(m_numNodes * sizeof(int*));
		m_numCategories = (int*) CATNET_MALLOC(m_numNodes * sizeof(int));
		m_catIndices = (int**) CATNET_MALLOC(m_numNodes * sizeof(int*));
		for (i = 0; i < m_numNodes; i++) {
			m_numParents[i] = 0;
			m_parents[i] = 0;
			m_numCategories[i] = m_maxCategories;
			m_catIndices[i] = 0;
		}
		m_nodeNames = (t_node**) CATNET_MALLOC(m_numNodes * sizeof(t_node*));
		if (nodes) {	
			for (i = 0; i < m_numNodes; i++) {
				m_nodeNames[i] = 0;
				if (!nodes[i])
					continue;
				nodenamelen = strlen(nodes[i]);
				m_nodeNames[i] = (t_node*) CATNET_MALLOC((nodenamelen+1) * sizeof(t_node));
				strcpy(m_nodeNames[i], nodes[i]);
			}
		}
		else {
			memset(m_nodeNames, 0, m_numNodes * sizeof(t_node*));
		}

		if (pnumpars != 0 && ppars != 0) {
			memcpy(m_numParents, pnumpars, m_numNodes * sizeof(int));
			for (i = 0; i < m_numNodes; i++) {
				m_parents[i] = (int*) CATNET_MALLOC(m_numParents[i] * sizeof(int));
				memset(m_parents[i], 0, m_numParents[i] * sizeof(int));
				if(ppars[i])
					memcpy(m_parents[i], ppars[i], m_numParents[i] * sizeof(int));
			}
		}
		if (pcats != 0)
			memcpy(m_numCategories, pcats, m_numNodes * sizeof(int));
		
		nodeparcats = (int*)CATNET_MALLOC((m_maxParents+1)*sizeof(int));

		m_pProbLists = (PROB_LIST<t_prob>**) CATNET_MALLOC(m_numNodes * sizeof(PROB_LIST<t_prob>*));
		memset(m_pProbLists, 0, m_numNodes * sizeof(PROB_LIST<t_prob>*));
		for (i = 0; i < m_numNodes; i++) {
			if (pprobs && pprobs[i]) {
				m_pProbLists[i] = new PROB_LIST<t_prob>;
				*m_pProbLists[i] = *pprobs[i];
			}
			else {
				for(j = 0; j < m_numParents[i]; j++)
					nodeparcats[j] = m_numCategories[m_parents[i][j]]; 
				m_pProbLists[i] = new PROB_LIST<t_prob>(m_numCategories[i], m_maxCategories, m_numParents[i], nodeparcats);
			}
		}
		
		CATNET_FREE(nodeparcats);

	}

	void setNodeNames(char **pnames, const int *porder) {
		int i;
		const char *str;
		if (!porder || !pnames) 
			return;
		if(!m_nodeNames) {
			m_nodeNames = (t_node**) CATNET_MALLOC(m_numNodes * sizeof(t_node*));
			memset(m_nodeNames, 0, m_numNodes * sizeof(t_node*));
		}
		for (i = 0; i < m_numNodes; i++) {
			m_nodeNames[i] = 0;
			if (porder[i] < 1 || porder[i] > m_numNodes)
				break;
			str = pnames[porder[i]-1];
			if(m_nodeNames[i])
				CATNET_FREE(m_nodeNames[i]);
			m_nodeNames[i] = (t_node*) CATNET_MALLOC((strlen(str)+1) * sizeof(char));
			strcpy((char*)m_nodeNames[i], str);
		}
	}

	void setNodesOrder(const int *porder) {
		int i;
		char str[256];
		if (!porder) 
			return;
		if(!m_nodeNames) {
			m_nodeNames = (t_node**) CATNET_MALLOC(m_numNodes * sizeof(t_node*));
			memset(m_nodeNames, 0, m_numNodes * sizeof(t_node*));
		}
		for (i = 0; i < m_numNodes; i++) {
			m_nodeNames[i] = 0;
			if (porder[i] < 1 || porder[i] > m_numNodes)
				break;
			sprintf(str, "N%d", (int)porder[i]);
			if(m_nodeNames[i])
				CATNET_FREE(m_nodeNames[i]);
			m_nodeNames[i] = (t_node*) CATNET_MALLOC((strlen(str)+1) * sizeof(char));
			strcpy((char*)m_nodeNames[i], str);
		}
	}

	void setCategoryIndices(int *pNumCats, int **pCatIndices) {
		int i;
		// reset node names only
		if (!pNumCats || !pCatIndices || !m_numCategories) 
			return;
		//if(m_catIndices)
		//	CATNET_FREE(m_catIndices);
		if(!m_catIndices)
			m_catIndices = (int**) CATNET_MALLOC(m_numNodes * sizeof(int*));
		for (i = 0; i < m_numNodes; i++) {
			if(pNumCats[i] != m_numCategories[i])
				break;
			if(!m_catIndices[i])
				m_catIndices[i] = (int*) CATNET_MALLOC(m_numCategories[i] * sizeof(int));
			memcpy(m_catIndices[i], pCatIndices[i], m_numCategories[i] * sizeof(int));
		}
	}

	void setCondProb(int node, t_prob *pcondprob, int condprobsize) {
		if (m_numNodes < 1 || node < 0 || node >= m_numNodes)
			return;
		if (!m_pProbLists)
			m_pProbLists = (PROB_LIST<t_prob>**) CATNET_MALLOC(m_numNodes
					* sizeof(PROB_LIST<t_prob>*));
		if (m_pProbLists && m_pProbLists[node])
			delete m_pProbLists[node];
		m_pProbLists[node] = 0;
		if (m_numParents[node] < 0 || m_numParents[node] > m_maxParents)
			return;
		int *parcats = (int*) CATNET_MALLOC(m_maxParents * sizeof(int));
		for (int i = 0; i < m_numParents[node]; i++) {
			parcats[i] = m_numCategories[m_parents[node][i]];
			//cout << "parcats[i] = " << parcats[i] << "\n";
		}
		m_pProbLists[node] = new PROB_LIST<t_prob> (m_numCategories[node],
				m_maxCategories, m_numParents[node], parcats, pcondprob,
				condprobsize);
		CATNET_FREE(parcats);
	}

	int numNodes() {
		return m_numNodes;
	}

	const t_node** nodeNames() {
		return (const t_node**)m_nodeNames;
	}

	int maxCategories() {
		return m_maxCategories;
	}

	const int* numCategories() {
		return (const int*)m_numCategories;
	}

	int numCategories(int node) {
		if (node < 0 || node >= m_numNodes)
			return 0;
		return m_numCategories[node];
	}

	int maxParents() {
		return m_maxParents;
	}

	const int** parents() {
		return (const int**)m_parents;
	}

	const int* numParents() {
		return (const int*)m_numParents;
	}

	int numParents(int node) {
		if (node < 0 || node >= m_numNodes)
			return 0;
		return m_numParents[node];
	}

	const int* getParents(int node) {
		if (node < 0 || node >= m_numNodes)
			return 0;
		return (const int*)m_parents[node];
	}

	int setParents(int node, int* parents, int numparents) {
		if (node < 0 || node >= m_numNodes)
			return 0;
		if (!m_numParents || !m_parents)
			return 0;
		if(m_numParents[node] != numparents) {
			m_numParents[node] = numparents;
			if(m_parents[node])
				CATNET_FREE(m_parents[node]);
			m_parents[node] = (int*) CATNET_MALLOC(m_numParents[node] * sizeof(int));
		}
		memcpy(m_parents[node], parents, m_numParents[node] * sizeof(int));

		if (!m_pProbLists)
			m_pProbLists = (PROB_LIST<t_prob>**) CATNET_MALLOC(m_numNodes * sizeof(PROB_LIST<t_prob>*));
		else {
			if(m_pProbLists[node])
				delete m_pProbLists[node];
		}
		int *parcats = (int*)CATNET_MALLOC(m_numParents[node]*sizeof(int));
		for(int i = 0; i < m_numParents[node]; i++) 
			parcats[i] = m_numCategories[parents[i]];	
		m_pProbLists[node] = new PROB_LIST<t_prob>(m_numCategories[node], m_maxCategories, m_numParents[node], parcats);
		CATNET_FREE(parcats);
		             
		if(m_maxParents < numparents)
			m_maxParents = numparents;

		/* need to be calculated next time */
		m_complexity = 0;
		m_loglik = 0;

		return m_numParents[node];
	}

	int complexity() {
		if(m_complexity < m_numNodes)
			return findComplexity();
		return m_complexity;
	}

	int findComplexity() {
		int i, j, c;
		m_complexity = 0;
		for (i = 0; i < m_numNodes; i++) {
			if (!m_parents || !m_parents[i]) {
				m_complexity += (m_numCategories[i]-1);
				continue;
			}
			c = 1;
			for (j = 0; j < m_numParents[i]; j++)
				c *= m_numCategories[m_parents[i][j]];
			m_complexity += c*(m_numCategories[i]-1);
		}
		return m_complexity;
	}

	int nodeComplexity(int nnode) {
		int j, c;
		if(nnode < 0 || nnode >= m_numNodes)
			return(0);
		if (!m_parents || !m_parents[nnode])
			return(m_numCategories[nnode]-1);
		c = (m_numCategories[nnode]-1);
		for (j = 0; j < m_numParents[nnode]; j++)
			c *= m_numCategories[m_parents[nnode][j]];
		return c;
	}

	t_prob loglik() {
		if(m_loglik == 0)
			return findLoglik();
		return m_loglik;
	}

	t_prob findLoglik() {
		int i;
		if(!m_pProbLists)
			return -FLT_MAX;
		m_loglik = 0;
		for (i = 0; i < m_numNodes; i++) {
			if(m_pProbLists[i])
				m_loglik += (m_pProbLists[i]->loglik + m_pProbLists[i]->priorlik);
		}
		return m_loglik;
	}
	
	const PROB_LIST<t_prob>** probLists() {
		return (const PROB_LIST<t_prob>**)m_pProbLists;
	}

	PROB_LIST<t_prob>* getNodeProb(int nnode) {
		if (!m_pProbLists || nnode < 0 || nnode >= m_numNodes)
			return 0;
		return m_pProbLists[nnode];
	}

	const PROB_LIST<t_prob>* setNodeProb(int nnode, const PROB_LIST<t_prob> *pprob) {
		if(!m_pProbLists || nnode < 0 || nnode >= m_numNodes)
			return 0;
		if(!m_pProbLists[nnode])
			m_pProbLists[nnode] = new PROB_LIST<t_prob>;
		if (m_pProbLists[nnode] && pprob) {
			*m_pProbLists[nnode] = *pprob;
			// normalize
			//m_pProbLists[nnode]->normalize();
		}
		return m_pProbLists[nnode];
	}

	void setNodePriorProb(int nnode, t_prob pprior) {
		if(!m_pProbLists || nnode < 0 || nnode >= m_numNodes || !m_pProbLists[nnode])
			return;
		m_pProbLists[nnode]->priorlik = pprior;
	}

	t_prob nodeSampleLoglik(int nnode, int *pnodepars, int nodepars,
			int *psamples, int nsamples) {
		int j, ipar;
		/* psamples have categories in the range [1, m_maxCategories] */
		t_prob *pnodeprob, floglik = 0;
		int *pnodesample, samp;

		if(!psamples || nsamples < 1)
			return 0;

		pnodesample = (int*) CATNET_MALLOC(m_maxParents * sizeof(int));
		pnodepars = m_parents[nnode];
		for (j = 0; j < nsamples; j++) {
			for (ipar = 0; ipar < nodepars; ipar++) {
				pnodesample[ipar] = psamples[j * m_numNodes + pnodepars[ipar]];
			}
			pnodeprob = m_pProbLists[nnode]->find_slot(0, pnodesample, 0);
			samp = psamples[j * m_numNodes + nnode];
			if (samp >= 0 && samp < m_numCategories[nnode] && pnodeprob[samp] > 0)
				floglik += (t_prob)log((double)pnodeprob[samp]);
		}
		CATNET_FREE(pnodesample);
		return(floglik);
	}

	t_prob sampleLoglik(int *psamples, int nsamples) {
		int i, j, ipar;
		/* psamples have categories in the range [1, m_maxCategories] */
		t_prob *pnodeprob;
		int nodepars;
		int *pnodepars, *pnodesample, samp;

		if(!psamples || nsamples < 1)
			return 0;

		pnodesample = (int*) CATNET_MALLOC(m_maxParents * sizeof(int));
		m_loglik = 0;
		for (i = 0; i < m_numNodes; i++) {
			if(!m_pProbLists[i])
				continue;
			pnodepars = m_parents[i];
			nodepars = m_numParents[i];
			for (j = 0; j < nsamples; j++) {
				for (ipar = 0; ipar < nodepars; ipar++) {
					pnodesample[ipar] = psamples[j * m_numNodes + pnodepars[ipar]];
				}
				pnodeprob = m_pProbLists[i]->find_slot(0, pnodesample, 0);
				samp = psamples[j * m_numNodes + i];
				if (pnodeprob && samp >= 0 && samp < m_numCategories[i] && pnodeprob[samp] > 0)
					m_loglik += (t_prob)log((double)pnodeprob[samp]);
			}
			//m_loglik += nsamples*m_pProbLists[i]->priorlik;
		}
		CATNET_FREE(pnodesample);
		return m_loglik;
	}

	t_prob* sampleLoglikVector(int *psamples, int nsamples) {
		int i, j, ipar;
		/* psamples have categories in the range [1, m_maxCategories] */
		t_prob *pnodeprob, *ploglik;
		int nodepars;
		int *pnodepars, *pnodesample, samp;

		if(!psamples || nsamples < 1)
			return 0;

		ploglik = (t_prob*) CATNET_MALLOC(nsamples * sizeof(t_prob));
		for (j = 0; j < nsamples; j++) 
			ploglik[j] = 0;

		pnodesample = (int*) CATNET_MALLOC(m_maxParents * sizeof(int));
		for (i = 0; i < m_numNodes; i++) {
			if(!m_pProbLists[i])
				continue;
			pnodepars = m_parents[i];
			nodepars = m_numParents[i];
			//ploglik[j] = m_pProbLists[i]->priorlik;
			for (j = 0; j < nsamples; j++) {
				for (ipar = 0; ipar < nodepars; ipar++) {
					pnodesample[ipar] = psamples[j * m_numNodes + pnodepars[ipar]];
				}
				pnodeprob = m_pProbLists[i]->find_slot(0, pnodesample, 0);
				samp = psamples[j * m_numNodes + i];
				if (pnodeprob && samp >= 0 && samp < m_numCategories[i] && pnodeprob[samp] > 0)
					ploglik[j] += (t_prob)log((double)pnodeprob[samp]);
			}
		}
		CATNET_FREE(pnodesample);
		return ploglik;
	}

	t_prob sampleNodeLoglik(int nnode, int *psamples, int nsamples) {
		int j, ipar;
		/* psamples have categories in the range [1, m_maxCategories] */
		t_prob *pnodeprob, loglik;
		int nodepars;
		int *pnodepars, *pnodesample=0, samp;

		if(!psamples || nsamples < 1 || nnode < 0 || nnode >= m_numNodes)
			return 0;

		if(!m_pProbLists || !m_pProbLists[nnode])
			return 0;

		loglik = 0;

		pnodesample = (int*) CATNET_MALLOC(m_maxParents * sizeof(int));	
		pnodepars = m_parents[nnode];
		nodepars = m_numParents[nnode];
		
		for (j = 0; j < nsamples; j++) {
			for (ipar = 0; ipar < nodepars; ipar++) {
				pnodesample[ipar] = psamples[j * m_numNodes + pnodepars[ipar]];
			}
			pnodeprob = m_pProbLists[nnode]->find_slot(0, pnodesample, 0);
			samp = psamples[j * m_numNodes + nnode];
			if (pnodeprob && samp >= 0 && samp < m_numCategories[nnode] && pnodeprob[samp] > 0)
				loglik += (t_prob)log((double)pnodeprob[samp]);
		}
		CATNET_FREE(pnodesample);

		return(loglik);
	}

	// sets sample conditional probability and returns its log-likelihood
	t_prob setNodeSampleProb(int nnode, int *psamples, int nsamples, int bNormalize = 0) {
		int i, j;
		/* psamples have categories in the range [1, m_maxCategories] */
		t_prob *pnodeprob, floglik;
		int *pnodesample, *pnodepars, samp;

		if (!m_pProbLists || !psamples || nsamples < 1)
			return 0;

		if(!m_pProbLists[nnode]) {
			// error
			return 0;
		}
		else
			m_pProbLists[nnode]->set_zero();

		pnodesample = (int*) CATNET_MALLOC(m_maxParents * sizeof(int));
		pnodepars = m_parents[nnode];
		for (j = 0; j < nsamples; j++) {
			for (i = 0; i < m_numParents[nnode]; i++) {
				if (pnodepars[i] < 0 || pnodepars[i] >= m_numNodes)
					break;
				pnodesample[i] = psamples[j * m_numNodes + pnodepars[i]];
			}
			pnodeprob = m_pProbLists[nnode]->find_slot(0, pnodesample, 0);
			samp = psamples[j * m_numNodes + nnode];
			if (pnodeprob && samp >= 0 && samp < m_numCategories[nnode])
				pnodeprob[samp]++;
		}

		CATNET_FREE(pnodesample);

		// at this point m_pProbLists[nnode] has counts not probabilities 
		floglik = m_pProbLists[nnode]->loglikelihood();

		if(bNormalize)
			m_pProbLists[nnode] -> normalize();

		/* need to be calculated next time */
		m_loglik = 0;

		return(floglik);
	}

	t_prob *findNodeSampleProbError(int nnode, int *psamples, int nsamples, int &nerrors) {
		int i, j, k, cx, cy, *ny, nerr;
		/* psamples have categories in the range [1, m_maxCategories] */
		t_prob *pnodeprob, floglik, *My, h, hath, *perror, dd, delta, faux, perr1, perr2;
		int *pnodesample, *pnodepars, samp;

		if(!psamples || nsamples < 1 || nnode < 0 || nnode >= m_numNodes)
			return 0;

		if(!m_pProbLists || !m_pProbLists[nnode])
			return 0;

		PROB_LIST<t_prob> *pProbList = new PROB_LIST<t_prob>;
		*pProbList = *m_pProbLists[nnode];
		pProbList->set_zero();

		pnodesample = (int*) CATNET_MALLOC(m_maxParents * sizeof(int));
		pnodepars = m_parents[nnode];

		for (j = 0; j < nsamples; j++) {
			for (i = 0; i < m_numParents[nnode]; i++) {
				if (pnodepars[i] < 0 || pnodepars[i] >= m_numNodes)
					break;
				pnodesample[i] = psamples[j * m_numNodes + pnodepars[i]];
			}
			pnodeprob = pProbList->find_slot(0, pnodesample, 0);
			samp = psamples[j * m_numNodes + nnode];
			if (pnodeprob && samp >= 0 && samp < m_numCategories[nnode])
				pnodeprob[samp]++;
		}

		CATNET_FREE(pnodesample);

		// the number of parent configurations
		cx = pProbList->numCats;
		cy = pProbList->nProbSize / pProbList->numCats;
		// number of samples per parent configuration
		ny = (int*)CATNET_MALLOC(cy*sizeof(int));
		// maximum entropy per configuration
		My = (t_prob*)CATNET_MALLOC(cy*sizeof(t_prob));

printf("probSize = %d, numCats = %d, cx = %d,  cy = %d, nsamples = %d\n", m_pProbLists[nnode]->nProbSize, m_pProbLists[nnode]->numCats, cx, cy, nsamples);

		// entropies
		h = hath = 0;

		k = 0;
		for(j = 0; j < cy; j++) {
			ny[j] = 0;
			My[j] = 0;
			perr1 = 0;
			perr2 = 0;
			for (i = 0; i < cx; i++) {
				floglik = m_pProbLists[nnode]->pProbs[k + i];
				if(floglik > 0) 
					floglik *= (t_prob)log((double)floglik);
				perr1 += floglik;

				ny[j] += pProbList->pProbs[k + i];
			}
			if(ny[j] <= 0)
				continue;
			for (i = 0; i < cx; i++) {
				pProbList->pProbs[k + i] /= ny[j];
				floglik = pProbList->pProbs[k + i];
				if(floglik > 0) 
					floglik *= (t_prob)log((double)floglik);
				perr2 += floglik;
				My[j] += floglik;
			}
			if(My[j] < 0) My[j] = -My[j];

			h += (ny[j]*perr1/nsamples);
			hath += (ny[j]*perr2/nsamples);

			k += cx;
printf("ny[%d] = %d,  My[%d] = %f\n", j, ny[j], j, My[j]);
		}

		delta = (h - hath);
		if(delta < 0) delta = -delta;
		h = -h;
		hath = -hath;

		printf("h = %f, hath = %f, delta = %f\n", h, hath, delta);
		
		nerrors = 10;//(int)(2*abs(h) / delta) + 1;
		perror = (t_prob*)CATNET_MALLOC(nerrors*sizeof(t_prob));

		for(nerr = 0; nerr < nerrors; nerr++) { 
			perror[nerr] = 0;
			dd = delta + (t_prob)(4*nerr*nerr)*h/(t_prob)(nerrors*nerrors);
			faux = 4*cy*cy/(nsamples*dd*dd);

			perr1 = perr2 = 0;
			k = 0;
			for(j = 0; j < cy; j++) {
				for (i = 0; i < cx; i++) {
					perr1 += _gamma_upper_bound(m_pProbLists[nnode]->pProbs[k + i], ny[j], (dd)/(2*cx));
				}
				perr2 += faux*My[j]*My[j]*ny[j]*(nsamples-ny[j])/(nsamples*nsamples);
				k += cx;
			}
printf("dd = %f,  faux = %f, perr1 = %f, perr2 = %f\n", dd, faux, perr1, perr2);

			perror[nerr] = perr1 + perr2;
		}

		CATNET_FREE(ny);
		CATNET_FREE(My);
		delete pProbList;

		return(perror);
	}

	int normalizeProbabilities() {
		int i;
		if(!m_pProbLists)
			return -1;
		for (i = 0; i < m_numNodes; i++) {
			if(m_pProbLists[i])
				m_pProbLists[i] -> normalize();
		}
		return 0;
	}

	int *findParentPool(int nnode, int &poolsize) {

		int ipar, par, i, j, bfound, parpoolsize, *parpool, *ppool, *paux;

		poolsize = 0;
		if (nnode < 0 || nnode >= m_numNodes || !m_parents || !m_parents[nnode] || !m_numParents[nnode])
			return 0;

		//cout << "nnode = " << nnode+1 << ", m_numParents[nnode] = " << m_numParents[nnode] << "\n";

		ppool = 0;
		for (ipar = 0; ipar < m_numParents[nnode]; ipar++) {
			par = m_parents[nnode][ipar];
			parpool = findParentPool(par, parpoolsize);

			i = 0;
			while (i < parpoolsize) {
				bfound = 0;
				for (j = 0; j < poolsize; j++) {
					if (parpool[i] == ppool[j]) {
						bfound = 1;
						break;
					}
				}
				if (!bfound) {
					i++;
					continue;
				}
				for (j = i + 1; j < parpoolsize; j++)
					parpool[j - 1] = parpool[j];
				parpoolsize--;
			}

			//cout << "par = " << par+1 << "poolsize = " << poolsize << ", parpoolsize = " << parpoolsize << "\n";
			paux = (int*) CATNET_MALLOC((poolsize + parpoolsize + 1) * sizeof(int));
			if (poolsize > 0 && ppool)
				memcpy(paux, ppool, poolsize * sizeof(int));
			if (parpoolsize > 0 && parpool)
				memcpy(paux + poolsize, parpool, parpoolsize * sizeof(int));
			// release that
			if(parpool)
				CATNET_FREE(parpool);

			// check par
			bfound = 0;
			for (j = 0; j < poolsize + parpoolsize; j++) {
				if (paux[i] == par) {
					bfound = 1;
					break;
				}
			}
			if (bfound) {
				poolsize += parpoolsize;
			} else {
				paux[poolsize + parpoolsize] = par;
				poolsize += (parpoolsize + 1);
			}

			if (ppool)
				CATNET_FREE(ppool);
			ppool = paux;
		}
		return ppool;
	}

	t_prob *findJointProb(int nnode, int &jointprobsize) {
		t_prob *jointprob, *parprob;
		int i, ii, i0, ic, ipar, j, par, ipool, poolnode;
		int *parpool, *paux, parpoolsize, *blocksize, *paridx, *pcats;
		PROB_LIST<t_prob> *probnode = 0;

		if (nnode < 0 || nnode >= m_numNodes || !m_pProbLists)
			return 0;

		parpool = findParentPool(nnode, parpoolsize);
		// add nnode to the parents list
		paux = (int*) CATNET_MALLOC((parpoolsize + 1) * sizeof(int));
		if (parpool && parpoolsize > 0) 
			memcpy(paux, parpool, parpoolsize * sizeof(int));
		if (parpool)
			CATNET_FREE(parpool);
		paux[parpoolsize] = nnode;
		parpoolsize++;
		parpool = paux;

		pcats = (int*) CATNET_MALLOC(m_maxParents * sizeof(int));
		paridx = (int*) CATNET_MALLOC(parpoolsize * sizeof(int));

		blocksize = (int*) CATNET_MALLOC(parpoolsize * sizeof(int));
		blocksize[parpoolsize - 1] = 1;
		for (i = parpoolsize - 2; i >= 0; i--) {
			blocksize[i] = blocksize[i + 1] * m_numCategories[parpool[i + 1]];
		}

		jointprobsize = 1;
		for (i = 0; i < parpoolsize; i++) {
			jointprobsize *= m_numCategories[parpool[i]];
		}

		jointprob = (t_prob*) CATNET_MALLOC(jointprobsize * sizeof(t_prob));
		if (!jointprob)
			return jointprob;
		for (i = 0; i < jointprobsize; i++)
			jointprob[i] = 1.0;

		for (ipool = 0; ipool < parpoolsize; ipool++) {
			poolnode = parpool[ipool];
			probnode = m_pProbLists[poolnode];

			if (m_numParents[poolnode] == 0) {
				for (ii = 0; ii < jointprobsize; ii += (blocksize[ipool]
						* m_numCategories[poolnode])) {
					for (ic = 0; ic < m_numCategories[poolnode]; ic++) {
						i0 = ic * blocksize[ipool];
						for (i = 0; i < blocksize[ipool]; i++) {
							jointprob[ii + i0 + i] *= probnode->pProbs[ic];
						}
					}
				}
				continue;
			}

			//cout << "Parents: ";
			memset(paridx, 0, parpoolsize * sizeof(int));
			for (ipar = 0; ipar < m_numParents[poolnode]; ipar++) {
				par = m_parents[poolnode][ipar];
				for (i = 0; i < ipool; i++)
					if (par == parpool[i])
						break;
				if (i >= ipool) {
					// bad
					break;
				}
				paridx[i] = ipar + 1;
			}

			memset(pcats, 0, m_maxParents * sizeof(int));
			for (ii = 0; ii < jointprobsize; ii += (blocksize[ipool] * m_numCategories[poolnode])) {
				for (j = 0; j < ipool; j++) {
					if (paridx[j] > 0) {
						ic = (int) (ii / blocksize[j]);
						ic -= m_numCategories[poolnode] * (int) (ic / m_numCategories[poolnode]);
						pcats[paridx[j] - 1] = ic;
					}
				}
				parprob = probnode->find_slot(0, pcats, 0);
				if (!parpool)
					continue;
				for (ic = 0; ic < m_numCategories[poolnode]; ic++) {
					i0 = ic * blocksize[ipool];
					for (i = 0; i < blocksize[ipool]; i++) {
						jointprob[ii + i0 + i] *= parprob[ic];
					}
				}
			}
		}

		CATNET_FREE(pcats);
		CATNET_FREE(paridx);
		CATNET_FREE(parpool);
		CATNET_FREE(blocksize);

		return jointprob;
	}

	t_prob *marginal_prob(int nnode) {

		t_prob *jointprob, *margprob;
		int jointprobsize = 0, nodecats, ic, k;
		nodecats = m_numCategories[nnode];
		margprob = (t_prob*) CATNET_MALLOC(nodecats * sizeof(t_prob));
		jointprob = findJointProb(nnode, jointprobsize);
		if (!jointprob)
			return 0;
		for (ic = 0; ic < nodecats; ic++) {
			k = ic;
			margprob[ic] = 0;
			while (k < jointprobsize) {
				margprob[ic] += jointprob[k];
				k += nodecats;
			}
		}
		CATNET_FREE(jointprob);

		return margprob;
	}

};

#endif /* CATNET_CLASS_H_ */
