/*
 * rcatnet.h
 *
 *  Created on: Sep 21, 2009
 *      Author: nbalov
 */

#ifndef RCATNET_H_
#define RCATNET_H_

#include "catnet_class.h"

#define MAX_NODE_NAME	16
class RCatnet : public CATNET<char, MAX_NODE_NAME, double> {
public:
	RCatnet();
	RCatnet(SEXP cnet);

	SEXP genRcatnet(const char * objectName);
	SEXP genProbList(int node, int paridx, int *pcats);

	RCatnet& operator =(CATNET<char, MAX_NODE_NAME, double> &cnet) {
		CATNET<char, MAX_NODE_NAME, double>::init(
 			cnet.numNodes(), cnet.maxParents(), cnet.maxCategories(),
			cnet.nodeNames(), cnet.numParents(), cnet.parents(),
			cnet.numCategories(), cnet.probLists());
		return *this;
	}
};

extern "C" {

char *gen_prob_string(int node, SEXP parlist, int paridx, SEXP catlist, SEXP problist, char *str);
SEXP prob_string(SEXP rnodes, SEXP rparents, SEXP rcatlist, SEXP rproblist);
void gen_prob_vector(int node, SEXP parlist, int paridx, SEXP catlist, SEXP problist, double *&pvec, int &nvec);
SEXP prob_vector(SEXP rnodes, SEXP rparents, SEXP rcatlist, SEXP rproblist); 

} // extern "C"


#endif /* RCATNET_H_ */
