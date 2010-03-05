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

#endif /* RCATNET_H_ */
