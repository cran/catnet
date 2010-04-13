
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <R_ext/RConverters.h>
#include <R_ext/Rdynload.h>

/* these are called by R-functions directly */

SEXP createRCatnet(SEXP cnet);
SEXP catnetMarginalProb(SEXP cnet, SEXP rnode);
SEXP catnetJointProb(SEXP cnet, SEXP rnode);
SEXP catnetFindParentPool(SEXP cnet, SEXP rnode);
SEXP showCatnet(SEXP cnet);
SEXP catnetOptimalNetsForOrder(SEXP rSamples, SEXP rPerturbations,
				SEXP rMaxParents, SEXP rMaxComplexity, SEXP rOrder,
				SEXP rParentsPool, SEXP rFixedParentsPool, SEXP rUseCache, SEXP rEcho);
SEXP catnetOptimalClass2NetsForOrder(SEXP rSamples1, SEXP rPerturbations1,
				SEXP rSamples2, SEXP rPerturbations2,
				SEXP rMaxParents, SEXP rMaxComplexity, SEXP rOrder,
				SEXP rParentsPool, SEXP rFixedParentsPool, SEXP rUseCache, SEXP rEcho);
SEXP catnetLoglik(SEXP cnet, SEXP rSamples, SEXP rPerturbations);
SEXP catnetSetProb(SEXP cnet, SEXP rSamples, SEXP rPerturbations);
SEXP catnetReleaseCache();
