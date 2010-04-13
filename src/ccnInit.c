
#include "catnet_rexport.h"

size_t g_memcounter = 0;

static const R_CallMethodDef R_CallDef[] = {
	{"ccnOptimalNetsForOrder", (DL_FUNC)&catnetOptimalNetsForOrder, 9},
	{"ccnOptimalClass2NetsForOrder", (DL_FUNC)&catnetOptimalClass2NetsForOrder, 11},
	{"ccnSetProb", (DL_FUNC)&catnetSetProb, 3},
	{"ccnLoglik", (DL_FUNC)&catnetLoglik, 3},
	{"ccnMarginalProb", (DL_FUNC)&catnetMarginalProb, 2},
	{"ccnReleaseCache", (DL_FUNC)&catnetReleaseCache, 0},
	{NULL, NULL, 0},
};

void R_init_catnet(DllInfo *info)
{
	R_registerRoutines(info,NULL,R_CallDef,NULL,NULL);
	//printf("ccnInit\n");
	g_memcounter = 0;
}

void R_unload_catnet(DllInfo *info)
{
	//printf("ccnUnload\n");
	catnetReleaseCache();
}
