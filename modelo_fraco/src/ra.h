#ifndef __RA_H__
#define __RA_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

int heur_RA(SCIP* scip, char *filename, double tempo, double *value, int* x, int *dc, int *df);

#ifdef __cplusplus
}
#endif

#endif
