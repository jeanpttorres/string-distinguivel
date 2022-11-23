#ifndef __BCPA_H__
#define __BCPA_H__

#include <time.h>
#include "scip/scip.h"
#include "utils.h"

#ifdef __cplusplus
extern "C" {
#endif

int heur_BCPA(SCIP* scip_o, char *filename, double *value, int *x, int *dc, int *df);
SCIP_RETCODE setupProblem(
    SCIP*                 scip                /**< SCIP data structure */,
    SCIP*                 scip_o                /**< SCIP data structure */,
    SCIP_VAR*** vars
    );

#ifdef __cplusplus
}
#endif

#endif
