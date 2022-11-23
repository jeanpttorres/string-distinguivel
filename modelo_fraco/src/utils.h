#ifndef __UTILS_H__
#define __UTILS_H__

#include <stdbool.h>
#include "scip/scip.h"
#include "probdata_dssp.h"

#ifdef __cplusplus
extern "C" {
#endif

int randomInteger (int low, int high);
void validSol(SCIP *, char *);
int valorSolucao(SCIP*, char *, int *, int *, int*);
void printRelaxedSol(SCIP *, char *, double);
int objective_value(int  *sum_sc, int *sum_sf, char **Sc, char **Sf, int sc_size, int sf_size, int kc, int kf, char *alpha, int alpha_size, int string_size, int *solution);
int dc_df_objective_value(int  *sum_sc, int *sum_sf, char **Sc, char **Sf, int sc_size, int sf_size, int kc, int kf, char *alpha, int alpha_size, int string_size, int *solution, int* dc, int* df);

#ifdef __cplusplus
}
#endif

#endif

