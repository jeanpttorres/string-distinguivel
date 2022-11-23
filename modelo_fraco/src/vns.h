#ifndef __VNS_H__
#define __VNS_H__

#include <time.h>
#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

int vns_main (SCIP* scip, char *filename, double tempo, double* value, double* initial, int *, int *, int*);
int neighbourhoodChange(int  *sum_sc, int *sum_sf, int *x, int *x2, int k, char **sc, char**sf, int sc_size, int sf_size, int kc, int kf, char *alphabet, int alpha_size, int string_size, int *fx, int *fx2);
void shaking (double **T, double **best_T, int tcount, int *x2, int k, double *ng, int string_size, int alpha_size, char* alphabet);
void lpImprovement(int  *sum_sc, int *sum_sf, double **T, double** best_T, double **Y, int *x, int *x2, int k, int ky, double *ng, double *ny, char **sc, char**sf, int sc_size, int sf_size, int kc, int kf, char *alphabet, int alpha_size, int string_size, SCIP_PROBDATA *probdata,  SCIP_SOL *bestSolution);
void lpShaking(int  *sum_sc, int *sum_sf, double **T, double** best_T, double **Y, int *x, int *x2, int k, int ky, double *ng, double *ny, char **sc, char**sf, int sc_size, int sf_size, int kc, int kf, char *alphabet, int alpha_size, int string_size, SCIP_PROBDATA *probdata,  SCIP_SOL *bestSolution);
void firstImprovement(int  *sum_sc, int *sum_sf, double **T, int *x2, int k, double *ng, char **sc, char**sf, int sc_size, int sf_size, int kc, int kf, char *alphabet, int alpha_size, int string_size);
void bestImprovement(int  *sum_sc, int *sum_sf, double **T, int *x2, int k, double *ng, char **sc, char**sf, int sc_size, int sf_size, int kc, int kf, char *alphabet, int alpha_size, int string_size, double time_max);
void imprimeSol(int* x, int* sum_sc, int* sum_sf, char** sc, char** sf, int sc_size, int sf_size, int kc, int kf, char* alphabet, int alpha_size, int string_size);

#ifdef __cplusplus
}
#endif

#endif
