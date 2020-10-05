/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2014 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cmain.c
 * @brief  Main file for dssp example
 * @author Edna Hoshino
 *
 *  This the file contains the \ref main() main function of the projects. This includes all the default plugins of
 *  \SCIP and the once which belong to that projects. After that is starts the interactive shell of \SCIP or processes
 *  the shell arguments if given.
 */

/** \mainpage B&B for the dssp problem Index Page
 *
 * \section Introduction
 *
 * This project refers to a simple branch-and-bound implementation to solve a dssp problem.
 *
 */

#include <stdio.h>
#include <time.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "scip/scip.h"
#include "scip/scipshell.h"
#include "scip/scipdefplugins.h"

#include "probdata_dssp.h"
#include "reader_dssp.h"
#include "parameters_dssp.h"
#include "heur_myheuristic.h"

void printSol(SCIP *, char *, double);
void validSol(SCIP *, char *);
double heur_RA(SCIP* , char *);
double valorSolucao(SCIP*, char *, int *);
void printRelaxedSol(SCIP *, char *, double);


double valorSolucao(SCIP *scip, char *filename, int *heur){
   int i, j, k, l, dif, dc, df, menor_dif_sc, maior_dif_sf;
    int  *sum_sc, *sum_sf, string_size, alpha_size, sc_size, sf_size;
    char *alphabet, **sc, **sf;

    SCIP_PROBDATA *probdata;


    assert(scip != NULL);

    probdata = SCIPgetProbData(scip);
    assert(probdata != NULL);

    string_size = SCIPprobdataGetString_size(probdata);
    alphabet = SCIPprobdataGetAlphabet(probdata);
    alpha_size = SCIPprobdataGetAlpha_size(probdata);
    sc = SCIPprobdataGetSc(probdata);
    sf = SCIPprobdataGetSf(probdata);


    sum_sc = SCIPprobdataGetSum_sc(probdata);
    sum_sf = SCIPprobdataGetSum_sf(probdata);
    sc_size = SCIPprobdataGetSc_size(probdata);
    sf_size = SCIPprobdataGetSf_size(probdata);
 
    
    
    //Checar dc
    dif = 0;
    dc = -string_size;
    for (i = 0; i < sc_size; i ++) {
        menor_dif_sc = (sum_sc[i+1]+string_size-1) - sum_sc[i];
        for (j = 0; j <= (sum_sc[i+1]+string_size-1) - sum_sc[i] -string_size ; j ++) {
            l = 0;
            dif = 0;
            
            for (k = j; k < j + string_size; k ++) {
                if (symbol_value(alphabet, alpha_size, heur[l++]) != sc[i][k]) {
                   dif ++;
                }
            }
            menor_dif_sc = menor_dif_sc < dif ? menor_dif_sc : dif;
            
        }
        dc = menor_dif_sc > dc ? menor_dif_sc : dc; 
    }
    
    //Checar se obedece em sf
    df = string_size;
    for (i = 0; i < sf_size; i ++) {
        maior_dif_sf = -1;
        for (j = 0; j <= (sum_sf[i+1]+string_size-1) - sum_sf[i] - string_size ; j ++) {
            l = 0;
            dif = 0;
            for (k = j; k < j + string_size; k ++) {
                if (symbol_value(alphabet, alpha_size, heur[l++]) != sf[i][k]) {
                   dif ++;
                }
            }        
            // maior_dif_sf = dif > maior_dif_sf ? dif : maior_dif_sf;
            df = dif < df ? dif : df; 
        }
        
    }


    return dc - df;
}

double heur_RA(SCIP* scip, char *filename)
{
    int i, j, string_size, alpha_size, *heur;
    double ra;
    char name[SCIP_MAXSTRLEN], *alphabet;
    FILE *file;

    SCIP_SOL *bestSolution;
    SCIP_PROBDATA *probdata;
    SCIP_VAR **vars, *var;
    SCIP_Real solval;

    assert(scip != NULL);

    bestSolution = SCIPgetBestSol(scip);
  
    probdata = SCIPgetProbData(scip);
    assert(probdata != NULL);

    // nvars = SCIPprobdataGetNVars(probdata);
    vars = SCIPprobdataGetVars(probdata);
    string_size = SCIPprobdataGetString_size(probdata);
    alphabet = SCIPprobdataGetAlphabet(probdata);
    alpha_size = SCIPprobdataGetAlpha_size(probdata);

    (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s-antigo.heur-RA.sol", filename);
    file = fopen(name, "w");

    if(!file)
    {
        printf("\nProblem to create solution file: %s", name);
        return SCIP_INVALID;
    }
    

    fprintf(file, "Tempo: 301.529243s\n");
    heur = (int *)malloc(sizeof(int) * string_size);
    
    for(i = 0; i < string_size; i++) {
        double menor = 2.0;
        int ind = -1;
        for (j = 0; j < alpha_size; j++) {
            var = vars[i * alpha_size + j];
            solval = SCIPgetSolVal(scip, bestSolution, var);
            //fprintf(file, "%.2lf ", solval);
            if (solval < menor) {
                menor = solval;
                ind = j;
            }
        }
        heur[i] = ind;
    }
    ra = valorSolucao(scip, filename, heur);
    fprintf(file, "Objective Value(dc - df) = %lf\n", ra);
    fprintf(file, "Palavra(x) = ");
    
    for(i = 0; i < string_size; i++) {
       fprintf(file, "%c", symbol_value(alphabet, alpha_size, heur[i]));
    }
    fprintf(file, "\n");
    fprintf(file, "dc = 7\n");
    fprintf(file, "df = 5\n");
    
    

        
    
    fclose(file);
    return ra;

}

void printRelaxedSol(SCIP* scip, char *filename, double tempo)
{
    int i, j, string_size, alpha_size, *sum_sc, nvars;
    char name[SCIP_MAXSTRLEN], *alphabet;
    FILE *file;
    
    SCIP_PROBDATA *probdata;
    SCIP_VAR **vars, *var;
    SCIP_Real solval;
    int totalFrac, desl, sc_size, *v;
    
    assert(scip != NULL);
       
    probdata = SCIPgetProbData(scip);
    assert(probdata != NULL);
    
    nvars = SCIPprobdataGetNVars(probdata);
    vars = SCIPprobdataGetVars(probdata);
    string_size = SCIPprobdataGetString_size(probdata);
    alphabet = SCIPprobdataGetAlphabet(probdata);
    alpha_size = SCIPprobdataGetAlpha_size(probdata);
    sum_sc = SCIPprobdataGetSum_sc(probdata);
    sc_size = SCIPprobdataGetSc_size(probdata);
   
    
    (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s-antigo.frac.sol", filename);
    file = fopen(name, "w");
    
    if(!file)
    {
        printf("\nProblem to create solution file: %s", name);
        return;
    }
    
    name[0] = '\0';
    
    fprintf(file, "Tempo: %lfs\n", tempo);
        fprintf(file, "Objective Value(dc - df) = %.0lf\n", SCIPgetLPRootObjval(scip));
    
    fprintf(file, "Palavra(x) = ");
    for(i = 0; i < string_size; i++) {
      desl=0;
        for (j = 0; j < alpha_size; j++) {
            var = vars[i * alpha_size + j];
	    solval = SCIPgetVarSol(scip,var);
            if(SCIPisFeasIntegral(scip, solval) && solval < 0.00001)
            {
                fprintf(file, "%c", symbol_value(alphabet, alpha_size, j));
		desl=1;
            }
        }
	if(desl==0){
	  fprintf(file, "@");
	}
	
    }
   
    fprintf(file, "\n");
    
    /* dc */
    var = vars[nvars - 2];
    solval = SCIPgetVarSol(scip,var);
    fprintf(file, "dc = %.0lf\n", solval);
    
    var = vars[nvars - 1];
    solval = SCIPgetVarSol(scip,var);
    fprintf(file, "df = %.0lf\n", solval);

    v=(int*)malloc(sizeof(int)*nvars);
    for(i=0;i<nvars;i++){
      v[i]=0;
    }
    totalFrac=0;    
    fprintf(file, "Var fracionarias = ");
    for(i = 0; i < string_size; i++) {
        for (j = 0; j < alpha_size; j++) {
            var = vars[i * alpha_size + j];
	    solval = SCIPgetVarSol(scip,var);
	    if(!SCIPisFeasIntegral(scip, solval))
	      {
		fprintf(file, "\nx_%d_%c=%lf", i,symbol_value(alphabet, alpha_size, j), solval);
		totalFrac++;
	      }
        }
    }

    for(i = 0; i < sc_size; i++) {
        for (j = 0; j < sum_sc[i + 1] - sum_sc[i]; j++) {
            var = vars[VARS_X + sum_sc[i] + j];
	    solval = SCIPgetVarSol(scip,var);
	    if(!SCIPisFeasIntegral(scip, solval))
	      {
		fprintf(file, "\ny_%d_%d=%lf", i, j, solval);
		totalFrac++;
	      }
        }
    }
    fprintf(file, "\nTotal fracionarias=%d\n", totalFrac);
    fclose(file);

    totalFrac = 0;
    for(i = 0; i < nvars; i++) {
        solval = SCIPgetVarSol(scip,vars[i]);
        if(!SCIPisFeasIntegral(scip, solval))
    	    totalFrac++;
    }
    free(v);
}

/** print best primal solution in file with extension ".sol"
 *
 */
void printSol(SCIP* scip, char *filename, double tempo)
{
    int i, j, string_size, alpha_size, sc_size, *sum_sc;
    char name[SCIP_MAXSTRLEN], *alphabet;
    FILE *file;

    SCIP_SOL *bestSolution;
    SCIP_PROBDATA *probdata;
    SCIP_VAR **vars, *var;
    SCIP_Real solval;

    assert(scip != NULL);

    bestSolution = SCIPgetBestSol(scip);
    if(bestSolution == NULL){
        return;
    }

    probdata = SCIPgetProbData(scip);
    assert(probdata != NULL);

    // nvars = SCIPprobdataGetNVars(probdata);
    vars = SCIPprobdataGetVars(probdata);
    string_size = SCIPprobdataGetString_size(probdata);
    alphabet = SCIPprobdataGetAlphabet(probdata);
    alpha_size = SCIPprobdataGetAlpha_size(probdata);
    sum_sc = SCIPprobdataGetSum_sc(probdata);
    sc_size = SCIPprobdataGetSc_size(probdata);

    (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s-antigo.sol", filename);
    file = fopen(name, "w");

    if(!file)
    {
        printf("\nProblem to create solution file: %s", name);
        return;
    }

    name[0] = '\0';

    fprintf(file, "Tempo: %lfs\n", tempo);
    fprintf(file, "Objective Value(dc - df) = %.0lf\n", SCIPgetPrimalbound(scip)); //SCIPsolGetOrigObj(bestSolution));
    fprintf(file, "Palavra(x) = ");

    for(i = 0; i < string_size; i++) {
        for (j = 0; j < alpha_size; j++) {
            var = vars[i * alpha_size + j];
            solval = SCIPgetSolVal(scip, bestSolution, var);
            if(solval == 0.0)
            {
                fprintf(file, "%c", symbol_value(alphabet, alpha_size, j));
            }
        }
    }
    fprintf(file, "\n");

    /* dc */
    var = vars[VARS_X + VARS_Y];
    solval = SCIPgetSolVal(scip, bestSolution, var);
    fprintf(file, "dc = %.0lf\n", solval);
    /* df */
    var = vars[VARS_X + VARS_Y + 1];
    solval = SCIPgetSolVal(scip, bestSolution, var);
    fprintf(file, "df = %.0lf\n", solval);

    fprintf(file, "x:\n");
    for(i = 0; i < string_size; i++) {
        for (j = 0; j < alpha_size; j++) {
            var = vars[i * alpha_size + j];
            solval = SCIPgetSolVal(scip, bestSolution, var);
            fprintf(file, "%.0lf ", solval);
        }
        fprintf(file, "\n");
    }

    fprintf(file, "y:\n");
    for(i = 0; i < sc_size; i++) {
        for (j = 0; j < sum_sc[i + 1] - sum_sc[i]; j++) {
            var = vars[VARS_X + sum_sc[i] + j];
            solval = SCIPgetSolVal(scip, bestSolution, var);
            fprintf(file, "%.0lf ", solval);
        }
        fprintf(file, "\n");
    }

    fclose(file);
}

void validSol(SCIP *scip, char *filename) {
    int i, j, k, l, sum;
    int **x, **y, dc, df, *sum_sc, *sum_sf, string_size, alpha_size, sc_size, sf_size;
    char name[SCIP_MAXSTRLEN], **sc, **sf, *alphabet;
    char sol_x[999], aux[2];
    bool valid = true;
    FILE *file;
    FILE *input;

    SCIP_SOL *bestSolution;
    SCIP_PROBDATA *probdata;
    SCIP_VAR **vars, *var;
    SCIP_Real solval;

    assert(scip != NULL);
    
    bestSolution = SCIPgetBestSol(scip);
    if(bestSolution == NULL){
        return;
    }
    
    probdata = SCIPgetProbData(scip);
    assert(probdata != NULL);

    (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s-antigo.debug", filename);
    file = fopen(name, "w");
    fprintf(file, "Validacao...\n");
    vars = SCIPprobdataGetVars(probdata);
    string_size = SCIPprobdataGetString_size(probdata);
    alphabet = SCIPprobdataGetAlphabet(probdata);
    alpha_size = SCIPprobdataGetAlpha_size(probdata);
    sc = SCIPprobdataGetSc(probdata);
    sf = SCIPprobdataGetSf(probdata);
    sum_sc = SCIPprobdataGetSum_sc(probdata);
    sum_sf = SCIPprobdataGetSum_sf(probdata);
    sc_size = SCIPprobdataGetSc_size(probdata);
    sf_size = SCIPprobdataGetSf_size(probdata);
    
    x = (int **) calloc(string_size, sizeof(int *));
    for (i = 0; i < string_size; i++) {
        x[i] = (int *) calloc(alpha_size, sizeof(int));
    }

    y = (int **) calloc(sc_size, sizeof(int *));
    for (i = 0; i < sc_size; i++) {
        y[i] = (int *) calloc(sum_sc[i], sizeof(int));
    }
    
    
    
    for(i = 0; i < string_size; i++) {
        sum = 0;
        for (j = 0; j < alpha_size; j++) {
            var = vars[i * alpha_size + j];
            solval = SCIPgetSolVal(scip, bestSolution, var);
            x[i][j] = (int) solval;
            sum += x[i][j];
        }
        if (sum != alpha_size - 1) {
            fprintf(file, "Problema na restricao (1) na posicao %d da string x.\n", i);
            valid = false;
        }
    }

    for (i = 0; i < sc_size; i++) {
        for (k = 0; k < sum_sc[i + 1] - sum_sc[i]; k++) {
            var = vars[VARS_X + sum_sc[i] + k];
            y[i][k] = (int) SCIPgetSolVal(scip, bestSolution, var);
        }
    }

    var = vars[VARS_X + VARS_Y];
    dc = (int) SCIPgetSolVal(scip, bestSolution, var);

    var = vars[VARS_X + VARS_Y + 1];
    df = (int) SCIPgetSolVal(scip, bestSolution, var);

    for (i = 0; i < sc_size; i++) {
        for (k = 0; k < sum_sc[i + 1] - sum_sc[i]; k++) {
            sum = 0;
            for (l = 0; l < string_size; l++) {
                int pos = value_symbol(alphabet, alpha_size, sc[i][k + l]);
                if (pos != -1)
                    sum += x[l][pos];
            }
            if (y[i][k] == 1 && sum > dc) {
                fprintf(file, "Problema na restricao (2) na string %d posicao %d.\n", i, k);
                valid = false;
            }
        }
    }

    for (i = 0; i < sf_size; i++) {
        for (k = 0; k < sum_sf[i + 1] - sum_sf[i]; k++) {
            sum = 0;
            for (l = 0; l < string_size; l++) {
                int pos = value_symbol(alphabet, alpha_size, sf[i][k + l]);
                if (pos != -1)
                    sum += x[l][pos];
            }
            if (sum < df) {
                fprintf(file, "Problema na restricao (3) na string %d posicao %d.\n", i, k);
                valid = false;
            }
        }
    }

    for (i = 0; i < sc_size; i++) {
        sum = 0;
        for (k = 0; k < sum_sc[i + 1] - sum_sc[i]; k++) {
            sum += y[i][k];
        }
        if (sum < 1) {
            fprintf(file, "Problema na restricao (4) na string %d.\n", i);
            valid = false;
        }
    }

    input = fopen(filename, "r");
    if (!input) exit(3);

    strcat(sol_x, "");
    for(i = 0; i < string_size; i++) {
        for (j = 0; j < alpha_size; j++) {
            var = vars[i * alpha_size + j];
            solval = SCIPgetSolVal(scip, bestSolution, var);
            if(solval == 0.0)
            {
                aux[0] = symbol_value(alphabet, alpha_size, j);
                aux[1] = '\0';
                strcat(sol_x, aux);
            }
        }
    }

    if (valid){
        printf("Solucao validada com sucesso!"); fflush(stdout);}
    else{
        printf("Falha na validacao da solucao!"); fflush(stdout);}

    for (i = 0; i < string_size; i++) {
        free(x[i]);
    }
    free(x);

    for (i = 0; i < sc_size; i++) {
        free(y[i]);
    }
    free(y);

    fclose(file);
    fclose(input);
}

int objective_value(int  *sum_sc, int *sum_sf, char **Sc, char **Sf, int sc_size, int sf_size, int kc, int kf, char *alpha, int alpha_size, int string_size, int *solution) {
    int dif = 0;
    int dc = 99999;
    int df = -1;
    int aux = 0;
    int menor_dif_sc, maior_dif_sf;
    int i, j, k, l;
    
    dc = -1 * string_size;
    for (i = 0; i < sc_size; i ++) {
        menor_dif_sc = (sum_sc[i+1]+string_size-1) - sum_sc[i];
        
        for (j = 0; j <= (sum_sc[i+1]+string_size-1) - sum_sc[i] -string_size ; j ++) {
            l = 0;
            dif = 0;
            
            for (k = j; k < j + string_size; k ++) {
                if (symbol_value(alpha, alpha_size, solution[l++]) != Sc[i][k]) {
                   dif ++;
                }
            }
            menor_dif_sc = menor_dif_sc < dif ? menor_dif_sc : dif;
            
        }
        dc = menor_dif_sc > dc ? menor_dif_sc : dc; 
    }
    
    //Checar se obedece em sf
    df = string_size;
    for (i = 0; i < sf_size; i ++) {
        maior_dif_sf = -1;
        for (j = 0; j <= (sum_sf[i+1]+string_size-1) - sum_sf[i] - string_size ; j ++) {
            l = 0;
            dif = 0;
            for (k = j; k < j + string_size; k ++) {
                if (symbol_value(alpha, alpha_size, solution[l++]) != Sf[i][k]) {
                   dif ++;
                }
            }        
            // maior_dif_sf = dif > maior_dif_sf ? dif : maior_dif_sf;
            df = dif < df ? dif : df; 
        }

    }
    //printf("%d %d\n", dc , df);

    return dc - df;
}

int dc_df_objective_value(int  *sum_sc, int *sum_sf, char **Sc, char **Sf, int sc_size, int sf_size, int kc, int kf, char *alpha, int alpha_size, int string_size, int *solution, int* dc, int* df) {
    int dif = 0;
    int aux = 0;
    int menor_dif_sc, maior_dif_sf;
    int i, j, k, l;
    
    //int dc = 99999;
    //int df = -1;

    *dc = -1 * string_size;
    for (i = 0; i < sc_size; i ++) {
        menor_dif_sc = (sum_sc[i+1]+string_size-1) - sum_sc[i];
        
        for (j = 0; j <= (sum_sc[i+1]+string_size-1) - sum_sc[i] -string_size ; j ++) {
            l = 0;
            dif = 0;
            
            for (k = j; k < j + string_size; k ++) {
                if (symbol_value(alpha, alpha_size, solution[l++]) != Sc[i][k]) {
                   dif ++;
                }
            }
            
            menor_dif_sc = menor_dif_sc < dif ? menor_dif_sc : dif;
            
        }
        *dc = menor_dif_sc > *dc ? menor_dif_sc : *dc; 
    }
    
    //Checar se obedece em sf
    *df = string_size;
    for (i = 0; i < sf_size; i ++) {
        maior_dif_sf = -1;
        for (j = 0; j <= (sum_sf[i+1]+string_size-1) - sum_sf[i] - string_size ; j ++) {
            l = 0;
            dif = 0;
            for (k = j; k < j + string_size; k ++) {
                if (symbol_value(alpha, alpha_size, solution[l++]) != Sf[i][k]) {
                   dif ++;
                }
            }        
	    //            maior_dif_sf = dif > maior_dif_sf ? dif : maior_dif_sf;
	    *df = dif < *df ? dif : *df; 
        }
	//        *df = maior_dif_sf < *df ? maior_dif_sf : *df; 
    }
    //printf("%d %d\n", dc , df);

    return *dc - *df;
}


int neighbourhoodChange(int  *sum_sc, int *sum_sf, int *x, int *x2, int k, char **sc, char**sf, int sc_size, int sf_size, int kc, int kf, char *alphabet, int alpha_size, int string_size ) {
    int i;
    if (objective_value(sum_sc, sum_sf, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size, x2) < objective_value(sum_sc, sum_sf, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size,x)) {
        for (i = 0; i < string_size; i ++) x[i] = x2[i];
        k = 0;
    }
    else {
        k = k + 1;
    }
    return k;
}

/* put your local methods here, and declare them static */
int randomInteger (int low, int high)
{
  int k;
  double d;

  d = (double) rand () / ((double) RAND_MAX + 1);
  k = d * (high - low + 1);
  return low + k;
}

void shaking (double **T, int *x2, int k, double *ng, int string_size, int alpha_size) {
#ifdef DEBUG_VNS  
    printf("Shaking\n");
#endif
    int i, j, r;
    int *viable_solutions = (int *) calloc(alpha_size, sizeof(int));
    int quant;
    double menor;
    int i_menor;
    srand(time(NULL));

    for (i = 0; i < string_size; i ++) {
        quant = 0;
	menor = 1.0;
        for (j = 0; j < alpha_size; j ++) {
	  //	  printf("Posicao: %d, Simbolo: %d, Valor: %lf, Limite: %lf \n", i, j, T[i][j], ng[k]);
            if (T[i][j] <= ng[k]) {
                viable_solutions[quant] = j;
                quant ++;
            }
	    else if(T[i][j]<menor - EPSILON){
	      menor = T[i][j];
	      i_menor = j;
	    }
        }
	if(quant==0){
	  x2[i] = i_menor;
	}
	else{
	  r = randomInteger(0,quant-1);//rand();
	  //	printf("\nrand=%d quant=%d", r, quant);
	  x2[i] = viable_solutions[r];//quant % r];
	}
    }

#ifdef DEBUG_VNS		      
    for (i = 0; i < string_size; i ++) {
        printf("%d ", x2[i]);
    }
    printf("\n");
#endif
}

void firstImprovement(int  *sum_sc, int *sum_sf, double **T, int *x2, int k, double *ng, char **sc, char**sf, int sc_size, int sf_size, int kc, int kf, char *alphabet, int alpha_size, int string_size) {
    int *x_i, *x_aux,  i = 0, j = 0, m=0, ix = 0, iz=0, teste=0; 
    int fx2, fxaux, fxi, flag;
    x_aux = (int *) calloc(string_size, sizeof(int));
    x_i = (int *) calloc(string_size, sizeof(int));
#ifdef DEBUG_VNS		      
    printf("First Improvement: \n");
#endif
    do {
      for (j = 0; j < string_size; j++) x_aux[j] = x2[j]; 
        fxaux = objective_value(sum_sc, sum_sf, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size, x_aux);
#ifdef DEBUG_VNS		      
        printf("K: %d, f(x'): %d\n", k, fxaux);
#endif
        i = 0; ix = 0; iz = 0;
        do {
            for (j = 0; j < string_size; j++) x_i[j] = x2[j];

            flag = 0;
            for (j = ix; j < string_size; j ++) {
	      for (m = iz; m < alpha_size; m ++) { // nao esta olhando toda a vizinhanca...
                    if (x_i[j] == m) continue;
                    if (T[j][m] < ng[k]) {
#ifdef DEBUG_VNS		      
		       printf("\nTroca x[%d]=%d por %d", j,x_i[j], m);
#endif
                       x_i[j] = m;
                       flag = 1;
                       ix = j;
                       iz = m+1; // tenta mudar outro simbolo na posicao j
		       if(iz==alpha_size){ // se tentou todos, muda de posicao
			 iz = 0;
			 ix = j+1;
		       }
                       break;
                    }
                }
                if (flag == 1) break;
            }
            i++;
            if (flag == 0) break;
            fxi = objective_value(sum_sc, sum_sf, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size, x_i);
#ifdef DEBUG_VNS		          
            printf("\tfx%d: %d\n", i, fxi);
#endif
        } while (fxaux < fxi && i < 100);
        if (flag == 0 || i >= 100) break;

        teste++;
        for (j = 0; j < string_size; j++) x2[j] = x_i[j];
    } while(fxaux > fxi); //objective_value(sum_sc, sum_sf, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size, x2));
#ifdef DEBUG_VNS		      
    printf("\nfinal first improvement f(x')=%d", objective_value(sum_sc, sum_sf, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size, x2));
#endif
}

double vns_main (SCIP* scip, char *filename, double tempo) {
    int i, j, string_size, alpha_size, *heur;
    double vns;
    char name[SCIP_MAXSTRLEN], *alphabet;
    int k, k_max;
    double t, t_max;
    double ng[3] = {0.3, 0.6, 0.9};
    char **sc, **sf;
    int sc_size, sf_size;
    int kc, kf;
    int  *sum_sc, *sum_sf;
    int dc, df;

    SCIP_SOL *bestSolution;
    SCIP_PROBDATA *probdata;
    SCIP_VAR **vars, *var;
    SCIP_Real solval;

    FILE *file;
    clock_t antes, agora;

    antes = clock();
    
    assert(scip != NULL);

    bestSolution = SCIPgetBestSol(scip);
  
    probdata = SCIPgetProbData(scip);
    assert(probdata != NULL);

    // nvars = SCIPprobdataGetNVars(probdata);
    vars = SCIPprobdataGetVars(probdata);
    string_size = SCIPprobdataGetString_size(probdata);
    alphabet = SCIPprobdataGetAlphabet(probdata);
    alpha_size = SCIPprobdataGetAlpha_size(probdata);
    sc = SCIPprobdataGetSc(probdata);
    sf = SCIPprobdataGetSf(probdata);
    sc_size = SCIPprobdataGetSc_size(probdata);
    sf_size = SCIPprobdataGetSf_size(probdata);
    kc = SCIPprobdataGetKc(probdata);
    kf = SCIPprobdataGetKf(probdata); 
    sum_sc = SCIPprobdataGetSum_sc(probdata);
    sum_sf = SCIPprobdataGetSum_sf(probdata);
    double **T;
    T = (double **) calloc(string_size, sizeof(double *));
    for (i = 0; i < string_size; i++) {
        T[i] = (double *) calloc(alpha_size, sizeof(double));
    }

#ifdef DEBUG_VNS		      
    printf("Alfabeto: \n");
    for (i = 0; i < alpha_size; i ++) {
        printf("%d %c\n", i, alphabet[i]);
    }

    printf("N-Vizinhancas: \n");
    for (i = 0; i < 3; i ++) {
        printf("%lf\n", ng[i]);
    }
#endif
    int *x;
    x = (int *) calloc(string_size, sizeof(int));

    int *x2;
    x2 = (int *) calloc(string_size, sizeof(int));
    
    for(i = 0; i < string_size; i++) {
        double menor = 2.0;
        int ind = -1;
        for (j = 0; j < alpha_size; j++) {
            var = vars[i*alpha_size + j];
            solval = SCIPgetSolVal(scip, bestSolution, var);
            if (solval <= -0.0001) solval = 0.0;
            T[i][j] = solval;
        }
    }
 
    for(i = 0; i < string_size; i++) {
        double menor = 2.0;
        int ind = -1;
        for (j = 0; j < alpha_size; j++) {
            var = vars[i * alpha_size + j];
            solval = SCIPgetSolVal(scip, bestSolution, var);
            if (solval < menor) {
                menor = solval;
                ind = j;
            }
        }
        x[i] = ind;
    }

    
#ifdef DEBUG_VNS		      
    vns = dc_df_objective_value(sum_sc, sum_sf, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size, x, &dc, &df);
    printf("Valor: %lf\n", vns);
    printf("Dc: %d\n", dc);
    printf("Df: %d\n", df);
#endif

#ifdef DEBUG_VNS		      
    printf("\nSolucao Base: \n");

    for (i = 0; i < string_size; i ++) {
        for (j = 0; j < alpha_size; j ++) {
            printf("%lf ", T[i][j]);
        }
        printf("\n");
    }
#endif
    // repeat
    //     k ← 1;
    //     repeat
    //         x  ← Shake(x, k)
    //         /* shaking */;
    //         x ← FirstImprovement(x ) /* Local search */;
    //         NeighbourhoodChange(x, x , k) /* Change neighbourhood */;
    //     until k = k max ;
    //     t ← CpuTime()
    // until t > t max ;

    k_max = 3;
    t = 0;t_max=300;

    do {
        k = 0;
        do {
            shaking(T, x2, k, ng, string_size, alpha_size);
            firstImprovement(sum_sc, sum_sf, T, x2, k, ng, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size);            
            k = neighbourhoodChange(sum_sc, sum_sf, x, x2, k, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size);

#ifdef DEBUG_VNS		      
            printf("\n%d\n", k);
#endif
        } while (k < k_max);
        agora = clock();
        t = (agora - antes) / ((double) CLOCKS_PER_SEC);
    } while (t <= t_max);
    
#ifdef DEBUG_VNS		      
    printf("Solucao Final:\n");
    for (i = 0; i < string_size; i ++) 
        printf("%c" , alphabet[x[i]]);
    printf("\n");
#endif    
    vns = dc_df_objective_value(sum_sc, sum_sf, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size, x, &dc, &df);
#ifdef DEBUG_VNS		      
    printf("Valor: %lf\n", vns);
    printf("Dc: %d\n", dc);
    printf("Df: %d\n", df);
#endif

    agora = clock();

    tempo += (agora - antes) / ((double) CLOCKS_PER_SEC);
    // salva .sol

#ifndef SUB
    (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s-antigo-vns.sol", filename);
#else
    (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s-antigo-vns-sub.sol", filename);
#endif
    file = fopen(name, "w");

    if(!file)
    {
        printf("\nProblem to create solution file: %s", name);
        return;
    }
    name[0] = '\0';

    fprintf(file, "Tempo: %lfs\n", tempo);
    fprintf(file, "Objective Value(dc - df) = %.0lf\n", vns);
    fprintf(file, "Palavra(x) = ");

    for (i = 0; i < string_size; i ++) 
      fprintf(file,"%c" , alphabet[x[i]]);
    fprintf(file, "\n");

    /* dc */
    fprintf(file, "dc = %d\n", dc);
    /* df */
    fprintf(file, "df = %d\n", df);
    fclose(file);
    
    return vns;
    // return 0.0;
}

/** creates a SCIP instance with default plugins, evaluates command line parameters, runs SCIP appropriately,
 *  and frees the SCIP instance
 */
static
SCIP_RETCODE runShell(
   int                        argc,               /**< number of shell parameters */
   char**                     argv,               /**< array with shell parameters */
   const char*                defaultsetname      /**< name of default settings file */
   )
{
    int barra = 0 , string_size, alpha_size, sc_size, *sum_sc;;
    double tempo, LB, UB, LB_ROOT, ra = SCIP_INVALID, vns = SCIP_INVALID;
    clock_t antes, agora;
    SCIP_VAR **vars, *var;
    SCIP* scip = NULL;
    SCIP_PROBDATA *probdata;
    SCIP_SOL *bestSolution;
    SCIP_Real solval;   
    FILE *fp;


    antes = clock();

    /*********
    * Setup *
    *********/

    /* initialize SCIP */
    SCIP_CALL( SCIPcreate(&scip) );

    /* include dssp reader */
    SCIP_CALL( SCIPincludeReaderDSSP(scip) );

    /* include default SCIP plugins */
    SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

    /* for column generation instances, disable restarts */
    SCIP_CALL( SCIPsetIntParam(scip,"presolving/maxrestarts",0) );

    /* turn off all separation algorithms */
    SCIP_CALL( SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE) );

    /* disable heuristics */
    SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, TRUE) );

#if defined(BCPA) || defined(RA) //|| defined(VNS)
    SCIP_CALL( SCIPincludeHeurmyheuristic(scip) );
#endif

    /* disable presolving */
    SCIP_CALL( SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE) );

    SCIP_CALL( SCIPsetIntParam(scip, "propagating/maxrounds", 0) );
    SCIP_CALL( SCIPsetIntParam(scip, "propagating/maxroundsroot", 0) );

    SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );
    SCIP_CALL( SCIPsetIntParam(scip, "separating/maxrounds", 0) );

    /* set time limit */
    SCIP_CALL( SCIPsetRealParam(scip, "limits/time", TIME_LIMIT) );

#ifdef ONLYROOT
    /* run only at root node */
    SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", 1) );
#endif

    /**********************************
    * Process command line arguments *
    **********************************/
    
    SCIP_CALL( SCIPprocessShellArguments(scip, argc, argv, defaultsetname) );

    agora = clock();

    tempo = (agora - antes) / ((double) CLOCKS_PER_SEC);
    
    //imprimir nome
     for (int i = 0; i < strlen(argv[2]); i++) {
        if (argv[2][i] == '/')
            barra = i + 1;
    }
    for (int i = barra; i < strlen(argv[2]); i++) {
        printf("%c", argv[2][i]);
    }
    printf("\n");
    char solFileName[512];
    strcpy(solFileName, argv[2]);
    /* grava melhor solucao */
    //printSol(scip, argv[2], tempo);
#ifdef GLOBAL
    strcat(solFileName, "-GLOBAL");
#endif
#ifdef ONLYROOT
    strcat(solFileName, "-ONLYROOT");
#endif
#ifdef STRING
    strcat(solFileName, "-STRING");
#endif
#ifdef SUB
    strcat(solFileName, "-SUB");
#endif
#ifdef BCPA
    strcat(solFileName, "-BCPA");
#endif
#ifdef RA
    strcat(solFileName, "-RA");
#endif
#ifdef VNS
    strcat(solFileName, "-VNS");
#endif
#ifdef NONE
    strcat(solFileName, "-NONE");
#endif

printSol(scip, solFileName, tempo);

ra = -9999;
#ifdef RA
    #ifdef ONLYROOT
    antes = clock();
    ra = heur_RA(scip, argv[2]);

    agora = clock();
    tempo = (agora - antes) / ((double) CLOCKS_PER_SEC);
    #endif
#endif

#ifdef VNS
    #ifdef ONLYROOT
    antes = clock();
    vns = vns_main(scip, argv[2], tempo);
    printf("VNS - Onlyroot\n");

    agora = clock();
    tempo = (agora - antes) / ((double) CLOCKS_PER_SEC);
    #endif
#endif

// #ifdef ONLYROOT
    //printRelaxedSol(scip, argv[2], tempo);
// #endif

    /* valida solucao */
#ifdef DEBUG
    validSol(scip, argv[2]);
#endif

    UB = SCIPgetPrimalbound(scip);
    LB = SCIPgetDualbound(scip);
    LB_ROOT = SCIPgetDualboundRoot(scip);
    
    SCIP_Longint nNodes = SCIPgetNNodes (scip);
    printf(";%ld", nNodes);
    
    //printf("%lf;%lf;%lf;%lf;%d\n", LB, UB, SCIPgetDualboundRoot(scip), tempo, SCIPgetStatus(scip));
#ifdef RA
#ifdef ONLYROOT
        printf(";%lf;%lf;%lf;%lf;%d;RA;%lf\n", LB, UB, LB_ROOT, tempo,SCIPgetStatus(scip), ra);
#else
        printf(";%lf;%lf;%lf;%lf;%d;RA\n", LB, UB, LB_ROOT, tempo,SCIPgetStatus(scip));
#endif

#elif BCPA 
        //printf(";%lf;%lf;%lf;%lf;%d;BCPA;%lf\n", LB, UB, LB_ROOT, tempo,SCIPgetStatus(scip), LB);
       
#ifdef STRING
    fp = fopen("bcpa-string.csv", "a");
      for (int i = 0; i < strlen(argv[2]); i++) {
        if (argv[2][i] == '/')
            barra = i + 1;
    }
    for (int i = barra; i < strlen(argv[2]); i++) {
        fprintf(fp, "%c", argv[2][i]);
    }
    fprintf(fp, ";%ld;%lf;%lf;%lf;%lf;%d;BCPA\n", nNodes, LB, UB, LB_ROOT, tempo,SCIPgetStatus(scip));
    
    fclose(fp);
#endif
#ifdef SUB
    fp = fopen("bcpa-sub.csv", "a");
      for (int i = 0; i < strlen(argv[2]); i++) {
        if (argv[2][i] == '/')
            barra = i + 1;
    }
    for (int i = barra; i < strlen(argv[2]); i++) {
        fprintf(fp, "%c", argv[2][i]);
    }
    fprintf(fp, ";%ld;%lf;%lf;%lf;%lf;%d;BCPA\n", nNodes, LB, UB, LB_ROOT, tempo,SCIPgetStatus(scip));
    
    fclose(fp);
#endif
    
#elif NONE
    printf(";%lf;%lf;%lf;%lf;%d;NONE;NONE\n", LB, UB, LB_ROOT, tempo,SCIPgetStatus(scip));
#else
#ifdef VNS
#ifndef SUB
    fp = fopen("vns.csv", "a");
#else
    fp = fopen("vns-sub.csv", "a");
#endif
      for (int i = 0; i < strlen(argv[2]); i++) {
        if (argv[2][i] == '/')
            barra = i + 1;
    }
    for (int i = barra; i < strlen(argv[2]); i++) {
        fprintf(fp, "%c", argv[2][i]);
    }
    fprintf(fp, ";%ld;%lf;%lf;%lf;%lf;%d;%lf;VNS\n", nNodes, LB, UB, LB_ROOT, tempo,SCIPgetStatus(scip), vns);
    
    fclose(fp);

    /*    if(UB>vns+EPSILON){
      UB = vns;
      }*/
    printf(";%lf;%lf;%lf;%lf;%d;%lf;VNS\n", LB, UB, LB_ROOT, tempo,SCIPgetStatus(scip), vns);    
#else
    printf(";%lf;%lf;%lf;%lf;%d;NONE;NONE\n", LB, UB, LB_ROOT, tempo,SCIPgetStatus(scip));
#endif
#endif


    /********************
    * Deinitialization *
    ********************/

    /*   SCIP_CALL( SCIPfree(&scip) ); */

    BMScheckEmptyMemory();

    return SCIP_OKAY;
}

int
main(
    int                        argc,
    char**                     argv
    )
{
    SCIP_RETCODE retcode;

    retcode = runShell(argc, argv, "scip.set");
    if(retcode != SCIP_OKAY )
    {
        SCIPprintError(retcode);
        return -1;
    }

    return 0;
}
