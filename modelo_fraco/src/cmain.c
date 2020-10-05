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
#include "reader_csp.h"
#include "parameters_dssp.h"
#include "heur_myheuristic.h"
#include "cons_xyz.h"



void printSol(SCIP *, char *, double);
void validSol(SCIP *, char *);
double heur_RA(SCIP* , char *);
double heur_BCPA(SCIP* , char *);
int valorSolucao(SCIP*, char *, int *);
void printRelaxedSol(SCIP *, char *, double);
SCIP_RETCODE setupProblem(
    SCIP*                 scip                /**< SCIP data structure */,
    SCIP*                 scip_o                /**< SCIP data structure */
    );


int valorSolucao(SCIP *scip, char *filename, int *heur){
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
    if (sc_size > 0) {
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
    }
    else
      dc = 0;
    
    //Checar se obedece em sf
    if (sf_size > 0) {
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
    }
    else
      df = 0;


    return dc - df;
}

// double heur_bcpa (SCIP* scip, char *filename) {
//     int i, j, k, l, string_size, alpha_size, *heur;
//     double ra;
//     int bcpa;
    
//     FILE *file;
//     int         sc_size;
//     int*        sum_sc;
//     char name[SCIP_MAXSTRLEN], **sc, *alphabet;

//     SCIP_SOL *bestSolution;
//     SCIP_PROBDATA *probdata;
//     SCIP_VAR **vars, *var;
//     SCIP_Real solval;

//     assert(scip != NULL);
    
//     bestSolution = SCIPgetBestSol(scip);
  
//     probdata = SCIPgetProbData(scip);
//     assert(probdata != NULL);

//     // nvars = SCIPprobdataGetNVars(probdata);
//     vars = SCIPprobdataGetVars(probdata);
//     string_size = SCIPprobdataGetString_size(probdata);
//     alphabet = SCIPprobdataGetAlphabet(probdata);
//     alpha_size = SCIPprobdataGetAlpha_size(probdata);
//     sc_size = SCIPprobdataGetSc_size(probdata) == 0 ? SCIPprobdataGetSf_size(probdata) : SCIPprobdataGetSc_size(probdata);
//     sum_sc = SCIPprobdataGetSum_sc(probdata);
//     sc = SCIPprobdataGetSc(probdata);

//     #ifndef NOVO  
//     (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s-antigo.heur-BCPA.sol", filename);
// #else
//     (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s-novo.heur-BCPA.sol", filename);
// #endif
//     file = fopen(name, "w");

//     if(!file)
//     {
//         printf("\nProblem to create solution file: %s", name);
//         return SCIP_INVALID;
//     }

//     SCIP*               scip_b;
//     SCIP_PROBDATA* probdata_b;
//     SCIP_CONS** conss_b;
//     SCIP_VAR** vars_b, *var_b;
//     char *probname;

//     int ncons, nvars;

//     SCIP_CALL( SCIPcreate(&scip_b) );

//     /* include dssp reader */
//     SCIP_CALL( SCIPincludeReaderDSSP(scip_b) );
//     SCIP_CALL( SCIPincludeReadercsp(scip_b) );


//     /* include default SCIP plugins */
//     SCIP_CALL( SCIPincludeDefaultPlugins(scip)_b );

//     /* for column generation instances, disable restarts */
//     SCIP_CALL( SCIPsetIntParam(scip_b,"presolving/maxrestarts",0) );

//     /* turn off all separation algorithms */
//     SCIP_CALL( SCIPsetSeparating(scip_b, SCIP_PARAMSETTING_OFF, TRUE) );

//     /* disable heuristics */
//     SCIP_CALL( SCIPsetHeuristics(scip_b, SCIP_PARAMSETTING_OFF, TRUE) );

// // #if defined(BCPA) || defined(RA) //|| defined(VNS)
// //     SCIP_CALL( SCIPincludeHeurmyheuristic(scip) );
// // #endif

//     /* disable presolving */
//     SCIP_CALL( SCIPsetPresolving(scip_b, SCIP_PARAMSETTING_OFF, TRUE) );

//     SCIP_CALL( SCIPsetIntParam(scip_b, "propagating/maxrounds", 0) );
//     SCIP_CALL( SCIPsetIntParam(scip_n, "propagating/maxroundsroot", 0) );

//     SCIP_CALL( SCIPsetIntParam(scip_b, "presolving/maxrounds", 0) );
//     SCIP_CALL( SCIPsetIntParam(scip_b, "separating/maxrounds", 0) );

//     /* set time limit */
//     SCIP_CALL( SCIPsetRealParam(scip_b, "limits/time", TIME_LIMIT) );


 
//     assert(scip_b != NULL);

//     /* create event handler if it does not exist yet */
//     if(SCIPfindEventhdlr(scip_b, EVENTHDLR_NAME) == NULL)
//     {
//         SCIP_CALL(SCIPincludeEventhdlrBasic(scip_b, NULL, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecAddedVar, NULL));
//     }

//     /* create problem in scip_b and add non-NULL callbacks via setter functions */
//     SCIP_CALL(SCIPcreateProbBasic(scip_b, probname));

//     SCIP_CALL(SCIPsetProbDelorig(scip_b, probdelorigdssp));
//     SCIP_CALL(SCIPsetProbTrans(scip_b, probtransdssp));
//     SCIP_CALL(SCIPsetProbDeltrans(scip_b, probdeltransdssp));
//     SCIP_CALL(SCIPsetProbInitsol(scip_b, probinitsoldssp));
//     SCIP_CALL(SCIPsetProbExitsol(scip_b, probexitsoldssp));

//     /* set objective sense */
//     SCIP_CALL(SCIPsetObjsense(scip_b, SCIP_OBJSENSE_MINIMIZE));

//     /* tell scip_b that the objective will be always integral */
//     SCIP_CALL(SCIPsetObjIntegral(scip_b));

//     /* Number of constraints */
//     SCIP_CALL(SCIPallocBufferArray(scip_b, &conss_b, CONS_1 + CONS_2 + CONS_3 + CONS_4));

//     /* Number of variables */
//     SCIP_CALL(SCIPallocBufferArray(scip_b, &vars_b, VARS_X + VARS_Y + VARS_D));
   
        
    
//     /* (1) constraint - only one symbol is chosen */
//     for (i = 0; i < CONS_1; i++) {
//         SCIP_CALL(SCIPcreateConsBasicLinear (scip_b, &conss_b[i], "(1)", 0, NULL, NULL, (double) alpha_size - 1.0, (double) alpha_size - 1.0));
//         SCIP_CALL(SCIPaddCons(scip_b, conss_b[i]));
//     }

//     /* (2) constraint - is sci's symbol matched? */
//     ncons = i;
//     for (; i < ncons + CONS_2; i++) {
//         SCIP_CALL(SCIPcreateConsBasicLinear (scip_b, &conss_b[i], "(2)", 0, NULL, NULL, -SCIPinfinity(scip_b), (double) string_size)); // string_size faz o papel de M(no modelo)
//         SCIP_CALL(SCIPaddCons(scip_b, conss_b[i]));
//     }

//     /* (3) constraint - is sfi's symbol matched? */
//     ncons = i;
//     for (; i < ncons + CONS_3; i++) {
//         SCIP_CALL(SCIPcreateConsBasicLinear (scip_b, &conss_b[i], "(3)", 0, NULL, NULL,  0.0, SCIPinfinity(scip_b)));
//         SCIP_CALL(SCIPaddCons(scip_b, conss_b[i]));
//     }

//     /* (4) cosntraint - disjunctive constraint */
//     ncons = i;
//     for (; i < ncons + CONS_4; i++) {
//         SCIP_CALL(SCIPcreateConsBasicLinear (scip_b, &conss_b[i], "(4)", 0, NULL, NULL, 1.0, SCIPinfinity(scip_b)));
//         SCIP_CALL(SCIPaddCons(scip_b, conss_b[i]));
//     }
    
//     // /* (4) cosntraint - disjunctive constraint */
//     // if (definedVarsSize > 0) {
//     //     ncons = i;
//     //     for (; i < ncons + CONS_5; i++) {
//     //         SCIP_CALL(SCIPcreateConsBasicLinear (scip_b, &conss_b[i], "(5)", 0, NULL, NULL, 1.0, SCIPinfinity(scip_b)));
//     //         SCIP_CALL(SCIPaddCons(scip_b, conss_b[i]));
//     //     }
//     // }

//     nvars = 0;
//     ncons = i;
    
//     //int bcpa = 0;

//     /**
//     * add variables xij
//     */
//     for(i = 0; i < string_size; i++) {
//         for (j = 0; j < alpha_size; j++) {
//             (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x_%d_%c", i, symbol_value(alphabet, alpha_size, j));
//             SCIPdebugMessage("create variable %s\n", name);

//             /* create a basic variable object */
//             SCIP_CALL(SCIPcreateVarBasic(scip_b, &var_b, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY));
//             assert(var_b != NULL);

//             vars_b[nvars++] = var_b;

//             /* add variable to the problem */
//             SCIP_CALL(SCIPaddVar(scip_b, var_b));


//             /* add variable to corresponding constraint */
//             /* add variable to constraint (1) */
//             SCIP_CALL(SCIPaddCoefLinear(scip_b, conss_b[ROWS_1(i)], var_b, 1.0));

//             /* add variable to constraint (2) */
//             for (k = 0; k < sc_size; k++) {
//                 for (l = 0; l < sum_sc[k + 1] - sum_sc[k]; l++) {
//                     if (j == value_symbol(alphabet, alpha_size, sc[k][l + i]))
//                         SCIP_CALL(SCIPaddCoefLinear(scip_b, conss_b[ROWS_2(k, l)], var_b, 1.0));
//                 }
//             }

            
//             // if (arrayOfDefinedVars[i * alpha_size + j] != -1) {
//             //     SCIP_CALL(SCIPaddCoefLinear(scip_b, conss_b[ROWS_5(bcpa)], var_b, 1.0));
//             //     bcpa++;
//             // }
//         }
//     }

//     /**
//     * add variables yij
//     */
//     for (i = 0; i < sc_size; i++) {
//         for (j = 0; j < sum_sc[i + 1] - sum_sc[i]; j++) {
//             (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "y_%d_%d", i, j);
//             SCIPdebugMessage("create variable %s\n", name);
            
            
            
//             /* create a basic variable object */
//             SCIP_CALL(SCIPcreateVarBasic(scip_b, &var_b, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY));
//             assert(var_b != NULL);

//             vars_b[nvars++] = var_b;

//             /* add variable to the problem */
//             SCIP_CALL(SCIPaddVar(scip_b, var_b));

//             /* add variable to corresponding constraint */
//             /* add variable to constraint (2) */
//             SCIP_CALL(SCIPaddCoefLinear(scip_b, conss_b[ROWS_2(i, j)], var_b, string_size));

//             /* add variable to constraint (4) */
//             SCIP_CALL(SCIPaddCoefLinear(scip_b, conss_b[ROWS_4(i)], var_b, 1.0));
//         }
//     }

//     /**
//     * add variable dc
//     */
//     (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "d_c");
//     SCIPdebugMessage("create variable %s\n", name);

//     /* create a basic variable object */
// #ifdef FSP
//     SCIP_CALL(SCIPcreateVarBasic(scip_b, &var_b, name, 0.0, 0.0, 1.0, SCIP_VARTYPE_INTEGER));
// #else
//     SCIP_CALL(SCIPcreateVarBasic(scip_b, &var_b, name, 0.0, (double) string_size, 1.0, SCIP_VARTYPE_INTEGER));
// #endif
    
//     assert(var_b != NULL);

//     vars_b[nvars++] = var_b;

//     /* add variable to the problem */
//     SCIP_CALL(SCIPaddVar(scip_b, var_b));

//     /* add variable to corresponding constraint */
//     /* add variable to constraint (2) */
//     for (i = 0; i < sc_size; i++) {
//         for (j = 0; j < sum_sc[i + 1] - sum_sc[i]; j++) {
//             SCIP_CALL(SCIPaddCoefLinear(scip_b, conss_b[ROWS_2(i, j)], var_b, -1.0));
//         }
//     }

//     /**
//     * add variable df
//     */
//     (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "d_f");
//     SCIPdebugMessage("create variable %s\n", name);

//     /* create a basic variable object */
// #ifdef CSP
//     SCIP_CALL(SCIPcreateVarBasic(scip_b, &var_b, name, 0.0, 0.0, -1.0, SCIP_VARTYPE_INTEGER));
// #else
//     SCIP_CALL(SCIPcreateVarBasic(scip_b, &var_b, name, (double) 0.0, (double) string_size, -1.0, SCIP_VARTYPE_INTEGER));
// #endif
//     assert(var_b != NULL);

//     vars_b[nvars++] = var_b;

//     /* add variable to the problem */
//     SCIP_CALL(SCIPaddVar(scip_b, var_b));

//     /* add variable to corresponding constraint */
//     /* add variable to constraint (3) */
//     for (i = 0; i < sf_size; i++) {
//         for (j = 0; j < sum_sf[i + 1] - sum_sf[i]; j++) {
//             SCIP_CALL(SCIPaddCoefLinear(scip_b, conss_b[ROWS_3(i, j)], var_b, -1.0));
//         }
//     }

//     /* create problem data */
//     SCIP_CALL(probdataCreate(scip_b, &probdata, vars_b, conss_b, nvars, ncons, sc_size, sf_size, alpha_size, string_size, string_size, 0, sc, sf, sum_sc, sum_sf, alphabet, probname));

//     SCIP_CALL(SCIPwriteOrigProblem(scip_b, "dssp.lp", "lp", FALSE)); /* grava na saida padrao ou em file */

//     /* set user problem data */
//     SCIP_CALL(SCIPsetProbData(scip_b, probdata));

//     // /* free local buffer arrays */
//     // SCIPfreeBufferArray(scip_b, &conss_b);
//     // SCIPfreeBufferArray(scip_b, &vars_b);


//     return bcpa;


    
// }
double heur_BCPA(SCIP* scip_o, char *filename) {
    SCIP* scip;
    char name[SCIP_MAXSTRLEN];
    printf("BCPA\n");
    SCIP_CALL( SCIPcreate(&scip) );
    SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
    #ifdef ROUNDING
    SCIP_CALL( SCIPsetRealParam(scip, "limits/time", TIME_LIMIT) );
    #else
    SCIP_CALL( SCIPsetRealParam(scip, "limits/time", 900) );
    #endif
    SCIPinfoMessage(scip, NULL, "\n");
    SCIPinfoMessage(scip, NULL, "************************************************\n");
    SCIPinfoMessage(scip, NULL, "* Running BCPA *\n");
    SCIPinfoMessage(scip, NULL, "************************************************\n");
    SCIPinfoMessage(scip, NULL, "\n");
 
    SCIP_CALL( setupProblem(scip, scip_o) );
 
    SCIPinfoMessage(scip, NULL, "Original problem:\n");
    SCIP_CALL( SCIPprintOrigProblem(scip, NULL, "cip", FALSE) );
 
    SCIPinfoMessage(scip, NULL, "\n");
    SCIP_CALL( SCIPpresolve(scip) );
 
    /* SCIPinfoMessage(scip, NULL, "Reformulated problem:\n");
    SCIP_CALL( SCIPprintTransProblem(scip, NULL, "cip", FALSE) );
    */
 
    SCIPinfoMessage(scip, NULL, "\nSolving...\n");
    SCIP_CALL( SCIPsolve(scip) );
 
    if( SCIPgetNSols(scip) > 0 )
    {
       SCIPinfoMessage(scip, NULL, "\nSolution:\n");
       SCIP_CALL( SCIPprintSol(scip, SCIPgetBestSol(scip), NULL, FALSE) );

       return SCIPgetPrimalbound(scip);
    } else {
        printf("NO SOLUTIONS\n");
    }


 
    // SCIP_CALL( SCIPfree(&scip) );
 
    return SCIP_OKAY;
}
/** sets up problem */
 
 SCIP_RETCODE setupProblem(
    SCIP*                 scip                /**< SCIP data structure */,
    SCIP*                 scip_o                /**< SCIP data structure */
    )
 {
    printf("Setup Problem\n");
    char name[SCIP_MAXSTRLEN];
    int i, j, k, l;
    SCIP_CONS** conss;
    SCIP_VAR** vars, *var;



    int string_size, alpha_size, nvars, ncons;
    char *alphabet;
    double t, t_max, t_melhor;
    char **sc, **sf;
    int sc_size, sf_size = 0;
    int kc, kf = 0;
    int  *sum_sc, *sum_sf;
    int dc, df;

    SCIP_SOL *bestSolution;
    SCIP_PROBDATA *probdata;
    SCIP_VAR **vars_o, *var_o;
    SCIP_Real solval;

    FILE *file;
    clock_t antes, agora;

    antes = clock();
    
    assert(scip_o != NULL);

    bestSolution = SCIPgetBestSol(scip_o);
  
    probdata = SCIPgetProbData(scip_o);
    assert(probdata != NULL);

    // nvars = SCIPprobdataGetNVars(probdata);
    vars_o = SCIPprobdataGetVars(probdata);
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

    SCIP_CALL( SCIPcreateProbBasic(scip, "BCPA") );
    

    /* set objective sense */
    SCIP_CALL(SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE));

    /* tell SCIP that the objective will be always integral */
    SCIP_CALL(SCIPsetObjIntegral(scip));

    /* Number of constraints */
    // SCIP_CALL(SCIPallocBufferArray(scip, &conss, CONS_1 + CONS_2 + CONS_3 + CONS_4));
    /* Number of constraints */
    SCIP_CALL(SCIPallocBufferArray(scip, &conss, CONS_1 + CONS_2 + CONS_3 + CONS_4));

    /* Number of variables */
    // SCIP_CALL(SCIPallocBufferArray(scip, &vars, VARS_X + VARS_Y + VARS_D));
    //SCIP_CALL(SCIPallocBufferArray(scip, &vars, VARS_X + VARS_Y + VARS_D));
    SCIP_CALL(SCIPallocBufferArray(scip, &vars, VARS_X + VARS_Y + VARS_D ));

    /* (1) constraint - only one symbol is chosen */
    for (i = 0; i < CONS_1; i++) {
        SCIP_CALL(SCIPcreateConsBasicLinear (scip, &conss[i], "(1)", 0, NULL, NULL, (double) alpha_size - 1.0, (double) alpha_size - 1.0));
        SCIP_CALL(SCIPaddCons(scip, conss[i]));
    }

    /* (2) constraint - is sci's symbol matched? */
    ncons = i;
    for (; i < ncons + CONS_2; i++) {
        SCIP_CALL(SCIPcreateConsBasicLinear (scip, &conss[i], "(2)", 0, NULL, NULL, -SCIPinfinity(scip), (double) string_size)); // string_size faz o papel de M(no modelo)
        SCIP_CALL(SCIPaddCons(scip, conss[i]));
    }

    /* (3) constraint - is sfi's symbol matched? */
    ncons = i;
    for (; i < ncons + CONS_3; i++) {
        SCIP_CALL(SCIPcreateConsBasicLinear (scip, &conss[i], "(3)", 0, NULL, NULL,  0.0, SCIPinfinity(scip)));
        SCIP_CALL(SCIPaddCons(scip, conss[i]));
    }

    /* (4) cosntraint - disjunctive constraint */
    ncons = i;
    for (; i < ncons + CONS_4; i++) {
        SCIP_CALL(SCIPcreateConsBasicLinear (scip, &conss[i], "(4)", 0, NULL, NULL, 1.0, SCIPinfinity(scip)));
        SCIP_CALL(SCIPaddCons(scip, conss[i]));
    }
   
    nvars = 0;
    ncons = i;

     for(i = 0; i < string_size; i++) {
        for (j = 0; j < alpha_size; j++) {
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x_%d_%c", i, symbol_value(alphabet, alpha_size, j));
            SCIPdebugMessage("create variable %s\n", name);

            /* create a basic variable object */
           
            

            SCIP_Real solval = SCIPgetSolVal(scip_o, bestSolution, vars_o[nvars]); // guardar numa variavel
                    
            if (solval > 0.999999) { 
                SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, 1.0, 1.0, 0.0, SCIP_VARTYPE_BINARY));
            } else if (solval < 0.000001) {
                SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, 0.0, 0.0, 0.0, SCIP_VARTYPE_BINARY));
            } else {
                SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY));
            }

            assert(var != NULL);

            vars[nvars++] = var;

            /* add variable to the problem */
            SCIP_CALL(SCIPaddVar(scip, var));


            /* add variable to corresponding constraint */
            /* add variable to constraint (1) */
            SCIP_CALL(SCIPaddCoefLinear(scip, conss[ROWS_1(i)], var, 1.0));

            /* add variable to constraint (2) */
            for (k = 0; k < sc_size; k++) {
                for (l = 0; l < sum_sc[k + 1] - sum_sc[k]; l++) {
                    if (j == value_symbol(alphabet, alpha_size, sc[k][l + i]))
                        SCIP_CALL(SCIPaddCoefLinear(scip, conss[ROWS_2(k, l)], var, 1.0));
                }
            }

            /* add variable to constraint (3) */
            for (k = 0; k < sf_size; k++) {
                for (l = 0; l < sum_sf[k + 1] - sum_sf[k]; l++) {
                    if (j == value_symbol(alphabet, alpha_size, sf[k][l + i]))
                        SCIP_CALL(SCIPaddCoefLinear(scip, conss[ROWS_3(k, l)], var, 1.0));
                }
            }
            
            // if (arrayOfDefinedVars[i * alpha_size + j] != -1) {
            //     SCIP_CALL(SCIPaddCoefLinear(scip, conss[ROWS_5(bcpa)], var, 1.0));
            //     bcpa++;
            // }
        }
    }

     /**
    * add variables yij
    */
    for (i = 0; i < sc_size; i++) {
        for (j = 0; j < sum_sc[i + 1] - sum_sc[i]; j++) {
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "y_%d_%d", i, j);
            SCIPdebugMessage("create variable %s\n", name);
            
            
            
            /* create a basic variable object */
            SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY));
            assert(var != NULL);

            vars[nvars++] = var;

            /* add variable to the problem */
            SCIP_CALL(SCIPaddVar(scip, var));

            /* add variable to corresponding constraint */
            /* add variable to constraint (2) */
            SCIP_CALL(SCIPaddCoefLinear(scip, conss[ROWS_2(i, j)], var, string_size));

            /* add variable to constraint (4) */
            SCIP_CALL(SCIPaddCoefLinear(scip, conss[ROWS_4(i)], var, 1.0));
        }
    }

    /**
    * add variable dc
    */
    (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "d_c");
    SCIPdebugMessage("create variable %s\n", name);


    SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, 0.0, (double) kc, 1.0, SCIP_VARTYPE_INTEGER));
    
    assert(var != NULL);

    vars[nvars++] = var;

    /* add variable to the problem */
    SCIP_CALL(SCIPaddVar(scip, var));

    /* add variable to corresponding constraint */
    /* add variable to constraint (2) */
    for (i = 0; i < sc_size; i++) {
        for (j = 0; j < sum_sc[i + 1] - sum_sc[i]; j++) {
            SCIP_CALL(SCIPaddCoefLinear(scip, conss[ROWS_2(i, j)], var, -1.0));
        }
    }

    SCIPfreeBufferArray(scip, &conss);
    SCIPfreeBufferArray(scip, &vars);

    return SCIP_OKAY;
 }

double heur_RA(SCIP* scip, char *filename)
{
    int i, j, k, l, string_size, alpha_size, *heur;
    double ra;
    char name[SCIP_MAXSTRLEN], *alphabet;
    FILE *file;
    int         sc_size;
    int*        sum_sc;

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
    sc_size = SCIPprobdataGetSc_size(probdata) == 0 ? SCIPprobdataGetSf_size(probdata) : SCIPprobdataGetSc_size(probdata);
    sum_sc = SCIPprobdataGetSum_sc(probdata);

#ifndef NOVO    
    (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s-antigo", filename);
#else
    (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s-novo", filename);
#endif
#ifdef SUB
    strcat(name, "-SUB");
#else
    strcat(name, "-STRING");
#endif

#ifdef CSP
    strcat(name, "-CSP");
#elif FSP
    strcat(name, "-FSP");
#else
    strcat(name, "-DSP");
#endif

    strcat(name, "-RA-F.sol");
    file = fopen(name, "w");

    if(!file)
    {
        //printf("\nProblem to create solution file: %s", name);
        return SCIP_INVALID;
    }
    

    heur = (int *)malloc(sizeof(int) * string_size);
#ifdef NOVO
    double *guarda_maior; int change=0, ind;
    guarda_maior = malloc (string_size * sizeof (double));
    for(i=0;i<string_size;i++)guarda_maior[i]=-1.0;
    for(i=0;i<string_size;i++)heur[i]=0;

    int kcounter = 0, kpos = 0;
    //printf("REG2\n"); fflush(stdout);
#ifdef MERGE    
    /**
    Encontra o Melhor alinhamento em cada String em SC
    Para estes alinhamentos, encontra a variável de melhor valor em cada posição para computar 
    o alinhamento final
    */

    // printf("\nMERGE %d\n", sc_size); 

    double menor = 2;
    int melhor_posicao, jj;

    int *melhor_alinhamento = (int*)malloc(sizeof(int)*sc_size);
    double *melhor_alinhamento_valor = (double*)malloc(sizeof(double)*string_size);
    kcounter = 0;
    double menor_soma = 0, soma=0;;
    
    for (i = 0; i < sc_size; i ++) {
      menor_soma = 0;//2*string_size;
	//        printf("%I: %d\n", i);fflush(stdout);
        for (j = 0; j < string_size; j ++) {
            soma = 0;
            for (k = 0; k < sum_sc[i+1] - sum_sc[i]; k++) {
	      menor = 0; //2;
                for (l = 0; l < alpha_size; l++) {
                    var = vars[kcounter++];
		    if(j+k > (3/4.0)*string_size){
		      continue;
		    }
                    solval = SCIPgetSolVal(scip, bestSolution, var); // guardar numa variavel
                    
		    //                    printf("%lf\n", solval);fflush(stdout);

                    if (solval > menor + 0.00001) { //< menor) {
                        // printf("MAIOR %d %.2lf\n", l, solval); fflush(stdout);
                        melhor_posicao = j;
                        menor = solval;
                        ind = kcounter;
                    }
                }
                soma += menor;
            }
            if (soma > menor_soma + 0.00001) {
                menor_soma = soma; 
                melhor_alinhamento[i] = j;
		//printf("\nMelhor soma da string %d eh na posicao %d com valor %lf", i, melhor_alinhamento[i], menor_soma);
            }
        }
    }


    for (i = 0; i < sc_size; i ++) //printf("%d ", melhor_alinhamento[i]);

    for (i = 0; i < string_size; i ++ ) melhor_alinhamento_valor[i] = 0;//2;

    kcounter=0;
    for (i = 0; i < sc_size; i ++) {
        for (j = 0; j < string_size; j ++) {
            for (k = 0; k < sum_sc[i+1] - sum_sc[i]; k++) {
	      menor = 0;//2;
                for (l = 0; l < alpha_size; l++) {
                    var = vars[kcounter++];
                    solval = SCIPgetSolVal(scip, bestSolution, var); // guardar numa variavel

                    if (solval > menor) {
                        // printf("MAIOR %d %.2lf\n", l, solval); fflush(stdout);
                        menor=solval;
                        ind = l;
                    }
                }

                if (j == melhor_alinhamento[i]) {
                    if (menor > melhor_alinhamento_valor[j]) {
                        heur[j] = l;
                        melhor_alinhamento_valor[j] = menor;
                    }
                }
            }
            
        }
    }
    // printf("\n"); 

    free(melhor_alinhamento);

#elif REDUCE
    /**
    Ordena as variáveis
    Seleciona a de melhor valor
    Torna 1 as variáveis que se tornam inviáveis quando a variável acima é escolhida
    */

    double menor = 2;
    int melhor_posicao, jj;

    int *melhor_alinhamento = (int*)malloc(sizeof(int)*sc_size);
    kcounter = 0;
    int menor_soma = 0;

    double *solutions, *positions;
     
    for (count = 0; count < string_size; count +=1) {
        for (i = 0; i < sc_size; i ++) {
            for (j = 0; j < string_size; j ++) {
                for (k = 0; k < sum_sc[i+1] - sum_sc[i]; k++) {
                    menor = 2;
                    for (l = 0; l < alpha_size; l++) {
                        var = vars[kcounter++];
                        solval = SCIPgetSolVal(scip, bestSolution, var);

                        if (solval < menor) {
                            
                            melhor_posicao = j;
                            maior=solval;
                            ind = l;
                        }
                    }

                    
                }
            }
        }


    }
  
#else

    int melhor_posicao, jj;

    int *melhor_alinhamento_valor = (double*)malloc(sizeof(double)*string_size);

    for (i = 0; i < string_size; i ++) melhor_alinhamento_valor[i] = 2;
    kcounter = 0;
    int menor_soma = 0; 
    double menor;
    for (i = 0; i < sc_size; i ++) {
        for (j = 0; j < string_size; j ++) {
            for (k = 0; k < sum_sc[i+1] - sum_sc[i]; k++) {
                menor = 2;
                for (l = 0; l < alpha_size; l++) {
                    var = vars[kcounter++];
                    solval = SCIPgetSolVal(scip, bestSolution, var); // guardar numa variavel

                    if (solval < menor) {
                        // printf("MAIOR %d %.2lf\n", l, solval); fflush(stdout);
                        melhor_posicao = j;
                        menor=solval;
                        ind = l;
                    }
                }

                if (menor < melhor_alinhamento_valor[j]) {
                    heur[j] = l;
                    melhor_alinhamento_valor[j] = menor;
                }
            }
        }
    }
  
    
#endif
    
    
#else
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
#endif

    // for (i = 0; i < string_size; i ++) {
    //     printf("%d ", heur[i]); fflush(stdout);
    // } printf("\n");

    ra = valorSolucao(scip, filename, heur);
    fprintf(file, "Objective Value(dc - df) = %d\n", ra);
    fprintf(file, "Palavra(x) = ");
    
    for(i = 0; i < string_size; i++) {
       fprintf(file, "%c", symbol_value(alphabet, alpha_size, heur[i]));
    }
    fprintf(file, "\n");        
    
    fprintf(file, "# ");        
    for (j = 0; j < alpha_size; j++) {
            fprintf(file, "%c ", symbol_value(alphabet, alpha_size, j));
        }
    fprintf(file, "\n");
    for(i = 0; i < string_size; i++) {
        fprintf(file, "%d ", i);
        for (j = 0; j < alpha_size; j++) {
            if (heur[i] != j) {
                fprintf(file, "0 ");
            } else {
                fprintf(file, "1 ");
            }
        }
        fprintf(file, "\n");        
    }
    //printf("#3\n"); fflush(stdout);
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
    int totalFrac, desl, sc_size, *v, kcounter, k, l;
    
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
    


#ifndef NOVO    
    (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s-ANTIGO", filename);
#else
    (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s-NOVO", filename);
#endif
#ifdef GLOBAL
    strcat(name, "-GLOBAL");
#else
    strcat(name, "-ONLYROOT");
#endif
#ifdef SUB
    strcat(name, "-SUB");
#else
    strcat(name, "-STRING");
#endif

#ifdef CSP
    strcat(name, "-CSP");
#elif FSP
    strcat(name, "-FSP");
#else
    strcat(name, "-DSP");
#endif

#ifdef NONE
    strcat(name, "-NONE");
#elif RA
    strcat(name, "-RA");
#elif BCPA
    strcat(name, "-BCPA");
#elif VNS
    strcat(name, "-VNS");
#endif

#ifdef VBPL
    strcat(name, "-VBPL");
#elif BLPL
    strcat(name, "-BLPL");
#else
    strcat(name, "-PBPL");
#endif

#ifdef TIMESSS
    strcat(name, "-5S");
#elif TIMESS
    strcat(name, "-30S");
#elif TIMES
    strcat(name, "-60S");
#elif TIMEM
    strcat(name, "-120S");
#elif TIMEL
    strcat(name, "-180S");
#elif TIMELL
    strcat(name, "-300S");
#endif

#ifdef NSIZES
    strcat(name, "-SIZE2");
#elif NSIZEM
    strcat(name, "-SIZE3");
#elif NSIZEL
    strcat(name, "-SIZE4");
#else
    strcat(name, "-SIZE4");
#endif

#ifdef BEST
    strcat(name, "-BEST");
#elif FIRST
    strcat(name, "-FIRST");
#endif


#ifdef OMEGAS
    strcat(name, "-OMEGAS");
#elif OMEGAM
    strcat(name, "-OMEGAM");
#else
    strcat(name, "-OMEGAL");
#endif

    strcat(name, ".sol");
    file = fopen(name, "w");
    
    if(!file)
    {
        //printf("\nProblem to create solution file: %s", name);
        return;
    }
    
    name[0] = '\0';
    
    fprintf(file, "Tempo: %lfs\n", tempo);
        fprintf(file, "Objective Value(dc - df) = %lf\n", SCIPgetPrimalbound(scip));
    
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

#ifdef NOVO
    fprintf(file, "Var fracionarias = ");
    for (i = 0, kcounter = 0; i < sc_size; i ++) { //printf("F1\n"); fflush(stdout);
	 for (j = 0; j < string_size; j ++) { //printf("   F2\n"); fflush(stdout);
	    for (k = 0; k < sum_sc[i+1] - sum_sc[i]; k++) {  //printf("      F3\n"); fflush(stdout);
	      for (l = 0; l < alpha_size; l++) { //printf("         F4\n"); fflush(stdout);
                var = vars[kcounter++];
                solval = SCIPgetVarSol(scip, var); 
		if(!SCIPisFeasIntegral(scip, solval)){
		  fprintf(file, "\nx_%d_%d_%c_%d=%lf", i,j+k,symbol_value(alphabet, alpha_size, l),j, solval);
		}
	      }	
            }                    
	 }
    }

    fprintf(file, "Var inteiras 1.0 = ");
    for (i = 0, kcounter = 0; i < sc_size; i ++) { //printf("F1\n"); fflush(stdout);
	 for (j = 0; j < string_size; j ++) { //printf("   F2\n"); fflush(stdout);
	    for (k = 0; k < sum_sc[i+1] - sum_sc[i]; k++) {  //printf("      F3\n"); fflush(stdout);
	      for (l = 0; l < alpha_size; l++) { //printf("         F4\n"); fflush(stdout);
                var = vars[kcounter++];
                solval = SCIPgetVarSol(scip, var); 
		if(SCIPisFeasIntegral(scip, solval) && solval >= 1.0 - 0.00001){
		  fprintf(file, "\nx_%d_%d_%c_%d=%lf", i,j+k,symbol_value(alphabet, alpha_size, l),j, solval);
		}
	      }	
            }                    
	 }
    }    
    
#else  
    fprintf(file, "Var fracionarias = ");
    for(i = 0; i < string_size; i++) {
        for (j = 0; j < alpha_size; j++) {
            var = vars[i * alpha_size + j];
	    solval = SCIPgetVarSol(scip,var);
	    if(!SCIPisFeasIntegral(scip, solval))
	      {
		fprintf(file, "\nx_%d_%c=%lf", i,symbol_value(alphabet, alpha_size, j), solval);
	      }
        }
    }


    fprintf(file, "Var inteiras 0.0 = ");
    for(i = 0; i < string_size; i++) {
        for (j = 0; j < alpha_size; j++) {
            var = vars[i * alpha_size + j];
	    solval = SCIPgetVarSol(scip,var);
	    if(SCIPisFeasIntegral(scip, solval) && solval <= 0.00001)
	      {
		fprintf(file, "\nx_%d_%c=%lf", i,symbol_value(alphabet, alpha_size, j), solval);
	      }
        }
    }

    fprintf(file, "Var fracionarias = ");
    for(i = 0; i < sc_size; i++) {
        for (j = 0; j < sum_sc[i + 1] - sum_sc[i]; j++) {
            var = vars[VARS_X + sum_sc[i] + j];
	    solval = SCIPgetVarSol(scip,var);
	    if(!SCIPisFeasIntegral(scip, solval))
	      {
		fprintf(file, "\ny_%d_%d=%lf", i, j, solval);
	      }
        }
    }
    fprintf(file, "Var inteiras 1.0 = ");
    for(i = 0; i < sc_size; i++) {
        for (j = 0; j < sum_sc[i + 1] - sum_sc[i]; j++) {
            var = vars[VARS_X + sum_sc[i] + j];
	    solval = SCIPgetVarSol(scip,var);
	    if(SCIPisFeasIntegral(scip, solval) && solval >= 1.0 - 0.00001)
	      {
		fprintf(file, "\ny_%d_%d=%lf", i, j, solval);
	      }
        }
    }
#endif
    

    totalFrac=0;int totalNFRAC=0;
    for(i=0;i<nvars; i++){
      solval = SCIPgetVarSol(scip,vars[i]);
      if(!SCIPisFeasIntegral(scip, solval))
	{
	  totalFrac++;
	}else totalNFRAC +=1;
    }

    fprintf(file, "\nTotal fracionarias=%d\n", totalFrac);
    fclose(file);

    printf(";%d;%d;", totalFrac, totalNFRAC);
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


#ifndef NOVO    
    (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s-ANTIGO", filename);
#else
    (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s-NOVO", filename);
#endif
#ifdef GLOBAL
    strcat(name, "-GLOBAL");
#else
    strcat(name, "-ONLYROOT");
#endif
#ifdef SUB
    strcat(name, "-SUB");
#else
    strcat(name, "-STRING");
#endif

#ifdef CSP
    strcat(name, "-CSP");
#elif FSP
    strcat(name, "-FSP");
#else
    strcat(name, "-DSP");
#endif

#ifdef NONE
    strcat(name, "-NONE");
#elif RA
    strcat(name, "-RA");
#elif BCPA
    strcat(name, "-BCPA");
#elif VNS
    strcat(name, "-VNS");
#endif

#ifdef VBPL
    strcat(name, "-VBPL");
#elif BLPL
    strcat(name, "-BLPL");
#else
    strcat(name, "-PBPL");
#endif

#ifdef TIMESSS
    strcat(name, "-5S");
#elif TIMESS
    strcat(name, "-30S");
#elif TIMES
    strcat(name, "-60S");
#elif TIMEM
    strcat(name, "-120S");
#elif TIMEL
    strcat(name, "-180S");
#elif TIMELL
    strcat(name, "-300S");
#else
    strcat(name, "-3000S");
#endif

#ifdef NSIZES
    strcat(name, "-SIZE2");
#elif NSIZEM
    strcat(name, "-SIZE3");
#elif NSIZEL
    strcat(name, "-SIZE4");
#else
    strcat(name, "-SIZE4");
#endif

#ifdef BEST
    strcat(name, "-BEST");
#elif FIRST
    strcat(name, "-FIRST");
#endif


#ifdef OMEGAS
    strcat(name, "-OMEGAS");
#elif OMEGAM
    strcat(name, "-OMEGAM");
#else
    strcat(name, "-OMEGAL");
#endif
    strcat(name, ".sol");

    file = fopen(name, "w");
    
    if(!file)
    {
        printf("\nProblem to create solution file: %s\n", name);
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
#ifndef NOVO
    (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s-antigo.debug", filename);
#else
    (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s-antigo.debug", filename);
#endif
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

    if (sc_size > 0) {
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
    }
    else dc = 0;
    
    //Checar se obedece em sf
    df = string_size;
    if (sf_size > 0) {
        for (i = 0; i < sf_size; i ++) {
            maior_dif_sf = -1;
            for (j = 0; j <= (sum_sf[i+1]+string_size-1) - sum_sf[i] - string_size ; j ++) { l = 0;
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
    }
    else df = 0;

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
    if (sc_size > 0) {

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
    }
    else *dc = 0;
    
    //Checar se obedece em sf
    *df = string_size;

    if (sf_size > 0) {
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
    }
    else *df = 0;
    

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

void shaking (double **T, double **best_T, int tcount, int *x2, int k, double *ng, int string_size, int alpha_size) {
#ifdef DEBUG_VNS  
    // printf("Shaking\n");
#endif
    int i, j, r, m;
    int *viable_solutions = (int *) calloc(alpha_size, sizeof(int));
    int quant;
    double menor;
    int i_menor;

    //    srand(1);//time(NULL));
#ifdef VBPL
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
#elif BLPL

    for (i = 0; i < string_size; i ++) {
	    menor = 2.0; i_menor=2.0;
        for (j = 0; j < alpha_size; j ++) {
    
	        if(T[i][j]<menor - EPSILON){
	            menor = T[i][j];
	            i_menor = j;
	        }
        }
	    x2[i] = i_menor;
    }
    
    int number_positions = 0.10 * string_size;
    number_positions = randomInteger(0,string_size-1);
    if (number_positions == 0) number_positions = 1;
    
    
    for (i = 0; i < number_positions; i ++) {
        m = randomInteger(0,string_size-1);//rand();
        r = randomInteger(0,alpha_size-1);//rand();
        menor = 2.0;
        for (j = 0; j < alpha_size; j ++) {
            if (T[m][j] < menor) {
                menor = T[m][j];
                i_menor = j;
            }
        }   
        menor = T[m][i_menor];
        T[m][i_menor] = T[m][r];
        T[m][r] = menor;
        x2[m] = r;
    }
    
#else
    
    
#endif
#ifdef DEBUG_VNS		      
    for (i = 0; i < string_size; i ++) {
        //printf("%d ", x2[i]);
    }
    //printf("\n");
#endif
    free(viable_solutions);
    return;
}


void lpImprovement(int  *sum_sc, int *sum_sf, double **T, double** best_T, double **Y, int *x, int *x2, int k, int ky, double *ng, double *ny, char **sc, char**sf, int sc_size, int sf_size, int kc, int kf, char *alphabet, int alpha_size, int string_size, SCIP_PROBDATA *probdata,  SCIP_SOL *bestSolution ) {
    int neighbourhood = k;
    SCIP* scip;
    char name[SCIP_MAXSTRLEN];
    int i, j, l, m;
    SCIP_CONS** conss;
    SCIP_VAR** vars, *var;
    int fxaux;
    // printf("LPIMPROVEMNT\n");fflush(stdout);
 
 
    fxaux = objective_value(sum_sc, sum_sf, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size, x2);

    SCIP_CALL( SCIPcreate(&scip) );
    SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
    SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", 1) );
    SCIP_CALL( SCIPsetRealParam(scip, "limits/time", 30.00) );
    
    SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", 0) );
    //SCIP_CALL( SCIPsetRealParam(scip, "limits/memory", 1000) );
    
    // SCIPinfoMessage(scip, NULL, "\n");
    // SCIPinfoMessage(scip, NULL, "************************************************\n");
    // SCIPinfoMessage(scip, NULL, "* Running BCPA *\n");
    // SCIPinfoMessage(scip, NULL, "************************************************\n");
    // SCIPinfoMessage(scip, NULL, "\n");
 

    //printf("Setup Problem\n");



    int nvars, ncons;
    double t, t_max, t_melhor;

    SCIP_Real solval;

    FILE *file;
    clock_t antes, agora;

    antes = clock();
    
    assert(scip_o != NULL);

    
    SCIP_CALL( SCIPcreateProbBasic(scip, "VNS - LP") );
    

    /* set objective sense */
    SCIP_CALL(SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE));

    /* tell SCIP that the objective will be always integral */
    SCIP_CALL(SCIPsetObjIntegral(scip));

    /* Number of constraints */
    // SCIP_CALL(SCIPallocBufferArray(scip, &conss, CONS_1 + CONS_2 + CONS_3 + CONS_4));
    /* Number of constraints */
    SCIP_CALL(SCIPallocBufferArray(scip, &conss, CONS_1 + CONS_2 + CONS_3 + CONS_4));
    /* Number of variables */
    // SCIP_CALL(SCIPallocBufferArray(scip, &vars, VARS_X + VARS_Y + VARS_D));
    //SCIP_CALL(SCIPallocBufferArray(scip, &vars, VARS_X + VARS_Y + VARS_D));
    SCIP_CALL(SCIPallocBufferArray(scip, &vars, VARS_X + VARS_Y + VARS_D));

    /* (1) constraint - only one symbol is chosen */
    for (i = 0; i < CONS_1; i++) {
        SCIP_CALL(SCIPcreateConsBasicLinear (scip, &conss[i], "(1)", 0, NULL, NULL, (double) alpha_size - 1.0, (double) alpha_size - 1.0));
        SCIP_CALL(SCIPaddCons(scip, conss[i]));
    }

    /* (2) constraint - is sci's symbol matched? */
    ncons = i;
    for (; i < ncons + CONS_2; i++) {
        SCIP_CALL(SCIPcreateConsBasicLinear (scip, &conss[i], "(2)", 0, NULL, NULL, -SCIPinfinity(scip), (double) string_size)); // string_size faz o papel de M(no modelo)
        SCIP_CALL(SCIPaddCons(scip, conss[i]));
    }

    /* (3) constraint - is sfi's symbol matched? */
    ncons = i;
    for (; i < ncons + CONS_3; i++) {
        SCIP_CALL(SCIPcreateConsBasicLinear (scip, &conss[i], "(3)", 0, NULL, NULL,  0.0, SCIPinfinity(scip)));
        SCIP_CALL(SCIPaddCons(scip, conss[i]));
    }

    /* (4) cosntraint - disjunctive constraint */
    ncons = i;
    for (; i < ncons + CONS_4; i++) {
        SCIP_CALL(SCIPcreateConsBasicLinear (scip, &conss[i], "(4)", 0, NULL, NULL, 1.0, SCIPinfinity(scip)));
        SCIP_CALL(SCIPaddCons(scip, conss[i]));
    }
   
    nvars = 0;
    ncons = i;nvars=0;
    int posRandom;
    int *randomTeste = (int *)malloc(alpha_size * sizeof(int));
     for(i = 0; i < string_size; i++) {
        int soma_vars = 0;
        int pos_melhor=0; double melhor=2.0;
        for (k = 0; k < alpha_size; k++) {
            if (T[i][k] < (ng[neighbourhood] - 0.0001) && T[i][k] > (neighbourhood > 0 ? ng[neighbourhood-1] : -0.001) ) {
                soma_vars +=1;
            }
            if (T[i][k] < melhor) {
                pos_melhor = k;
                melhor = T[i][k];
            }
        }
        
        for (j = 0; j < alpha_size; j++) {
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x_%d_%c", i, symbol_value(alphabet, alpha_size, j));
            // SCIPdebugMessage("create variable %s\n", name);

            /* create a basic variable object */
            // int countRandom = 0;
            // if (j == 0) {
            //     ////printf("Comparando .....\n"); 
            //     for (l = 0; l < alpha_size; l ++) {
            //         //printf("T[i][l] %lf n %lf ", T[i][l], ng[k]);
            //         if (T[i][l] < ng[neighbourhood] && T[i][l] > (neighbourhood > 0 ? ng[neighbourhood-1] : -0.001) ) {
            //             randomTeste[countRandom] = l;
            //             //printf("YAY %d\n", countRandom);
            //             countRandom+=1;
            //         } else {
            //             //printf("NOP %d \n", countRandom);
            //         }
            //     }
            //     //printf("CountRandom %d\n", countRandom); 
            //     if (countRandom >= 0) {
            //         posRandom = randomInteger(0, countRandom-1);
            //     }
            //     else posRandom = -1;
            //     //printf("PosRandom definida como %d", posRandom);
            // }
            // //printf("J: %d  RandomTeste %d \n", j, randomTeste[posRandom]);
            // if (j == randomTeste[posRandom] && randomInteger(0,100) > 25)  { 
            //     SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, 0.0, 0.0, 0.0, SCIP_VARTYPE_BINARY));
            //     //printf("ZEROU %d %d ", i, j);
            // } else if (x2[i] == j){
            //     SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, 0.0, 0.5, 0.0, SCIP_VARTYPE_BINARY));
            // } else {
            //     SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY));
            //     //printf("Criou normal");
            // }
            if (soma_vars == 0) {
                    SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY));
            }else  if (T[i][j] < ng[neighbourhood] && T[i][j] > (neighbourhood > 0 ? ng[neighbourhood-1] : -0.001) ) {
                SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, 0.0, T[i][j], 0.0, SCIP_VARTYPE_BINARY));
            } else {
                SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, 1.0, 1.0, 0.0, SCIP_VARTYPE_BINARY));  
            }
            assert(var != NULL);

            vars[nvars++] = var;

            /* add variable to the problem */
            SCIP_CALL(SCIPaddVar(scip, var));


            /* add variable to corresponding constraint */
            /* add variable to constraint (1) */
            SCIP_CALL(SCIPaddCoefLinear(scip, conss[ROWS_1(i)], var, 1.0));

            /* add variable to constraint (2) */
            for (k = 0; k < sc_size; k++) {
                for (l = 0; l < sum_sc[k + 1] - sum_sc[k]; l++) {
                    if (j == value_symbol(alphabet, alpha_size, sc[k][l + i]))
                        SCIP_CALL(SCIPaddCoefLinear(scip, conss[ROWS_2(k, l)], var, 1.0));
                }
            }

            /* add variable to constraint (3) */
            for (k = 0; k < sf_size; k++) {
                for (l = 0; l < sum_sf[k + 1] - sum_sf[k]; l++) {
                    if (j == value_symbol(alphabet, alpha_size, sf[k][l + i]))
                        SCIP_CALL(SCIPaddCoefLinear(scip, conss[ROWS_3(k, l)], var, 1.0));
                }
            }

            
            // if (arrayOfDefinedVars[i * alpha_size + j] != -1) {
            //     SCIP_CALL(SCIPaddCoefLinear(scip, conss[ROWS_5(bcpa)], var, 1.0));
            //     bcpa++;
            // }
        }
    }

     /**
    * add variables yij
    */
     
    for (i = 0; i < sc_size; i++) {
        int soma_vars = 0;
        int pos_melhor = -1; double melhor = -1;
        
        for (k = 0; k < sum_sc[i + 1] - sum_sc[i]; k++) {
            if (Y[i][k] >= ny[ky] &&  Y[i][k] < ky == 0 ? 1.111 : ny[ky -1]) {
                soma_vars +=1;
            }
            if (Y[i][k] > melhor) {
                pos_melhor = k;
                melhor = Y[i][k];
            }
        }
        
        for (j = 0; j < sum_sc[i + 1] - sum_sc[i]; j++) {
            #ifdef BLPL
                (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "y_%d_%d", i, j);
                
                /* create a basic variable object */
                SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY));
                assert(var != NULL);
    
                vars[nvars++] = var;
    
                /* add variable to the problem */
                SCIP_CALL(SCIPaddVar(scip, var));
    
                /* add variable to corresponding constraint */
                /* add variable to constraint (2) */
                SCIP_CALL(SCIPaddCoefLinear(scip, conss[ROWS_2(i, j)], var, string_size));
    
                /* add variable to constraint (4) */
                SCIP_CALL(SCIPaddCoefLinear(scip, conss[ROWS_4(i)], var, 1.0));
            #else
                
                (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "y_%d_%d", i, j);
                // SCIPdebugMessage("create variable %s\n", name);
                
                /* create a basic variable object */
                // if (soma_vars == 0) {
                    SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY));
                // }else  if (Y[i][j] >= ny[ky] &&  Y[i][j] < ky == 0 ? 1.111 : ny[ky -1]) {
                    // SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, Y[i][j], 1.0, 0.0, SCIP_VARTYPE_BINARY));
                // } else {
                    // SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, 0.0, 0.0, 0.0, SCIP_VARTYPE_BINARY));    
                // }
                
                assert(var != NULL);
    
                vars[nvars++] = var;
    
                /* add variable to the problem */
                SCIP_CALL(SCIPaddVar(scip, var));
    
                /* add variable to corresponding constraint */
                /* add variable to constraint (2) */
                SCIP_CALL(SCIPaddCoefLinear(scip, conss[ROWS_2(i, j)], var, string_size));
    
                /* add variable to constraint (4) */
                SCIP_CALL(SCIPaddCoefLinear(scip, conss[ROWS_4(i)], var, 1.0));
            #endif
        }
    }

    /**
    * add variable dc
    */
    (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "d_c");
    // SCIPdebugMessage("create variable %s\n", name);


    SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, 0.0, (double) kc, 1.0, SCIP_VARTYPE_INTEGER));
    
    assert(var != NULL);

    vars[nvars++] = var;

    /* add variable to the problem */
    SCIP_CALL(SCIPaddVar(scip, var));

    /* add variable to corresponding constraint */
    /* add variable to constraint (2) */
    for (i = 0; i < sc_size; i++) {
        for (j = 0; j < sum_sc[i + 1] - sum_sc[i]; j++) {
            SCIP_CALL(SCIPaddCoefLinear(scip, conss[ROWS_2(i, j)], var, -1.0));
        }
    }


    // SCIPinfoMessage(scip, NULL, "Original problem:\n");
    // SCIP_CALL( SCIPprintOrigProblem(scip, NULL, "cip", FALSE) );
 
    // SCIPinfoMessage(scip, NULL, "\n");
    SCIP_CALL( SCIPpresolve(scip) );
 
    /* SCIPinfoMessage(scip, NULL, "Reformulated problem:\n");
    SCIP_CALL( SCIPprintTransProblem(scip, NULL, "cip", FALSE) );
    */
 
    // SCIPinfoMessage(scip, NULL, "\nSolving...\n");
    SCIP_CALL( SCIPsolve(scip) );

    
    bestSolution = SCIPgetBestSol(scip);
  
    probdata = SCIPgetProbData(scip);
   
    //    SCIPinfoMessage(scip, NULL, "\nSolution:\n");
    //    SCIP_CALL( SCIPprintSol(scip, SCIPgetBestSol(scip), NULL, FALSE) );
    //     

            for(i = 0; i < string_size; i++) {
            double menor = 2.0;
            int ind = -1;
            for (j = 0; j < alpha_size; j++) {
                var = vars[i * alpha_size + j];
                solval = SCIPgetSolVal(scip, bestSolution, var);
                // printf("%f ", solval);
                if (solval < menor) {
                    menor = solval;
                    ind = j;
                }
            }
            // printf("\n");
            x2[i] = ind;
            // printf("INDICE %d\n", ind);
        } 
        
        //printf("###################################\n");
        // for(i = 0; i < string_size; i++) {
            
        //     printf("%c ",symbol_value(alphabet, alpha_size,  x2[i]));
            
        // }printf("\n");
        
        // printf("%d\n", objective_value(sum_sc, sum_sf, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size, x2));
        if (objective_value(sum_sc, sum_sf, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size, x2) < objective_value(sum_sc, sum_sf, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size,x)){
            int nz = 0;
            for(i = 0; i < string_size; i++) {
                for (j = 0; j < alpha_size; j++) {
                    var = vars[nz++];
                    solval = SCIPgetSolVal(scip, bestSolution, var);
                    if (solval < 0.0000) solval = 0.0;
                    best_T[i][j] = solval;
                }
            }
            for (i = 0; i < sc_size; i ++) {
                for (j = 0; j < sum_sc[i + 1] - sum_sc[i]; j ++) {
                    var = vars[nz++];
                    solval = SCIPgetSolVal(scip, bestSolution, var);
                    if (solval <= -0.0001) solval = 0.0;
                    Y[i][j] = solval;
                }
            }
        }
       FILE *fp;
    //     #ifdef BLPL
    //   fp = fopen("relaxed1.txt", "a+");
    //      for(i = 0; i < string_size; i++) {
            
    //         fprintf(fp, "%c ",symbol_value(alphabet, alpha_size,  x2[i]));
            
    //     }fprintf(fp, " : ");
        
    //     fprintf(fp, "%d\n", objective_value(sum_sc, sum_sf, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size, x2));
    //   fclose (fp);
    //   #else
    //     fp = fopen("relaxed2.txt", "a+");
    //      fprintf(fp, "rand-20-5-5-250-1.dssp;%d;", k );
    //      for(i = 0; i < string_size; i++) {
            
    //         fprintf(fp, "%c",symbol_value(alphabet, alpha_size,  x2[i]));
            
    //     }fprintf(fp, ";");
        
    //     fprintf(fp, "%d;\n", objective_value(sum_sc, sum_sf, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size, x2));
    //   fclose (fp);
    //   #endif
    
    for(i = 0; i < ncons; i++) {
      SCIP_CALL( SCIPreleaseCons(scip, &(conss[i]) ));  
    }
    for(i = 0; i < nvars; i++ ){
      SCIP_CALL( SCIPreleaseVar(scip, &(vars[i]) ));
   }
    SCIPfreeBufferArray(scip, &conss);
    SCIPfreeBufferArray(scip, &vars);
    //free(scip);
    SCIP_CALL( SCIPfree(&scip) );
    free(randomTeste);
   
    // printf("Fim LP Improvement\n"); fflush(stdout);

    return;
}

void lpShaking(int  *sum_sc, int *sum_sf, double **T, double** best_T, double **Y, int *x, int *x2, int k, int ky, double *ng, double *ny, char **sc, char**sf, int sc_size, int sf_size, int kc, int kf, char *alphabet, int alpha_size, int string_size, SCIP_PROBDATA *probdata,  SCIP_SOL *bestSolution ) {
    int neighbourhood = k;
    SCIP* scip;
    char name[SCIP_MAXSTRLEN];
    int i, j, l, m;
    SCIP_CONS** conss;
    SCIP_VAR** vars, *var;
    int fxaux;
   
    fxaux = objective_value(sum_sc, sum_sf, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size, x2);
    // printf("lpShaking\n");fflush(stdout);

    SCIP_CALL( SCIPcreate(&scip) );
    SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
    SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", 1) );
    SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", 0) );
    SCIP_CALL( SCIPsetRealParam(scip, "limits/time", 60) );
     

    //printf("Setup Problem\n");

    int quantRandom = randomInteger(1, k);
    //    int quantRandom = randomInteger(1, string_size*0.5);
    int nvars, ncons;
    double t, t_max, t_melhor;

    SCIP_Real solval;

    FILE *file;
    clock_t antes, agora;

    antes = clock();
    
    assert(scip_o != NULL);

    
    SCIP_CALL( SCIPcreateProbBasic(scip, "VNS - LP") );
    
    /* set objective sense */
    SCIP_CALL(SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE));

    /* tell SCIP that the objective will be always integral */
    SCIP_CALL(SCIPsetObjIntegral(scip));

    /* Number of constraints */
    // SCIP_CALL(SCIPallocBufferArray(scip, &conss, CONS_1 + CONS_2 + CONS_3 + CONS_4));
    /* Number of constraints */
    SCIP_CALL(SCIPallocBufferArray(scip, &conss, CONS_1 + CONS_2 + CONS_3 + CONS_4 + quantRandom));
    /* Number of variables */
    // SCIP_CALL(SCIPallocBufferArray(scip, &vars, VARS_X + VARS_Y + VARS_D));
    //SCIP_CALL(SCIPallocBufferArray(scip, &vars, VARS_X + VARS_Y + VARS_D));
    SCIP_CALL(SCIPallocBufferArray(scip, &vars, VARS_X + VARS_Y + VARS_D));

    /* (1) constraint - only one symbol is chosen */
    for (i = 0; i < CONS_1; i++) {
        SCIP_CALL(SCIPcreateConsBasicLinear (scip, &conss[i], "(1)", 0, NULL, NULL, (double) alpha_size - 1.0, (double) alpha_size - 1.0));
        SCIP_CALL(SCIPaddCons(scip, conss[i]));
    }

    /* (2) constraint - is sci's symbol matched? */
    ncons = i;
    for (; i < ncons + CONS_2; i++) {
        SCIP_CALL(SCIPcreateConsBasicLinear (scip, &conss[i], "(2)", 0, NULL, NULL, -SCIPinfinity(scip), (double) string_size)); // string_size faz o papel de M(no modelo)
        SCIP_CALL(SCIPaddCons(scip, conss[i]));
    }

    /* (3) constraint - is sfi's symbol matched? */
    ncons = i;
    for (; i < ncons + CONS_3; i++) {
        SCIP_CALL(SCIPcreateConsBasicLinear (scip, &conss[i], "(3)", 0, NULL, NULL,  0.0, SCIPinfinity(scip)));
        SCIP_CALL(SCIPaddCons(scip, conss[i]));
    }

    /* (4) cosntraint - disjunctive constraint */
    ncons = i;
    for (; i < ncons + CONS_4; i++) {
        SCIP_CALL(SCIPcreateConsBasicLinear (scip, &conss[i], "(4)", 0, NULL, NULL, 1.0, SCIPinfinity(scip)));
        SCIP_CALL(SCIPaddCons(scip, conss[i]));
    }
    ncons = i;
    int kcons = i;

    /* (5) shaking constraints */
    for (; i < ncons + quantRandom; i++) {
      SCIP_CALL(SCIPcreateConsBasicLinear (scip, &conss[i], "(5)", 0, NULL, NULL, -SCIPinfinity(scip), 0.0));
        SCIP_CALL(SCIPaddCons(scip, conss[i]));
    }
    nvars = 0;
    ncons = i;nvars=0;
    int posRandom;
    int *randomTeste = (int *)malloc(alpha_size * sizeof(int));

    int* posicaoTrocar = (int *) malloc(quantRandom * sizeof(int));
    int* simboloTrocar = (int *) malloc(quantRandom * sizeof(int));
    int* posicoes = (int*) malloc(string_size * sizeof(int));
    int totalPosicoes = string_size, posicaoSorteada;

    for (i = 0; i < string_size; i++) {
      posicoes[i]=i;
    }

    for (i = 0; i < quantRandom; i ++) {
      // sortea uma posicao da string alvo a ser setada pelo shaking
      posicaoSorteada = randomInteger(0, totalPosicoes-1);
      posicaoTrocar[i] = posicoes[posicaoSorteada];
      posicoes[posicaoSorteada] = posicoes[--totalPosicoes];

      // sortea um simbolo para ser preferencial na posicao sorteada
      simboloTrocar[i] = randomInteger(0, alpha_size-1);
    }
    
    int countPosTrocar = 0;
     for(i = 0; i < string_size; i++) {
        int soma_vars = 0;
        int pos_melhor=0; double melhor=2.0;
        // for (k = 0; k < alpha_size; k++) {
        //     if (T[i][k] < (ng[neighbourhood] - 0.0001) && T[i][k] > (neighbourhood > 0 ? ng[neighbourhood-1] : -0.001) ) {
        //         soma_vars +=1;
        //     }
        //     if (T[i][k] < melhor) {
        //         pos_melhor = k;
        //         melhor = T[i][k];
        //     }
        // }
        int flag_trocar_pos = -1; // inicializa com -1
        for (k = 0; k < quantRandom; k ++) {
          if (posicaoTrocar[k] == i)
            flag_trocar_pos = k;    // seta com a posicao k
        }
        
        for (j = 0; j < alpha_size; j++) {
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x_%d_%c", i, symbol_value(alphabet, alpha_size, j));
            // SCIPdebugMessage("create variable %s\n", name);

            /* create a basic variable object */
            // int countRandom = 0;
            // if (j == 0) {
            //     ////printf("Comparando .....\n"); 
            //     for (l = 0; l < alpha_size; l ++) {
            //         //printf("T[i][l] %lf n %lf ", T[i][l], ng[k]);
            //         if (T[i][l] < ng[neighbourhood] && T[i][l] > (neighbourhood > 0 ? ng[neighbourhood-1] : -0.001) ) {
            //             randomTeste[countRandom] = l;
            //             //printf("YAY %d\n", countRandom);
            //             countRandom+=1;
            //         } else {
            //             //printf("NOP %d \n", countRandom);
            //         }
            //     }
            //     //printf("CountRandom %d\n", countRandom); 
            //     if (countRandom >= 0) {
            //         posRandom = randomInteger(0, countRandom-1);
            //     }
            //     else posRandom = -1;
            //     //printf("PosRandom definida como %d", posRandom);
            // }
            // //printf("J: %d  RandomTeste %d \n", j, randomTeste[posRandom]);
            // if (j == randomTeste[posRandom] && randomInteger(0,100) > 25)  { 
            //     SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, 0.0, 0.0, 0.0, SCIP_VARTYPE_BINARY));
            //     //printf("ZEROU %d %d ", i, j);
            // } else if (x2[i] == j){
            //     SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, 0.0, 0.5, 0.0, SCIP_VARTYPE_BINARY));
            // } else {
            //     SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY));
            //     //printf("Criou normal");
            // }
            
            if (flag_trocar_pos >= 0) // se posicao i eh uma das posicoes a ser setada pelo shaking
                SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY));
            else
              SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, best_T[i][j], 1.0, 0.0, SCIP_VARTYPE_BINARY));// fixar limite superior como best_T também?
            assert(var != NULL);


            vars[nvars++] = var;

            /* add variable to the problem */
            SCIP_CALL(SCIPaddVar(scip, var));


            /* add variable to corresponding constraint */
            /* add variable to constraint (1) */
            SCIP_CALL(SCIPaddCoefLinear(scip, conss[ROWS_1(i)], var, 1.0));

            /* add variable to constraint (2) */
            for (k = 0; k < sc_size; k++) {
                for (l = 0; l < sum_sc[k + 1] - sum_sc[k]; l++) {
                    if (j == value_symbol(alphabet, alpha_size, sc[k][l + i]))
                        SCIP_CALL(SCIPaddCoefLinear(scip, conss[ROWS_2(k, l)], var, 1.0));
                }
            }

            /* add variable to constraint (3) */
            for (k = 0; k < sf_size; k++) {
                for (l = 0; l < sum_sf[k + 1] - sum_sf[k]; l++) {
                    if (j == value_symbol(alphabet, alpha_size, sf[k][l + i]))
                        SCIP_CALL(SCIPaddCoefLinear(scip, conss[ROWS_3(k, l)], var, 1.0));
                }
            }

            /*(5)*/
            // for (k = 0; k < quantRandom; k ++) {
            //     if (i == posicaoTrocar[k] && j == simboloTrocar[k]) {
            //         for (l = 0; l < alpha_size; l ++)
            //             SCIP_CALL(SCIPaddCoefLinear(scip, conss[kcons+(i*(alpha_size) + l)], var, 1.0));
            //     } else if (i == posicaoTrocar[k]) {
            //         SCIP_CALL(SCIPaddCoefLinear(scip, conss[kcons+(i*(alpha_size) + j)], var, -1.0));
            //     }
            // }

            if(flag_trocar_pos >= 0){ // se >=0 , isto eh, i eh uma posicao sorteada
              if(simboloTrocar[flag_trocar_pos]==j){ // se j eh o simbolo sorteado para a posicao i
                //
                SCIP_CALL(SCIPaddCoefLinear(scip, conss[kcons + flag_trocar_pos], var, (double) alpha_size-1.0));
              }
              else{
                SCIP_CALL(SCIPaddCoefLinear(scip, conss[kcons + flag_trocar_pos], var, -1.0));
              }
            }
        
            /*            if (countPosTrocar < quantRandom && i == posicaoTrocar[countPosTrocar] && j == simboloTrocar[countPosTrocar]) {
                SCIP_CALL(SCIPaddCoefLinear(scip, conss[kcons + countPosTrocar], var, 1.0));
                flag = 1;
            } else if (i == posicaoTrocar[countPosTrocar] && countPosTrocar < quantRandom) {
                SCIP_CALL(SCIPaddCoefLinear(scip, conss[kcons + countPosTrocar], var, -1.0));
                flag = 1;
                } */
            
            // if (arrayOfDefinedVars[i * alpha_size + j] != -1) {
            //     SCIP_CALL(SCIPaddCoefLinear(scip, conss[ROWS_5(bcpa)], var, 1.0));
            //     bcpa++;
            // }
        }
        //        if (flag == 1) countPosTrocar++;
    }

     /**
    * add variables yij
    */
     
     for (i = 0; i < sc_size; i++) { // para cada string s em Sc
        int soma_vars = 0;
        int pos_melhor = -1; double melhor = -1;
        
        for (k = 0; k < sum_sc[i + 1] - sum_sc[i]; k++) { // verifica se y_sk estah em Nk e atualiza o maior y_sk
            if (Y[i][k] >= ny[ky] &&  Y[i][k] < ky == 0 ? 1.111 : ny[ky -1]) {
                soma_vars +=1;
            }
            if (Y[i][k] > melhor) {
                pos_melhor = k;
                melhor = Y[i][k];
            }
        }
        
        for (j = 0; j < sum_sc[i + 1] - sum_sc[i]; j++) {
            #ifdef SBPL_SEMY
                (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "y_%d_%d", i, j);
                
                /* create a basic variable object */
                SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY));
                assert(var != NULL);
    
                vars[nvars++] = var;
    
                /* add variable to the problem */
                SCIP_CALL(SCIPaddVar(scip, var));
    
                /* add variable to corresponding constraint */
                /* add variable to constraint (2) */
                SCIP_CALL(SCIPaddCoefLinear(scip, conss[ROWS_2(i, j)], var, string_size));
    
                /* add variable to constraint (4) */
                SCIP_CALL(SCIPaddCoefLinear(scip, conss[ROWS_4(i)], var, 1.0));
            #else
                
                (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "y_%d_%d", i, j);
                // SCIPdebugMessage("create variable %s\n", name);
                
                /* create a basic variable object */
                if (soma_vars == 0 && Y[i][j] > (melhor-0.05)) { // se nao existe nenhum y_ik em Nk => fixa o limite inferior do melhor (maior) dos y_ik com o valor da solucao atual
                    SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, Y[i][j], 1.0, 0.0, SCIP_VARTYPE_BINARY));
                }else  if (Y[i][j] >= ny[ky] &&  Y[i][j] < ky == 0 ? 1.111 : ny[ky -1]) { // se y_ij esta na vizinha Nk => fixa o limite inferior
                    SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, Y[i][j], 1.0, 0.0, SCIP_VARTYPE_BINARY));
                } else { // os demais sao descartados, isto eh, fixados em 0
                    SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, 0.0, 0.0, 0.0, SCIP_VARTYPE_BINARY));    
                }
                
                assert(var != NULL);
    
                vars[nvars++] = var;
    
                /* add variable to the problem */
                SCIP_CALL(SCIPaddVar(scip, var));
    
                /* add variable to corresponding constraint */
                /* add variable to constraint (2) */
                SCIP_CALL(SCIPaddCoefLinear(scip, conss[ROWS_2(i, j)], var, string_size));
    
                /* add variable to constraint (4) */
                SCIP_CALL(SCIPaddCoefLinear(scip, conss[ROWS_4(i)], var, 1.0));
            #endif
        }
    }
    
    /**
    * add variable dc
    */
    (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "d_c");
    // SCIPdebugMessage("create variable %s\n", name);


    SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, 0.0, (double) kc, 1.0, SCIP_VARTYPE_INTEGER));
    
    assert(var != NULL);

    vars[nvars++] = var;

    /* add variable to the problem */
    SCIP_CALL(SCIPaddVar(scip, var));

    /* add variable to corresponding constraint */
    /* add variable to constraint (2) */
    for (i = 0; i < sc_size; i++) {
        for (j = 0; j < sum_sc[i + 1] - sum_sc[i]; j++) {
            SCIP_CALL(SCIPaddCoefLinear(scip, conss[ROWS_2(i, j)], var, -1.0));
        }
    }


    // SCIPinfoMessage(scip, NULL, "Original problem:\n");
    // SCIP_CALL( SCIPprintOrigProblem(scip, NULL, "cip", FALSE) );
 
    // SCIPinfoMessage(scip, NULL, "\n");
    SCIP_CALL( SCIPpresolve(scip) );
 
    /* SCIPinfoMessage(scip, NULL, "Reformulated problem:\n");
    SCIP_CALL( SCIPprintTransProblem(scip, NULL, "cip", FALSE) );
    */
 
    // SCIPinfoMessage(scip, NULL, "\nSolving...\n");
    SCIP_CALL( SCIPsolve(scip) );

    SCIP_CALL(SCIPwriteOrigProblem(scip, "dssp_lpShaking.lp", "lp", FALSE));
    
    bestSolution = SCIPgetBestSol(scip);
  
    probdata = SCIPgetProbData(scip);
   
    //    SCIPinfoMessage(scip, NULL, "\nSolution:\n");
    //    SCIP_CALL( SCIPprintSol(scip, SCIPgetBestSol(scip), NULL, FALSE) );
    //     

            for(i = 0; i < string_size; i++) {
            double menor = 2.0;
            int ind = -1;
            for (j = 0; j < alpha_size; j++) {
                var = vars[i * alpha_size + j];
                solval = SCIPgetSolVal(scip, bestSolution, var);
                // printf("%f ", solval);
                if (solval < menor) {
                    menor = solval;
                    ind = j;
                }
            }
            // printf("\n");
            x2[i] = ind;
            // printf("INDICE %d\n", ind);
        } 
        
        //printf("###################################\n");
        // for(i = 0; i < string_size; i++) {
            
        //     printf("%c ",symbol_value(alphabet, alpha_size,  x2[i]));
            
        // }printf("\n");
        
        // printf("%d\n", objective_value(sum_sc, sum_sf, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size, x2));
        if (objective_value(sum_sc, sum_sf, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size, x2) < objective_value(sum_sc, sum_sf, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size,x)){
            int nz = 0;
            for(i = 0; i < string_size; i++) {
                for (j = 0; j < alpha_size; j++) {
                    var = vars[nz++];
                    solval = SCIPgetSolVal(scip, bestSolution, var);
                    if (solval < 0.0000) solval = 0.0;
                    best_T[i][j] = solval;
                }
            }
            for (i = 0; i < sc_size; i ++) {
                for (j = 0; j < sum_sc[i + 1] - sum_sc[i]; j ++) {
                    var = vars[nz++];
                    solval = SCIPgetSolVal(scip, bestSolution, var);
                    if (solval <= -0.0001) solval = 0.0;
                    Y[i][j] = solval;
                }
            }
        }
       FILE *fp;
    //     #ifdef BLPL
    //   fp = fopen("relaxed1.txt", "a+");
    //      for(i = 0; i < string_size; i++) {
            
    //         fprintf(fp, "%c ",symbol_value(alphabet, alpha_size,  x2[i]));
            
    //     }fprintf(fp, " : ");
        
    //     fprintf(fp, "%d\n", objective_value(sum_sc, sum_sf, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size, x2));
    //   fclose (fp);
    //   #else
    //     fp = fopen("relaxed2.txt", "a+");
    //      fprintf(fp, "rand-20-5-5-250-1.dssp;%d;", k );
    //      for(i = 0; i < string_size; i++) {
            
    //         fprintf(fp, "%c",symbol_value(alphabet, alpha_size,  x2[i]));
            
    //     }fprintf(fp, ";");
        
    //     fprintf(fp, "%d;\n", objective_value(sum_sc, sum_sf, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size, x2));
    //   fclose (fp);
    //   #endif
    
    for(i = 0; i < ncons; i++) {
      SCIP_CALL( SCIPreleaseCons(scip, &(conss[i]) ));  
    }
    for(i = 0; i < nvars; i++ ){
      SCIP_CALL( SCIPreleaseVar(scip, &(vars[i]) ));
   }
    SCIPfreeBufferArray(scip, &conss);
    SCIPfreeBufferArray(scip, &vars);
    //    free(scip);
    SCIP_CALL( SCIPfree(&scip) );    
    free(randomTeste);
    free(posicaoTrocar);
    free(simboloTrocar);
   
    //  printf("Fim Shaking\n"); fflush(stdout);

    return;
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
            int old_ix = ix;
            for (j = ix; j < string_size; j ++) {
	            for (m = iz; m < alpha_size; m ++) { // nao esta olhando toda a vizinhanca...
                    if (x_i[j] == m) continue;

#ifdef VNSFIRST                    
                    if (T[j][m] < ng[k] && T[j][m] > (k > 0 ? ng[k-1] : 0.0)) {
#ifdef DEBUG_VNS		      
                        printf("\nTroca x[%d]=%d por %d", j,x_i[j], m);
#endif
                        x_i[j] = m;
                        flag = 1;
                        if(old_ix == ix) ix = j;
                        iz = m+1; // tenta mudar outro simbolo na posicao j
                        if(iz==alpha_size){ // se tentou todos, muda de posicao
                            iz = 0;
                            ix = j+1;
                        }
                        break;
                    }
#endif
#ifdef VNSSECOND
                    if (T[j][m] < ng[k] && T[j][m] > 0.0) {
#ifdef DEBUG_VNS		      
                        printf("\nTroca x[%d]=%d por %d", j,x_i[j], m);
#endif
                        x_i[j] = m;
                        flag = 1;
                        if(old_ix == ix)  ix = j;
                        iz = m+1; // tenta mudar outro simbolo na posicao j
                        if(iz==alpha_size){ // se tentou todos, muda de posicao
                            iz = 0;
                            ix = j+1;
                        }
                        break;
                    }
#endif
                }
                    if (flag == 1) break;
                }
            i++;
            if (flag == 0) break;
            fxi = objective_value(sum_sc, sum_sf, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size, x_i);
#ifdef DEBUG_VNS		          
            printf("\tfx%d: %d\n", i, fxi);
#endif
        } while (fxaux <= fxi && i < 10000);
        if (flag == 0 || i >= 100) break;

        teste++;
        for (j = 0; j < string_size; j++) x2[j] = x_i[j];
    } while(fxaux >= fxi); //objective_value(sum_sc, sum_sf, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size, x2));
#ifdef DEBUG_VNS		      
    printf("\nfinal first improvement f(x')=%d", objective_value(sum_sc, sum_sf, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size, x2));
#endif

    free(x_aux);
    free(x_i);
    return;
}

double vns_main (SCIP* scip, char *filename, double tempo) {
    int i, j, string_size, alpha_size, *heur;
    double vns;
    char name[SCIP_MAXSTRLEN], *alphabet;
    int k, k_max;
    double t, t_max, t_melhor;
    #ifdef NSIZEM
    double ng[2] = {0.333, 0,666 };
    double ny[2] = {0.666, 0,333 };
    int ky = 0, ky_max = 3;
    #elif NSIZES
    double ng[1] = { 0.49999999999};
    double ny[1] = { 0.49999999999 };
    int ky = 0, ky_max = 2;
    #elif NSIZEL
    double ng[3] = {0.25555555, 0.49999999999, 0.745555555};
    double ny[3] = {0.745555555, 0.49999999999, 0.25555555 };
    int ky = 0, ky_max = 4;
    #elif NSIZELL
    double ng[4] = {0.2, 0.4, 0.6, 0.8};
    double ny[4] = {0.8, 0.4, 0.6, 0.8};
    int ky = 0, ky_max = 4;
    #endif
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

    srand(1);//time(NULL));    
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
    double **T, **best_T, **Y;
    T = (double **) calloc(string_size, sizeof(double *));
    for (i = 0; i < string_size; i++) {
        T[i] = (double *) calloc(alpha_size, sizeof(double));
    }
    best_T = (double **) calloc(string_size, sizeof(double *));
    for (i = 0; i < string_size; i++) {
        best_T[i] = (double *) calloc(alpha_size, sizeof(double));
    }
    Y = (double **) calloc(sc_size, sizeof(double *));
    for (i = 0; i < sc_size; i++) {
        Y[i] = (double *) calloc(sum_sc[i + 1] - sum_sc[i], sizeof(double));
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
    int nz = 0;
    for(i = 0; i < string_size; i++) {
        double menor = 2.0;
        int ind = -1;
        for (j = 0; j < alpha_size; j++) {
            var = vars[nz++];
            solval = SCIPgetSolVal(scip, bestSolution, var);
            if (solval <= -0.0001) solval = 0.0;
            T[i][j] = solval;
            best_T[i][j] = solval;
        }
    }
    for (i = 0; i < sc_size; i ++) {
        for (j = 0; j < sum_sc[i + 1] - sum_sc[i]; j ++) {
            var = vars[nz++];
            solval = SCIPgetSolVal(scip, bestSolution, var);
            if (solval <= -0.0001) solval = 0.0;
            Y[i][j] = solval;
            // printf("%lf ", Y[i][j]);fflush(stdout);
        }
    }

    // printf("\nT\n");
    //  for(i = 0; i < string_size; i++) {
    //     for (j = 0; j < alpha_size; j++) {
    //         printf("%.3lf ", T[i][j]);
    //     }
    //     printf("\n");
    // }
    // printf("\nY\n");
    // for(i = 0; i < sc_size; i++) {
    //     for (j = 0; j < string_size; j++) {
    //         printf("%.3lf ", Y[i][j]);
    //     }
    //     printf("\n");
    // } 
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

    t = 0;t_max=180; t_melhor = 0;
    #ifdef TIMESSS
    t_max = 5;
    #elif TIMESS
    t_max = 30;
    #elif TIMES
    t_max = 60;
    #elif TIMEM
    t_max = 120;
    #elif TIMEL
    t_max = 180;
    #elif TIMELL
    t_max = 300;
    #endif
    // printf("VNS_MAIN\n");fflush(stdout);
    int icount = 0, tcount = 0;;
    do {
        k = 0; icount = 0;
        do {
            
            #ifdef VBPL
                firstImprovement(sum_sc, sum_sf, T, x2, k, ng, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size);            
            #else
                lpImprovement(sum_sc, sum_sf, T, best_T, Y, x, x2, k, ky, ng, ny, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size, probdata, bestSolution);            
            #endif
            k = neighbourhoodChange(sum_sc, sum_sf, x, x2, k, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size);
    	    if(k==0){
    	        icount ++;
    	      // achou solucao melhor
    	      t_melhor = t;	      
    	    }
#ifdef DEBUG_VNS		      
            // printf("\n%d\n", k);
#endif
        } while (k < k_max);
        if (icount == 0) tcount = 0;
        tcount += 1;
        #ifdef SBPL
        
        lpShaking(sum_sc, sum_sf, T, best_T, Y, x, x2, k, ky, ng, ny, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size, probdata, bestSolution);            
        
        #else
        shaking(T, best_T, tcount, x2, k, ng, string_size, alpha_size);
        #endif
        ky += 1; 
        if (ky >= ky_max) ky = 0;
        agora = clock();
        t = (agora - antes) / ((double) CLOCKS_PER_SEC);
        // printf("%lf\n", t);fflush(stdout);

        
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
    fprintf(file, "Tempo da melhor: %lfs\n", t_melhor);
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
    
    for (i = 0; i < string_size; i++) {
        free (T[i]);
    }
    free (T);    
    
    for (i = 0; i < sc_size; i++) {
        free (Y[i]);
    }
    free (Y); 
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
    int barra = 0 , string_size, alpha_size, sc_size;//, *sum_sc;;
    double tempo, LB, UB, LB_ROOT, ra = SCIP_INVALID;//, vns = SCIP_INVALID;
    clock_t antes, agora;
    //SCIP_VAR **vars, *var;
    SCIP* scip = NULL;
    SCIP_PROBDATA *probdata;
    int sf_size;
    FILE *file; 
    double vns, bcpa;
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
    SCIP_CALL( SCIPincludeReadercsp(scip) );


    /* include default SCIP plugins */
    SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
    
     /* include Xyz constraint handler */
     //SCIP_CALL( SCIPincludeConshdlrXyz(scip) );   


    /* for column generation instances, disable restarts */
    SCIP_CALL( SCIPsetIntParam(scip,"presolving/maxrestarts",0) );

    /* turn off all separation algorithms */
    SCIP_CALL( SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE) );

    /* disable heuristics */
    SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, TRUE) );

// #if defined(BCPA) //|| defined(VNS)
//     SCIP_CALL( SCIPincludeHeurmyheuristic(scip) );
// #endif

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

       /* recupera dados do problema*/
    probdata = SCIPgetProbData(scip);
    string_size = SCIPprobdataGetString_size(probdata);
    sf_size = SCIPprobdataGetSf_size(probdata);
    sc_size = SCIPprobdataGetSc_size(probdata);
    alpha_size = SCIPprobdataGetAlpha_size(probdata);

    printf(";%d;%d;%d;%d",alpha_size, sc_size,sf_size,string_size);

    char solFileName[512];
    strcpy(solFileName, argv[2]);
    /* grava melhor solucao */
    //printSol(scip, argv[2], tempo);
#ifdef DSP
    strcat(solFileName, "-DSP");
#endif
#ifdef CSP
    strcat(solFileName, "-CSP");
#endif
#ifdef FSP
    strcat(solFileName, "-FSP");
#endif
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

strcat(solFileName, ".resumo");


    file = fopen(solFileName, "w");
    fprintf(file, "filename;alpha_size;sc_size;sf_size;string_size;frac;inteiras;lb;ub;root_lb,nodes,time;status;heuristica;valor\n");
    //imprimir nome
    for (int i = 0; i < strlen(argv[2]); i++) {
        if (argv[2][i] == '/')
            barra = i + 1;
    }
    for (int i = barra; i < strlen(argv[2]); i++) {
        fprintf(file, "%c", argv[2][i]);
    }

fprintf(file, ";%d;%d;%d;%d",alpha_size, sc_size,sf_size,string_size);
double tempo2=0.0;
ra = -9999;
#ifdef RA
    #ifdef ONLYROOT
    antes = clock();
    ra = heur_RA(scip, argv[2]);

    agora = clock();
    tempo = (agora - antes) / ((double) CLOCKS_PER_SEC);
    #endif
#endif
    // int bcpa;

#ifdef BCPA
    #ifdef ONLYROOT
    antes = clock();
    bcpa = heur_BCPA(scip, argv[2]);
    char name[SCIP_MAXSTRLEN];
    agora = clock();
    tempo2 = (agora - antes) / ((double) CLOCKS_PER_SEC);
    FILE *file2;
    (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "BCPA.csv");
        
    file2 = fopen(name, "w");
    fprintf(file2, "%s;%lf;%d;\n", argv[2], tempo+tempo2, bcpa);
    fclose(file2);
    #endif
#endif

#ifdef VNS
    #ifdef ONLYROOT
    antes = clock();    
    vns = vns_main(scip, argv[2], tempo);

    agora = clock();
    tempo2 = (agora - antes) / ((double) CLOCKS_PER_SEC);
    char name[SCIP_MAXSTRLEN];

    FILE *file2;
    (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "VNS.csv");
        
    file2 = fopen(name, "a");
    fprintf(file2, "%s;%lf;%d;\n", argv[2], tempo+tempo2, bcpa);
    fclose(file2);
    #endif
#endif

#ifdef ONLYROOT
    printRelaxedSol(scip, argv[2], tempo);
#endif

    /* valida solucao */
#ifdef DEBUG
    validSol(scip, argv[2]);
#endif

    // int kc = SCIPprobdataGetKc(probdata);
    // int kf = SCIPprobdataGetKf(probdata);
    
    UB = SCIPgetPrimalbound(scip);
    LB = SCIPgetDualbound(scip);
    LB_ROOT = SCIPgetDualboundRoot(scip);
    
    SCIP_Longint nNodes = SCIPgetNNodes (scip);
    printf("%lf;%lf;%lf;%lld;%lf;%d", LB, UB, LB_ROOT, nNodes, tempo,SCIPgetStatus(scip));
    fprintf(file, "%lf;%lf;%lf;%lld;%lf;%d", LB, UB, LB_ROOT, nNodes, tempo,SCIPgetStatus(scip));
    
    //printf("%lf;%lf;%lf;%lf;%d\n", LB, UB, SCIPgetDualboundRoot(scip), tempo, SCIPgetStatus(scip));
#ifdef RA
#ifdef ONLYROOT
        printf(";%lf;RA\n", ra);
#else
         printf(";%lf;GLOBAL+RA\n", UB);
        fprintf(file, ";RA\n", ra);
#endif

#elif BCPA 
        
    printf(";%lf;%lf;%lf;%lf;%d;BCPA;%d\n", LB, UB, LB_ROOT, tempo,SCIPgetStatus(scip), bcpa);
       
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
    printf(";%lf;%lf;NONE;ROOT\n", tempo,LB_ROOT);
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
    printf(";%lf;%lf;VNS\n", tempo2, vns);    
    fprintf(file, ";%lf;%lf;%lf;%lf;%d;%lf;VNS\n", LB, UB, LB_ROOT, tempo,SCIPgetStatus(scip), vns);    
#else
    printf(";%lf;%lf;NONE;NONE\n", tempo, LB_ROOT);
    fprintf(file, ";%lf;%lf;%lf;%lf;%d;NONE;NONE\n", LB, UB, LB_ROOT, tempo,SCIPgetStatus(scip));
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
