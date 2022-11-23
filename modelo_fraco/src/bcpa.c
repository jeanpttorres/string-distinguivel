#include "bcpa.h"
#include "parameters_dssp.h"

int heur_BCPA(SCIP* scip_o, char *filename, double *value, int *x, int *dc, int *df) {
    SCIP* scip;
    char name[SCIP_MAXSTRLEN];
    SCIP_SOL *bestSolution, **sols, *sol;
    SCIP_VAR** vars, *var;
    int string_size, alpha_size, sc_size, *sum_sc, i, j;
    double solval;
    
    printf("BCPA\n");
    SCIP_CALL( SCIPcreate(&scip) );
    SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
    if(param.bcpa_rounding){
       //    #ifdef ROUNDING
       SCIP_CALL( SCIPsetRealParam(scip, "limits/time", param.bcpatime_limit) );//TIME_LIMIT) );
    }
    else{
       //#else
       SCIP_CALL( SCIPsetRealParam(scip, "limits/time", param.bcpatime_limit) );//900) );
    }
    //    #endif
    SCIPinfoMessage(scip, NULL, "\n");
    SCIPinfoMessage(scip, NULL, "************************************************\n");
    SCIPinfoMessage(scip, NULL, "* Running BCPA *\n");
    SCIPinfoMessage(scip, NULL, "************************************************\n");
    SCIPinfoMessage(scip, NULL, "\n");
 
    SCIP_CALL( setupProblem(scip, scip_o, &vars) );
 
    SCIPinfoMessage(scip, NULL, "Original problem:\n");
    //    SCIP_CALL( SCIPprintOrigProblem(scip, NULL, "cip", FALSE) );
#ifdef DEBUG_VNS
     SCIP_CALL(SCIPwriteOrigProblem(scip, "dssp-bcpa.lp", "lp", FALSE)); /* grava na saida padrao ou em file */
#endif 
    SCIPinfoMessage(scip, NULL, "\n");
    //    SCIP_CALL( SCIPpresolve(scip) );
 
    /* SCIPinfoMessage(scip, NULL, "Reformulated problem:\n");
    SCIP_CALL( SCIPprintTransProblem(scip, NULL, "cip", FALSE) );
    */
 
    SCIPinfoMessage(scip, NULL, "\nSolving...\n");
    SCIP_CALL( SCIPsolve(scip) );
 
    if( SCIPgetNSols(scip) > 0 )
    {
       sols = SCIPgetSols(scip);
       
       SCIPinfoMessage(scip, NULL, "\nSolution:\n");
#ifdef DEBUG_VNS
       SCIP_CALL( SCIPprintSol(scip, SCIPgetBestSol(scip), NULL, FALSE) );
#endif
       //*value=SCIPgetPrimalbound(scip);
       //       bestSolution = SCIPgetBestSol(scip);

       *value = SCIPgetSolOrigObj(scip, sols[0]);
       bestSolution = sols[0];
       // falta salvar a solucao em filename
       string_size = SCIPprobdataGetString_size(SCIPgetProbData(scip_o));
       sc_size = SCIPprobdataGetSc_size(SCIPgetProbData(scip_o));
       //       alphabet = SCIPprobdataGetAlphabet(probdata);
       alpha_size = SCIPprobdataGetAlpha_size(SCIPgetProbData(scip_o));
       sum_sc = SCIPprobdataGetSum_sc(SCIPgetProbData(scip_o));

       var = vars[VARS_X + VARS_Y];
       *dc = (int) (SCIPgetSolVal(scip, bestSolution, var) + 0.5);
       
       var = vars[VARS_X + VARS_Y + 1];
       *df = (int) (SCIPgetSolVal(scip, bestSolution, var) + 0.5);
       //       printf("\nvalue=%lf dc (%s)=%d [%lf] df (%s)=%d [%lf]\n", *value, SCIPvarGetName(vars[VARS_X+VARS_Y]), *dc, SCIPgetSolVal(scip, bestSolution, vars[VARS_X+VARS_Y]), SCIPvarGetName(vars[VARS_X+VARS_Y+1]), *df, SCIPgetSolVal(scip, bestSolution, vars[VARS_X+VARS_Y+1]));
       //       SCIP_CALL( SCIPprintSol(scip, sols[0], NULL, FALSE) );
       for(i = 0; i < string_size; i++)
       {
          for (j = 0; j < alpha_size; j++)
          {
             var = vars[i * alpha_size + j];
             solval = SCIPgetSolVal(scip, bestSolution, var);
             if(solval < EPSILON){
                x[i]=j;
             }
          }
       }
       SCIPfreeBufferArray(scip, &vars);
       return 1;
    }
    else {
        printf("NO SOLUTIONS\n");
    }


    SCIPfreeBufferArray(scip, &vars);
    // SCIP_CALL( SCIPfree(&scip) );
 
    return 0;
}
/** sets up problem */
 
 SCIP_RETCODE setupProblem(
    SCIP*                 scip                /**< SCIP data structure */,
    SCIP*                 scip_o                /**< SCIP data structure */,
    SCIP_VAR*** pvars
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
           
            

            solval = SCIPgetSolVal(scip_o, bestSolution, vars_o[nvars]); // guardar numa variavel
                    
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
            SCIP_CALL(SCIPaddCoefLinear(scip, conss[ROWS_2(i, j)], var, (double) string_size));

            /* add variable to constraint (4) */
            SCIP_CALL(SCIPaddCoefLinear(scip, conss[ROWS_4(i)], var, 1.0));
        }
    }

    /**
    * add variable dc
    */
    (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "d_c");
    SCIPdebugMessage("create variable %s\n", name);


    /* create a basic variable object */
    if(param.problem==FSP){
       //#ifdef FSP
       SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, 0.0, 0.0, 1.0, SCIP_VARTYPE_INTEGER));
    }
    else{
       //#else
       SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, 0.0, (double) kc, 1.0, SCIP_VARTYPE_INTEGER));
       //#endif
    }

    //    SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, 0.0, (double) kc, 1.0, SCIP_VARTYPE_INTEGER));
    
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

     /**
     * add variable df
     */
     (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "d_f");
     SCIPdebugMessage("create variable %s\n", name);

     /* create a basic variable object */
     if(param.problem==CSP){
        //#ifdef CSP
        SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, 0.0, 0.0, -1.0, SCIP_VARTYPE_INTEGER));
     }
     else{
        //#else
        SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, (double) kf, (double) string_size, -1.0, SCIP_VARTYPE_INTEGER));
        //#endif
     }
    
     assert(var != NULL);

     vars[nvars++] = var;

     /* add variable to the problem */
     SCIP_CALL(SCIPaddVar(scip, var));

     /* add variable to corresponding constraint */
     /* add variable to constraint (3) */
     for (i = 0; i < sf_size; i++) {
        for (j = 0; j < sum_sf[i + 1] - sum_sf[i]; j++) {
           SCIP_CALL(SCIPaddCoefLinear(scip, conss[ROWS_3(i, j)], var, -1.0));
        }
     }


    SCIPfreeBufferArray(scip, &conss);
    //    SCIPfreeBufferArray(scip, &vars);

    *pvars = vars;
    return SCIP_OKAY;
 }
