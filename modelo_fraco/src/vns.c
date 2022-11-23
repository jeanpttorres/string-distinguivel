#include "vns.h"
#include "utils.h"
#include "parameters_dssp.h"

int neighbourhoodChange(int  *sum_sc, int *sum_sf, int *x, int *x2, int k, char **sc, char**sf, int sc_size, int sf_size, int kc, int kf, char *alphabet, int alpha_size, int string_size, int *fx, int *fx2 ) {
    int i;

    *fx = objective_value(sum_sc, sum_sf, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size,x);
    *fx2 = objective_value(sum_sc, sum_sf, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size, x2);
    //    if (objective_value(sum_sc, sum_sf, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size, x2) < objective_value(sum_sc, sum_sf, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size,x)) {
    if(*fx2 < *fx){
#ifdef DEBUG_VNS
      printf("\nSolucao melhor encontrada... troca de %d por %d\n", *fx, *fx2);
#endif
        for (i = 0; i < string_size; i ++) x[i] = x2[i];
        k = 0;
    }
    else {
        k = k + 1;
    }
    return k;
}

/*
 * sorteia L posicoes e L simbolos. O tamanho L tambem eh sorteado (estava sorteando ate um limite de string_size, mudei para 10% de string_size, se string_size <=100, cc, usa-se 5% de string_size).
 * para cada posicao i sorteada, sorteia um simbolo k. 
 * Usa-se T.x* para setar x2 e
 * troca-se x2_i para ser k. Essa operacao tambem muda T.x*
 * Nenhum efeito eh causado em best_T e ele tambem nao eh usado. Nem ng
 * IN: T
 * OUT: T, x2
 */
void shaking (double **T, double **best_T, int tcount, int *x2, int k, double *ng, int string_size, int alpha_size, char* alphabet) {
#ifdef DEBUG_VNS  
    // printf("Shaking\n");
#endif
    int i, j, r, m;
    int *viable_solutions = (int *) calloc(alpha_size, sizeof(int));
    int quant;
    double menor;
    int i_menor;

    //    srand(1);//time(NULL));
    if(param.heur==VBPL){
       //#ifdef VBPL
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
    }
    else if(param.heur==BLPL || param.heur==SBPL){
       //#elif BLPL
       
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
       
       int number_positions = string_size<100?0.10 * string_size: 0.05 * string_size;//0.10 * string_size;
       number_positions = randomInteger(0,number_positions);//string_size-1);
       if (number_positions == 0) number_positions = 1;
       
       //    printf("\nSelecionado: %d posicoes para shaking!", number_positions);
       
       for (i = 0; i < number_positions; i ++) {
          m = randomInteger(0,string_size-1);//rand();
          r = randomInteger(0,alpha_size-1);//rand();
          menor = 2.0;
          for (j = 0; j < alpha_size; j ++) {
             if (T[m][j] < menor - EPSILON) {
                menor = T[m][j];
                i_menor = j;
             }
          }
          //        printf("\nShaking T[%d][%d]=%c por T[%d][%d]=%c", m,i_menor,alphabet[i_menor],m,r,alphabet[r]);
          menor = T[m][i_menor];
          T[m][i_menor] = T[m][r];
          T[m][r] = menor;
          x2[m] = r;
       }
    }
    else{
       //#else
       //#endif
    }
#ifdef DEBUG_VNS2
    printf("\nApos shaking: x2 = ");
    for (i = 0; i < string_size; i ++) {
       printf("%d ", x2[i]);
    }
    printf("\n");
#endif
    free(viable_solutions);
    return;
}

/*
 * resolve um PL com xij livre, se Nik={} ou se j \in Nik na solucao dada por T.x*. Cc, fixa xij=0 (nao selecionado, na implementacao do modelo estah 1)
 * x2 := solucao do RA deste PL guiado por T.
 * Se f(x2) < f(x) entao best_T.x* := solucao deste PL.
 * IN: T (solucao do PL a ser melhorada pela vizinhanca), x (solucao de referencia para atualizar best_T ou nao)
 * OUT: x2 (solucao vizinhanca), best_T (atualizado, somente se x2 melhora x)
 * EDNA: faltava T mudar como best_T
 */

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
    // SCIPinfoMessage(scip, NULL, "* Running LPIMPROV *\n");
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
          //            if (T[i][k] < (ng[neighbourhood] - 0.0001) && T[i][k] > (neighbourhood > 0 ? ng[neighbourhood-1] : -0.001) ) {
            if (T[i][k] < (ng[neighbourhood] + EPSILON) && T[i][k] > (neighbourhood > 0 ? ng[neighbourhood-1] - EPSILON : -0.001) ) {
                soma_vars +=1;
            }
            if (T[i][k] < melhor - EPSILON) {
                pos_melhor = k;
                melhor = T[i][k];
            }
        }
        
        for (j = 0; j < alpha_size; j++) {
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x_%d_%c", i, symbol_value(alphabet, alpha_size, j));
            //            printf("\nCria var %s. soma_vars=%d", name, soma_vars);
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
              SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY));//CONTINUOUS));//BINARY));
            }else if (T[i][j] < (ng[neighbourhood] + EPSILON) && T[i][j] > (neighbourhood > 0 ? ng[neighbourhood-1] - EPSILON : -0.001) ){ //if (T[i][j] < ng[neighbourhood] && T[i][j] > (neighbourhood > 0 ? ng[neighbourhood-1] : -0.001) ) {
              SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY));//CONTINUOUS));//BINARY));
              //              SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, 0.0, T[i][j], 0.0, SCIP_VARTYPE_CONTINUOUS));//BINARY));
            } else {
              SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, 1.0, 1.0, 0.0, SCIP_VARTYPE_BINARY));//CONTINUOUS));//BINARY));  
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
           if (Y[i][k] >= ny[ky] &&  Y[i][k] < (ky == 0 ? 1.111 : ny[ky -1])) {
                soma_vars +=1;
            }
            if (Y[i][k] > melhor) {
                pos_melhor = k;
                melhor = Y[i][k];
            }
        }
        
        for (j = 0; j < sum_sc[i + 1] - sum_sc[i]; j++) {
           if(param.heur==BLPL){
              //            #ifdef BLPL
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
           }
           else{
              //            #else
                
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
                //#endif
           }
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

    /**
    * add variable df
    */
    (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "d_f");
    // SCIPdebugMessage("create variable %s\n", name);


    SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, (double) kf, (double) string_size, -1.0, SCIP_VARTYPE_INTEGER));
    
    assert(var != NULL);

    vars[nvars++] = var;

    /* add variable to the problem */
    SCIP_CALL(SCIPaddVar(scip, var));

    /* add variable to corresponding constraint */
    /* add variable to constraint (3) */
    for (i = 0; i < sf_size; i++) {
        for (j = 0; j < sum_sc[i + 1] - sum_sc[i]; j++) {
            SCIP_CALL(SCIPaddCoefLinear(scip, conss[ROWS_3(i, j)], var, -1.0));
        }
    }


    // SCIPinfoMessage(scip, NULL, "Original problem:\n");
    // SCIP_CALL( SCIPprintOrigProblem(scip, NULL, "cip", FALSE) );

#ifdef DEBUG_VNS
   SCIP_CALL( SCIPwriteOrigProblem(scip, "vns-lpimprove.lp", "lp", FALSE) );
#endif
   // SCIPinfoMessage(scip, NULL, "\n");
    SCIP_CALL( SCIPpresolve(scip) );
 
    /* SCIPinfoMessage(scip, NULL, "Reformulated problem:\n");
    SCIP_CALL( SCIPprintTransProblem(scip, NULL, "cip", FALSE) );
    */

#ifdef DEBUG_VNS
    printf("\nGravado\n");
    //    getchar();
#endif
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
                if (solval < menor - EPSILON) {
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
            //            printf("\nLpimprovement achou melhor...");
            for(i = 0; i < string_size; i++) {
                for (j = 0; j < alpha_size; j++) {
                    var = vars[nz++];
                    solval = SCIPgetSolVal(scip, bestSolution, var);
                    if (solval < EPSILON) solval = 0.0;
                    best_T[i][j] = solval;
                    //                    T[i][j] = solval; // @EDNA: faltava isso!
                }
            }
            for (i = 0; i < sc_size; i ++) {
                for (j = 0; j < sum_sc[i + 1] - sum_sc[i]; j ++) {
                    var = vars[nz++];
                    solval = SCIPgetSolVal(scip, bestSolution, var);
                    if (solval < EPSILON) solval = 0.0;
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

/*
 * sorteia de 1 a k+1 posicoes i e simbolos j, R={(i,k)} e deixa o LP definir o valor das variaveis xij, para todo simbolo j, quando existe algum par (i,k), cc, deixa xij fixado com best_T.x*ij.
 * x2 := solucao perturbada 
 * T.x* := solucao do PL da solucao perturbada 
 * se f(x2) < f(x) entao atualiza best_T.
 * nao usados: ng
 * IN: best_T, x (dados da melhor solucao)
 * OUT: T, x2 (dados da solucao perturbada) + best_T (atualizado, se perturbado fica melhor)
 */
void lpShaking(int  *sum_sc, int *sum_sf, double **T, double** best_T, double **Y, int *x, int *x2, int K, int ky, double *ng, double *ny, char **sc, char**sf, int sc_size, int sf_size, int kc, int kf, char *alphabet, int alpha_size, int string_size, SCIP_PROBDATA *probdata,  SCIP_SOL *bestSolution ) {
    int neighbourhood = K;
    SCIP* scip;
    char name[SCIP_MAXSTRLEN];
    int i, j, l, m, k;
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

    int quantRandom = randomInteger(1, K+1); // (1,K) masKk pode ser 0
    //    int quantRandom = randomInteger(1, string_size*0.5);
    int nvars, ncons;
    double t, t_max, t_melhor;

    SCIP_Real solval;

    FILE *file;
    clock_t antes, agora;

    antes = clock();
    

    
    SCIP_CALL( SCIPcreateProbBasic(scip, "VNS - LP") );
    
    /* set objective sense */
    SCIP_CALL(SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE));

    /* tell SCIP that the objective will be always integral */
    SCIP_CALL(SCIPsetObjIntegral(scip));

    /* Number of constraints */
    // SCIP_CALL(SCIPallocBufferArray(scip, &conss, CONS_1 + CONS_2 + CONS_3 + CONS_4));
    /* Number of constraints */
    //    SCIP_CALL(SCIPallocBufferArray(scip, &conss, CONS_1 + CONS_2 + CONS_3 + CONS_4 + quantRandom));
    SCIP_CALL(SCIPallocBufferArray(scip, &conss, CONS_1 + CONS_2 + CONS_3 + CONS_4 + quantRandom*alpha_size));
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
    for (; i < ncons + quantRandom*alpha_size; i++) {
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
              SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS));//BINARY));
            else
              SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, best_T[i][j], 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS));//BINARY));// fixar limite superior como best_T também?
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
                //                  SCIP_CALL(SCIPaddCoefLinear(scip, conss[kcons + flag_trocar_pos], var, (double) alpha_size-1.0));
                for (int jj = 0; jj < alpha_size; jj++) {
                  if(jj==j)
                    SCIP_CALL(SCIPaddCoefLinear(scip, conss[kcons + flag_trocar_pos*alpha_size + jj], var, -1.0));
                  SCIP_CALL(SCIPaddCoefLinear(scip, conss[kcons + flag_trocar_pos*alpha_size + jj], var, 1.0));
                }
              }
              else{
                SCIP_CALL(SCIPaddCoefLinear(scip, conss[kcons + flag_trocar_pos*alpha_size + j], var, -1.0));
                //                SCIP_CALL(SCIPaddCoefLinear(scip, conss[kcons + flag_trocar_pos*quantRandom + j], var, (double) alpha_size-1.0));
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
           if (Y[i][k] >= ny[ky] &&  (Y[i][k] < (ky == 0 ? 1.111 : ny[ky -1]))) {
                soma_vars +=1;
            }
            if (Y[i][k] > melhor) {
                pos_melhor = k;
                melhor = Y[i][k];
            }
        }
        
        for (j = 0; j < sum_sc[i + 1] - sum_sc[i]; j++) {
           if(param.semy){
              //            #ifdef SBPL_SEMY
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
           }
           else{
           //            #else
                
                (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "y_%d_%d", i, j);
                // SCIPdebugMessage("create variable %s\n", name);
                
                /* create a basic variable object */
                if (soma_vars == 0 && Y[i][j] > (melhor-0.05)) { // se nao existe nenhum y_ik em Nk => fixa o limite inferior do melhor (maior) dos y_ik com o valor da solucao atual
                    SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, Y[i][j], 1.0, 0.0, SCIP_VARTYPE_BINARY));
                }else  if (Y[i][j] >= ny[ky] &&  Y[i][j] < (ky == 0 ? 1.111 : ny[ky -1])) { // se y_ij esta na vizinha Nk => fixa o limite inferior
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
                //#endif
           }
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

    /**
    * add variable df
    */
    (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "d_f");
    // SCIPdebugMessage("create variable %s\n", name);


    SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, (double) kf, (double) string_size, -1.0, SCIP_VARTYPE_INTEGER));
    
    assert(var != NULL);

    vars[nvars++] = var;

    /* add variable to the problem */
    SCIP_CALL(SCIPaddVar(scip, var));

    /* add variable to corresponding constraint */
    /* add variable to constraint (3) */
    for (i = 0; i < sf_size; i++) {
        for (j = 0; j < sum_sc[i + 1] - sum_sc[i]; j++) {
            SCIP_CALL(SCIPaddCoefLinear(scip, conss[ROWS_3(i, j)], var, -1.0));
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

#ifdef DEBUG_VNS
    SCIP_CALL(SCIPwriteOrigProblem(scip, "dssp_lpShaking.lp", "lp", FALSE));
#endif    
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
                if (solval < EPSILON) solval = 0.0;
                T[i][j] = solval;
                if (solval < menor - EPSILON) {
                    menor = solval;
                    ind = j;
                }
            }
            // printf("\n");
            x2[i] = ind;
            //            printf("INDICE %d\n", ind);
        } 
        
        //printf("###################################\n");
    //    for(i = 0; i < string_size; i++) {
            
    //      printf("%c ",symbol_value(alphabet, alpha_size,  x2[i]));
            
    //    }printf("\n");
        
    //            printf("%d\n", objective_value(sum_sc, sum_sf, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size, x2));
        if (objective_value(sum_sc, sum_sf, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size, x2) < objective_value(sum_sc, sum_sf, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size,x)){
            int nz = 0;
            for(i = 0; i < string_size; i++) {
                for (j = 0; j < alpha_size; j++) {
                    var = vars[nz++];
                    solval = SCIPgetSolVal(scip, bestSolution, var);
                    if (solval < EPSILON) solval = 0.0;
                    best_T[i][j] = solval;
                }
            }
            for (i = 0; i < sc_size; i ++) {
                for (j = 0; j < sum_sc[i + 1] - sum_sc[i]; j ++) {
                    var = vars[nz++];
                    solval = SCIPgetSolVal(scip, bestSolution, var);
                    if (solval <= EPSILON) solval = 0.0;
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
                    if(!param.vnssecond){
                       //#ifdef VNSFIRST                    
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
                       //#endif
                    }
                    else{
                       //#ifdef VNSSECOND
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
                       //#endif
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
        } while (fxaux <= fxi && i < 10000);
        if (flag == 0 || i >= 100) break;

        teste++;
        for (j = 0; j < string_size; j++) x2[j] = x_i[j];
    } while(fxaux >= fxi); //objective_value(sum_sc, sum_sf, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size, x2));
#ifdef DEBUG_VNS		      
    printf("\nfinal first improvement f(x')=%d\n", objective_value(sum_sc, sum_sf, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size, x2));
#endif

    free(x_aux);
    free(x_i);
    return;
}

void bestImprovement(int  *sum_sc, int *sum_sf, double **T, int *x2, int K, double *ng, char **sc, char**sf, int sc_size, int sf_size, int kc, int kf, char *alphabet, int alpha_size, int string_size, double time_max) {
    int *x_i, *x_aux,  i = 0, j = 0, m=0, ix = 0, iz=0, teste=0; 
    int fx2, fxaux, fxi, flag, best, mbest;
    double time = 0;
    
    x_aux = (int *) calloc(string_size, sizeof(int));
    x_i = (int *) calloc(string_size, sizeof(int));
#ifdef DEBUG_VNS		      
    printf("Best Improvement: \n");
#endif
    for (j = 0; j < string_size; j++)
      x_aux[j] = x2[j];
    for(i=0; i < string_size && time < time_max; i++) {
      // x_aux, f_xaux := melhor solucao atual
      fxaux = objective_value(sum_sc, sum_sf, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size, x_aux);
#ifdef DEBUG_VNS2		      
      printf("K: %d, f(x'): %d i=%d\n", K, fxaux, i);
#endif
      best = fxaux + 1;
      //      i = 0; ix = 0; iz = 0;
      //      do {
      for (j = 0; j < string_size; j++)
        x_i[j] = x_aux[j];//x2[j];
        
      //        flag = 0;
      int old;
      //      for (j = ix; j < string_size; j ++) {
      //      for (m = iz; m < alpha_size; m ++) { // nao esta olhando toda a vizinhanca...
      for (m = 0; m < alpha_size; m ++) { // para cada possivel vizinho
        if (x_i[i] == m) continue;
        old = x_i[i];
        flag = 0;
        //#ifdef VNSFIRST
        if(!param.vnssecond){
           if (T[i][m] < ng[K] +EPSILON && T[i][m] > (K > 0 ? ng[K-1] : 0.0) - EPSILON) {
#ifdef DEBUG_VNS2
              printf("\nAvalia vizinho - Troca x[%d]=%d por %d", i,x_i[i], m);
#endif
              x_i[i] = m;
              
              flag = 1;
              //              if(old_ix == ix) ix = j;
              //              iz = m+1; // tenta mudar outro simbolo na posicao j
              //              if(iz==alpha_size){ // se tentou todos, muda de posicao
              //                iz = 0;
              //                ix = j+1;
              //              }
              //              break;
           }
           //#endif
        }
        else{
           //#ifdef VNSSECOND
        if (T[i][m] < ng[K] && T[i][m] > 0.0) {
#ifdef DEBUG_VNS2
          printf("\nTroca x[%d]=%d por %d", i,x_i[i], m);
#endif
          x_i[i] = m;
          flag = 1;
              //              if(old_ix == ix)  ix = j;
              //              iz = m+1; // tenta mudar outro simbolo na posicao j
              //              if(iz==alpha_size){ // se tentou todos, muda de posicao
              //                iz = 0;
              //                ix = j+1;
              //              }
              //              break;
        }
        //#endif
        }
            //          }
            //          if (flag == 1) break;
            //        }
            //        i++;
        if (flag == 0) continue;
        fxi = objective_value(sum_sc, sum_sf, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size, x_i);
#ifdef DEBUG_VNS2		          
        printf("\tfx%d: %d\n", i, fxi);
#endif
        // guarda o melhor vizinho em Nik
        if(fxi < best){
          best = fxi;
          mbest = m;
        }
        // desfaz a mudanca feita pelo vizinho para tentar outro vizinho
        x_i[i] = old;
      } // for each Nik

      // atualiza se for melhor
      if(best < fxaux){
        //        printf("\nMuda para vizinho (troca x_aux[%d]=%c por %c)", i, alphabet[x_aux[i]], alphabet[mbest]);
        //        printf("\nMelhor Nik(i=%d,k=%d) = %d",i,K,best);
        fxaux = best;
        x_aux[i] = mbest;
        //        getchar();
      }
    } // muda para a proxima posicao da string alvo
    //while(fxaux >= best); //objective_value(sum_sc, sum_sf, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size, x2));
    for (j = 0; j < string_size; j++)
      x2[j] = x_aux[j];
#ifdef DEBUG_VNS		      
    printf("\nfinal best improvement f(x')=%d\n", objective_value(sum_sc, sum_sf, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size, x2));
#endif
    
    free(x_aux);
    free(x_i);
    return;
}

void imprimeSol(int* x, int* sum_sc, int* sum_sf, char** sc, char** sf, int sc_size, int sf_size, int kc, int kf, char* alphabet, int alpha_size, int string_size)
{
    double vns;
    int i, dc, df;
    printf("Solucao:\n");
    for (i = 0; i < string_size; i ++) 
        printf("%c" , alphabet[x[i]]);
    printf("\n");
    vns = dc_df_objective_value(sum_sc, sum_sf, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size, x, &dc, &df);
    printf("Valor: %lf\n", vns);
    printf("Dc: %d\n", dc);
    printf("Df: %d\n", df);
}
int vns_main (SCIP* scip, char *filename, double tempo, double* value, double* initial, int *x, int *pdc, int *pdf) {
  int i, j, string_size, alpha_size, *heur, fx, fx2;
  double vns, vns0;
    char name[SCIP_MAXSTRLEN], *alphabet;
    int k, k_max;
    double t, t_max, t_melhor;
    double ng[5], ny[5];
    int ky, ky_max;
    /*    #ifdef NSIZEM
    double ng[2] = {0.333, 0.666 }; // 0,666 };
    double ny[2] = {0.666, 0.333 }; // 0,333 };
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
    */
    char **sc, **sf;
    int sc_size, sf_size;
    int kc, kf;
    int  *sum_sc, *sum_sf;
    int dc, df;

    SCIP_SOL *bestSolution;
    SCIP_PROBDATA *probdata;
    SCIP_VAR **vars, *var;
    SCIP_Real solval;

    //    FILE *file;
    clock_t antes, agora;

    srand(1);//time(NULL));    
    antes = clock();
    
    assert(scip != NULL);

    bestSolution = SCIPgetBestSol(scip);
  
    probdata = SCIPgetProbData(scip);
    assert(probdata != NULL);


    k_max = param.ngsize;
    ky = 0;
    ky_max = k_max + 1;
    switch(param.ngsize){
    case 1:
       //    double ng[1] = { 0.49999999999};
       ng[0]=0.5;
       ny[0]=0.5;
       break;
    case 2:
       //    double ng[2] = {0.333, 0.666 }; // 0,666 };
       ng[0]=0.333; ng[1]=0.666;
       ny[0]=0.666; ny[1]=0.333;
       break;
    case 3:
       //    double ng[3] = {0.25555555, 0.49999999999, 0.745555555};
       ng[0]=0.25; ng[1]=0.5; ng[2]=0.75;
       ny[3]=0.25; ny[1]=0.5; ny[0]=0.75;
       break;
    case 4:
       //    double ng[4] = {0.2, 0.4, 0.6, 0.8};
       ng[0]=0.2; ng[1]=0.4; ng[2]=0.6; ng[3]=0.8;
       ny[3]=0.2; ny[2]=0.4; ny[1]=0.6; ny[0]=0.8;
       break;
    }
    
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
    for (i = 0; i < ky_max-1; i ++) {
      printf("%d:%lf ", i, ng[i]);
    }
#endif
    //    int *x; // best solution to be constructed
    //    x = (int *) calloc(string_size, sizeof(int));

    int *x2; // work solution that try to improve x
    x2 = (int *) calloc(string_size, sizeof(int));
    int nz = 0;
    for(i = 0; i < string_size; i++) {
        double menor = 2.0;
        int ind = -1;
        for (j = 0; j < alpha_size; j++) {
            var = vars[nz++];
            solval = SCIPgetSolVal(scip, bestSolution, var);
            if (solval < EPSILON) //<= -0.0001)
              solval = 0.0;
            T[i][j] = solval;
            best_T[i][j] = solval;
        }
    }
    for (i = 0; i < sc_size; i ++) {
        for (j = 0; j < sum_sc[i + 1] - sum_sc[i]; j ++) {
            var = vars[nz++];
            solval = SCIPgetSolVal(scip, bestSolution, var);
            if (solval < EPSILON) //<= -0.0001)
              solval = 0.0;
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
            if (solval < menor){// - EPSILON) {
                menor = solval;
                ind = j;
            }
        }
        x[i] = ind; // best solution
        x2[i] = x[i]; // first work solution is equal to x
    }

    
    vns0 = dc_df_objective_value(sum_sc, sum_sf, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size, x, &dc, &df);
    *initial = vns0;
#ifdef DEBUG_VNS2
    printf("Valor: %lf\n", vns0);
    printf("Dc: %d\n", dc);
    printf("Df: %d\n", df);
#endif
#ifdef DEBUG_VNS
    imprimeSol(x, sum_sc, sum_sf, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size);
#endif
    
#ifdef DEBUG_VNS2
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


    // (x, bestT) <- RA()
    // (x2, T) <- (x, bestT) // correto seria apenas T <- bestT
    // repeat
    //     k ← 1;
    //     repeat
    //         x2 ← LpImprovement(x, &bestT, &T) /* Local search (T). T muda somente se melhorar */;
    //         NeighbourhoodChange(&x, x2, &k) /* Change neighbourhood. So troca x, pois bestT e T ja alterados */;
    //     until k = k max ;
    //     /* shaking */;
    //     (x2,T)  ← LpShake(x, &bestT, k).
    //     t ← CpuTime()
    // until t > t max ;

    //    k_max = 3;

    t = 0;
    //    t_max=180;
    t_melhor = 0;
    /*    #ifdef TIMESSS
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
    */
    t_max = param.vnstime;
    
    // printf("VNS_MAIN\n");fflush(stdout);
    int icount = 0, tcount = 0;;
    do {
        k = 0; icount = 0;
        do {
            
           //            #ifdef VBPL
           switch (param.heur){
           case VBPL:
                firstImprovement(sum_sc, sum_sf, T, x2, k, ng, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size);
                break;
           case BLPL:           
           //            #else
           //               #ifdef BLPL
                bestImprovement(sum_sc, sum_sf, T, x2, k, ng, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size, 1800);
                break;
           case SBPL:
                //               #else
                //void lpImprovement(int  *sum_sc, int *sum_sf, double **T, double** best_T, double **Y, int *x, int *x2, int k, int ky, double *ng, double *ny, char **sc, char**sf, int sc_size, int sf_size, int kc, int kf, char *alphabet, int alpha_size, int string_size, SCIP_PROBDATA *probdata,  SCIP_SOL *bestSolution ) {
                lpImprovement(sum_sc, sum_sf, T, best_T, Y, x, x2, k, ky, ng, ny, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size, probdata, bestSolution);            
                //                lpImprovement(sum_sc, sum_sf, T, best_T, Y, x, x2, k, k, k_max, k_max, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size, probdata, bestSolution);
                //                #endif
                //            #endif
                break;
           }
           k = neighbourhoodChange(sum_sc, sum_sf, x, x2, k, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size, &fx, &fx2);
           if(k==0){
    	        icount ++;
    	      // achou solucao melhor
    	      t_melhor = t;
    	    }
#ifdef DEBUG_VNS
            printf("k=%d\nApos neighbourhoodChange x = ", k);
#endif
#ifdef DEBUG_VNS            
            imprimeSol(x, sum_sc, sum_sf, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size);
            if(fx2 < fx)
              getchar();
#endif
        } while (k < k_max);
        if (icount == 0) tcount = 0;
        tcount += 1;
        if(param.heur==SBPL){//param.shaking==LP_SHAKING){
           //        #ifdef SBPL
        
           lpShaking(sum_sc, sum_sf, T, best_T, Y, x, x2, k, ky, ng, ny, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size, probdata, bestSolution);            
        }
        else{
           //        #else
           shaking(T, best_T, tcount, x2, k, ng, string_size, alpha_size, alphabet);
           //        shaking(T, best_T, tcount, x2, k, k_max, string_size, alpha_size);
           //void shaking (double **T, double **best_T, int tcount, int *x2, int k, double *ng, int string_size, int alpha_size) {
        }
        //        #endif
        k += 1; 
        if (k >= k_max) k = 0;

#ifdef DEBUG_VNS
        printf("\nApos shaking. x2 = ");
        imprimeSol(x2, sum_sc, sum_sf, sc, sf, sc_size, sf_size, kc, kf, alphabet, alpha_size, string_size);
#endif

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

    /*
    if(!param.substring){
    //#ifndef SUB
       (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s-antigo-vns.sol", filename);
    }
    else{
    //#else
       (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s-antigo-vns-sub.sol", filename);
    //#endif
    }
    */
    
    /*
    (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s-vns.sol", filename);
    
    file = fopen(name, "w");

    if(!file)
    {
        printf("\nProblem to create solution file: %s", name);
        return 0;
    }
    name[0] = '\0';

    fprintf(file, "Time: %lfs\n", tempo);
    fprintf(file, "Best Solution Time: %lfs\n", t_melhor);
    fprintf(file, "Objective Value(dc - df) = %.0lf\n", vns);
    fprintf(file, "Target String(x) = ");

    for (i = 0; i < string_size; i ++) 
      fprintf(file,"%c" , alphabet[x[i]]);
    fprintf(file, "\n");

    / * dc * /
    fprintf(file, "dc = %d\n", dc);
    / * df * /
    fprintf(file, "df = %d\n", df);
    fclose(file);
    */
       
    for (i = 0; i < string_size; i++) {
       free (T[i]); free(best_T[i]);
    }
    free (T); free(best_T);
    
    for (i = 0; i < sc_size; i++) {
        free (Y[i]);
    }
    free (Y);
    //    free (x);
    free (x2);
    *value = vns;
    *pdc = dc;
    *pdf = df;
    return 1;
    // return 0.0;
}
