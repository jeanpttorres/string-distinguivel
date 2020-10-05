/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_myheuristic.c
 * @brief  myheuristic primal heuristic
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <limits.h>

#include "heur_myheuristic.h"
#include "probdata_dssp.h"

#ifdef BCPA

#define HEUR_NAME             "myheuristic"
#define HEUR_DESC             "primal heuristic template"
#define HEUR_DISPCHAR         '?'
#define HEUR_PRIORITY         0
#define HEUR_FREQ             0     //Heuristica chamada apenas no nó raiz
#define HEUR_FREQOFS          0     //Heuristica chamada apenas no nó raiz
#define HEUR_MAXDEPTH         -1    //Sem limite
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */

#elif RA

#define HEUR_NAME             "myheuristic"
#define HEUR_DESC             "primal heuristic template"
#define HEUR_DISPCHAR         '?'
#define HEUR_PRIORITY         0
#define HEUR_FREQ             1     //Heuristica chamada apenas no nó raiz
#define HEUR_FREQOFS          0     //Heuristica chamada apenas no nó raiz
#define HEUR_MAXDEPTH         -1    //Sem limite
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */

#else

#define HEUR_NAME             "myheuristic"
#define HEUR_DESC             "primal heuristic template"
#define HEUR_DISPCHAR         '?'
#define HEUR_PRIORITY         0
#define HEUR_FREQ             0     //Heuristica chamada apenas no nó raiz
#define HEUR_FREQOFS          0    //Heuristica chamada apenas no nó raiz
#define HEUR_MAXDEPTH         -1    //Sem limite
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */

#endif


int flag = 0;
/*
 * Data structures
 */

/* TODO: fill in the necessary primal heuristic data */

/** primal heuristic data */
struct SCIP_HeurData
{
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */
static int objective_value() {
    return 0;
}

static int neighbourhood_change(int k) {
    return k+1;
}
static char * shaking(double **newX, int k, int alpha_size, char *alphabet, int string_size) {
    int i, j;
    char *xx;
    xx = (char *)malloc(sizeof(char) * string_size);
    for (i = 0; i < string_size; i ++) {
        for (j = 0; j < alpha_size; j ++) {
            printf("%d %c %lf\n", j, alphabet[j], newX[i][j]);
        }
    }
    return ("ASDJIASJD");
}
static char * first_improvement() {
    return "ASHDUASHD";
}
    






/*
 * Callback methods of primal heuristic
 */

/* TODO: Implement all necessary primal heuristic methods. The methods with an #if 0 ... #else #define ... are optional */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_HEURCOPY(heurCopymyheuristic)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of myheuristic primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurCopymyheuristic NULL
#endif

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_HEURFREE(heurFreemyheuristic)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of myheuristic primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurFreemyheuristic NULL
#endif


/** initialization method of primal heuristic (called after problem was transformed) */
#if 0
static
SCIP_DECL_HEURINIT(heurInitmyheuristic)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of myheuristic primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurInitmyheuristic NULL
#endif


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_HEUREXIT(heurExitmyheuristic)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of myheuristic primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurExitmyheuristic NULL
#endif


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_HEURINITSOL(heurInitsolmyheuristic)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of myheuristic primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurInitsolmyheuristic NULL
#endif


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_HEUREXITSOL(heurExitsolmyheuristic)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of myheuristic primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurExitsolmyheuristic NULL
#endif


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecmyheuristic)
{  /*lint --e{715}*/
   // Conta variaveis inteiras e fracionarias
   SCIP_VAR *var, *var_antiga;
   int **x, **y, sum, xInd, yInd, valid;
   int i, j, l, k, setvar = 0, ind, kc, kf;
   int dif, dc, menor_dif_sc, *hResult, df, maior_dif_sf,  sc_size, sf_size, string_size ; 
   int  *sum_sc, *sum_sf;
   SCIP_SOL*             sol;                /**< solution to round */
   SCIP_Real solval;    
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandssol;
   SCIP_Bool stored;
   int nlpcands; 
   int c;
   SCIP_Real oldsolval;
   SCIP_Real newsolval;
   SCIP_PROBDATA* probdata = SCIPgetProbData(scip);
   int alpha_size;
   char *alphabet, **sc, **sf;
    FILE *fp;

   int nvars = SCIPprobdataGetNVars(probdata);
   SCIP_VAR **vars = SCIPprobdataGetVars(probdata);
   SCIP_SOL *bestSolution;
   
   int totalFrac = 0;
   
   assert(scip != NULL);
   assert(probdata != NULL);

#ifdef BCPA
   int *fixed;
   int nFixed;

   fixed = (int*)malloc(sizeof(int)*nvars);
   assert(fixed!=NULL);
   for(i=0;i<nvars;i++){
     fixed[i]=0;
   }
   
#endif
   
   printf("\n***** HEURISTICA ***\n"); fflush(stdout);


   for(i = 0; i < nvars; i++) {
      solval = SCIPgetVarSol(scip, vars[i]);
      if(!SCIPisFeasIntegral(scip, solval)) {
         totalFrac++;
      }
#ifdef BCPA
      else {
        	//            SCIP_CALL(SCIPaddVarLocks(scip, vars[i], solval, solval));
        	//	SCIP_CALL(SCIPchgVarUbGlobal(scip, vars[i], 1.0));
        	//        SCIP_CALL(SCIPchgVarLbGlobal(scip, vars[i], 1.0));
        	if(solval>0.1){
        	  fixed[i]=2; // fixado em 1
        	}
        	else{
        	  fixed[i]=1; // fixado em 0
        	}
      }
#endif   
   }
   //
   if (flag == 0) {
        fp = fopen("bcpa-string.csv", "a");
        fprintf(fp, "%d;%d;", totalFrac, nvars - totalFrac);
        fclose(fp);       
   }
    
    flag = 1;
    

   /*
SCIP_RETCODE SCIPcopy 	( 	SCIP *  	sourcescip,
		SCIP *  	targetscip,
		SCIP_HASHMAP *  	varmap,
		SCIP_HASHMAP *  	consmap,
		const char *  	suffix,
		SCIP_Bool  	global,
		SCIP_Bool  	enablepricing,
		SCIP_Bool *  	valid 
	) 		

copies source SCIP to target SCIP; the copying process is done in the following order: 1) copy the plugins 2) create problem data in target-SCIP and copy the problem data of the source-SCIP 3) copy all active variables 4) copy all constraints 5) copy the settings

Note:
    all variables and constraints which are created in the target-SCIP are not (user) captured 

Parameters:
    sourcescip	source SCIP data structure
    targetscip	target SCIP data structure
    varmap	a hashmap to store the mapping of source variables corresponding target variables, or NULL
    consmap	a hashmap to store the mapping of source constraints to the corresponding target constraints, or NULL
    suffix	suffix which will be added to the names of the target SCIP, might be empty
    global	create a global or a local copy?
    enablepricing	should pricing be enabled in copied SCIP instance? If TRUE, pricer plugins will be copied and activated, and the modifiable flag of constraints will be respected. If FALSE, valid will be set to FALSE, when there are pricers present
    valid	pointer to store whether the copying was valid or not
   */


#ifdef BCPA

   SCIP* scip_copy;
   SCIP_Bool validBCPA;
   SCIP_RETCODE resultado;

   SCIP_CALL( SCIPcreate(&scip_copy) );
   
   resultado = SCIPcopyOrig(scip, scip_copy, NULL, NULL, "copy", FALSE, FALSE, &validBCPA );
   #ifdef ONLYROOT
    /* run only at root node */
    //SCIP_CALL( SCIPsetLongintParam(scip_copy, "limits/nodes", LONG_MAX) );
    #endif

   if(validBCPA){
     assert(scip_copy != NULL);
     /*     SCIP_PROBDATA *probdata_copy = SCIPgetProbData(scip_copy);
     assert(probdata_copy != NULL);
     
     int nvars_copy = SCIPprobdataGetNVars(probdata_copy);
     SCIP_VAR **vars_copy = SCIPprobdataGetVars(probdata_copy);
     */

     int nvars_copy = SCIPgetNVars(scip_copy);
     SCIP_VAR **vars_copy = SCIPgetVars(scip_copy);
     
     
     for(i=0;i<nvars_copy;i++){
       if(fixed[i]){
	 //printf("\nFixa var %d em %lf", i, fixed[i]-1.0);
	 SCIP_CALL(SCIPchgVarUbGlobal(scip_copy, vars_copy[i], fixed[i]-1.0));
	 SCIP_CALL(SCIPchgVarLbGlobal(scip_copy, vars_copy[i], fixed[i]-1.0));
       }
     }
     SCIP_CALL(SCIPwriteOrigProblem(scip_copy, "dssp2.lp", "lp", FALSE)); /* grava na saida padrao ou em file */

     SCIP_CALL( SCIPsolve(scip_copy) );

     SCIP_SOL* bestsol;
     //make READLINE=false ZIMPL=false LPS=cpx EXECUTE=ONLYROOT PROBLEM=STRING HEUR=NONE
     // make READLINE=false ZLIB=false GMP=false
     
     bestsol = SCIPgetBestSol(scip_copy);

    #ifdef STRING
    fp = fopen("bcpa-string.csv", "a");
    fprintf(fp, "BCPA;%lf;", SCIPgetDualbound(scip_copy));
    fprintf(fp, "BCPA2;%lf;", SCIPgetPrimalbound(scip_copy));
    
    fclose(fp);
    #endif

    #ifdef SUB
    fp = fopen("bcpa-sub.csv", "a");
    fprintf(fp, "BCPA;%lf;", SCIPgetDualbound(scip_copy));
    
    fclose(fp);
    #endif
    
    #ifdef GLOBAL

     if ( bestsol == NULL )
       SCIPinfoMessage(scip, NULL, "no solution available\n");
     else
       {
         SCIP_SOL* origsol;
	 
         SCIP_CALL( SCIPcreateSolCopy(scip_copy, &origsol, bestsol) );
         SCIP_CALL( SCIPretransformSol(scip_copy, origsol) );
         SCIP_CALL( SCIPprintSol(scip_copy, origsol, NULL, FALSE) );
         SCIP_CALL( SCIPfreeSol(scip_copy, &origsol) );
       }

	 /* ver o status da otimizacao da copia e se OK => fornecer essa solucao pro SCIP */
          

   SCIP_SOL*             sol;                //< solution to return
	// cria uma estrutura para guardar uma  solucao 
	SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );
   // seta o novo valor da variavel na solucao 
   

    // printf("################ PREPARATOION ##############\n");fflush(stdout);

    probdata=SCIPgetProbData(scip);
    assert(probdata != NULL);

    // printf("################ GET SOL ##############\n");fflush(stdout);


    // printf("################ PDATA ##############\n");fflush(stdout);

    SCIP_VAR **cp_vars = SCIPprobdataGetVars(probdata);
    string_size = SCIPprobdataGetString_size(probdata);
    alpha_size = SCIPprobdataGetAlpha_size(probdata);
    sc = SCIPprobdataGetSc(probdata);
    sf = SCIPprobdataGetSf(probdata);
    sc_size = SCIPprobdataGetSc_size(probdata);
    sf_size = SCIPprobdataGetSf_size(probdata);
    sum_sc = SCIPprobdataGetSum_sc(probdata);
    sum_sf = SCIPprobdataGetSum_sf(probdata);
    alphabet = SCIPprobdataGetAlphabet(probdata);

        bestSolution = SCIPgetBestSol(scip_copy);



    setvar = 0;
     for(i = 0; i < string_size; i++) {
        for (j = 0; j < alpha_size; j++) {
            
            var = vars[i*alpha_size + j];
            solval = SCIPgetSolVal(scip_copy, bestSolution, var);
            SCIP_CALL( SCIPsetSolVal(scip, sol, var, solval) );
            setvar++;
        }
    }

    for (i = 0; i < sc_size; i ++) {
        for (j = 0; j <= (sum_sc[i+1]+string_size-1) - sum_sc[i] -string_size ; j ++) {
            var = vars[setvar++ ];
            solval = SCIPgetSolVal(scip_copy, bestSolution, var);
            SCIP_CALL( SCIPsetSolVal(scip, sol, var, solval) ); 
            
        }
    }
    var = vars[setvar++ ];
    solval = SCIPgetSolVal(scip_copy, bestSolution, var);
    SCIP_CALL( SCIPsetSolVal(scip, sol, var, solval) ); 
    
    var = vars[setvar++ ];
    solval = SCIPgetSolVal(scip_copy, bestSolution, var);
    SCIP_CALL( SCIPsetSolVal(scip, sol, var, solval) ); 

    // verificar se a solucao eh viavel e armazena 
    SCIP_CALL( SCIPtrySol(scip_copy, sol, FALSE, TRUE, FALSE, TRUE, &stored) );
    if( stored )
    {
        printf("found feasible rounded solution:\n");fflush(stdout);
        SCIP_CALL( SCIPprintSol(scip_copy, sol, NULL, FALSE) );
        *result = SCIP_FOUNDSOL;
    }
    else {
        printf("Solução Não Adicionada\n");fflush(stdout);
    }
    
    #endif
    }
    SCIP_CALL( SCIPfree(&scip_copy) );

#endif   

#ifdef VNS
    printf("########################################");
    printf("########################################");
    // #ifdef GLOBAl
    /* continua somente se LP concluido */
    if ( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

    /* continua somente se o valor do LP for menor que o limitante */
    if( SCIPisGE(scip, SCIPgetLPObjval(scip), SCIPgetCutoffbound(scip)) )
      return SCIP_OKAY;


//    /* recupera as variaveis fracionarias que deveriam ser inteiras */
//    SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, NULL, &nlpcands, NULL, NULL) );
     
    /* cria uma estrutura para guardar uma  solucao */
    SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );  
   
    /* copia a solucao do LP atual para a estrutura sol */
    SCIP_CALL( SCIPlinkLPSol(scip, sol) );

//    /* recupera os dados do dssp */
    bestSolution = SCIPgetBestSol(scip);

    probdata = SCIPgetProbData(scip);
    assert(probdata != NULL);

    vars = SCIPprobdataGetVars(probdata);
    string_size = SCIPprobdataGetString_size(probdata);
    alpha_size = SCIPprobdataGetAlpha_size(probdata);
    sc = SCIPprobdataGetSc(probdata);
    sf = SCIPprobdataGetSf(probdata);
    sc_size = SCIPprobdataGetSc_size(probdata);
    sf_size = SCIPprobdataGetSf_size(probdata);
    sum_sc = SCIPprobdataGetSum_sc(probdata);
    sum_sf = SCIPprobdataGetSum_sf(probdata);
    alphabet = SCIPprobdataGetAlphabet(probdata);
    kc = SCIPprobdataGetKc(probdata);
    kf = SCIPprobdataGetKf(probdata);
    
    hResult = (int *)malloc(sizeof(int) * string_size);
    setvar = 0;
    
    double **newX;
    newX = (double **) calloc(string_size, sizeof(double *));
    for (i = 0; i < string_size; i++) {
        newX[i] = (double *) calloc(alpha_size, sizeof(double));
    }

    y = (int **) calloc(sc_size, sizeof(int *));
    for (i = 0; i < sc_size; i++) {
        y[i] = (int *) calloc(sum_sc[i], sizeof(int));
    }
    valid=1;
    
    for(i = 0; i < string_size; i++) {
        double menor = 2.0;
        int ind = -1;
        for (j = 0; j < alpha_size; j++) {
            
            
            var = vars[i*alpha_size + j];
            solval = SCIPgetSolVal(scip, bestSolution, var);
            SCIP_CALL( SCIPsetSolVal(scip, sol, var, 1.0) );
            
            newX[i][j] = solval;
        }
      
    }

    for (i = 0; i < string_size; i ++) {
        for (j = 0; j < alpha_size; j ++) {
            printf("%lf ", newX[i][j]);
        }
        printf("\n");
    }
    
    char *xx;
    xx = (char *)malloc(sizeof(char) * string_size);
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
    int tExec = 0, tTotal = 100, kMax = 3;
    
    do {
        k = 1;
        do {
            //x*  ← Shake(x)
            strcpy(xx, shaking(newX, k, alpha_size, alphabet, string_size));
            printf("Pos-shaking: %s\n", xx);
            //x ← FirstImprovement(x* ) /* Local search */;
            //strcpy(xx, first_improvement());
            
            //NeighbourhoodChange(x, x , k) /* Change neighbourhood */;
            k = neighbourhood_change(k);
            
        } while (k < kMax);
        tExec++;
    } while (tExec < tTotal);

    printf("########################################");
    printf("########################################");
#endif

#ifdef RA
    #ifdef GLOBAl
    /* continua somente se LP concluido */
    if ( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

    /* continua somente se o valor do LP for menor que o limitante */
    if( SCIPisGE(scip, SCIPgetLPObjval(scip), SCIPgetCutoffbound(scip)) )
      return SCIP_OKAY;


//    /* recupera as variaveis fracionarias que deveriam ser inteiras */
//    SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, NULL, &nlpcands, NULL, NULL) );
     
    /* cria uma estrutura para guardar uma  solucao */
    SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );  
   
    /* copia a solucao do LP atual para a estrutura sol */
    SCIP_CALL( SCIPlinkLPSol(scip, sol) );

//    /* recupera os dados do dssp */
    bestSolution = SCIPgetBestSol(scip);

    probdata = SCIPgetProbData(scip);
    assert(probdata != NULL);

    vars = SCIPprobdataGetVars(probdata);
    string_size = SCIPprobdataGetString_size(probdata);
    alpha_size = SCIPprobdataGetAlpha_size(probdata);
    sc = SCIPprobdataGetSc(probdata);
    sf = SCIPprobdataGetSf(probdata);
    sc_size = SCIPprobdataGetSc_size(probdata);
    sf_size = SCIPprobdataGetSf_size(probdata);
    sum_sc = SCIPprobdataGetSum_sc(probdata);
    sum_sf = SCIPprobdataGetSum_sf(probdata);
    alphabet = SCIPprobdataGetAlphabet(probdata);
    kc = SCIPprobdataGetKc(probdata);
    kf = SCIPprobdataGetKf(probdata);
    
    hResult = (int *)malloc(sizeof(int) * string_size);
    setvar = 0;
    
    
    x = (int **) calloc(string_size, sizeof(int *));
    for (i = 0; i < string_size; i++) {
        x[i] = (int *) calloc(alpha_size, sizeof(int));
    }

    y = (int **) calloc(sc_size, sizeof(int *));
    for (i = 0; i < sc_size; i++) {
        y[i] = (int *) calloc(sum_sc[i], sizeof(int));
    }
    valid=1;
    
    for(i = 0; i < string_size; i++) {
        double menor = 2.0;
        int ind = -1;
        for (j = 0; j < alpha_size; j++) {
            x[i][j] = 1;
            
            var = vars[i*alpha_size + j];
            solval = SCIPgetSolVal(scip, bestSolution, var);
            SCIP_CALL( SCIPsetSolVal(scip, sol, var, 1.0) );
            
            setvar++;
            if (solval < menor) {
                menor = solval;
                ind = i*alpha_size + j;//setvar - 1;
                xInd = j;
            }
        }
        SCIP_CALL( SCIPsetSolVal(scip, sol, vars[ind], 0.0) );
        hResult[i] = xInd;
        x[i][xInd] = 0;
    }
    
    //Checar dc
    dif = 0;
    dc = -string_size;
    for (i = 0; i < sc_size; i ++) {
        menor_dif_sc = (sum_sc[i+1]+string_size-1) - sum_sc[i];
        for (j = 0; j <= (sum_sc[i+1]+string_size-1) - sum_sc[i] -string_size ; j ++) {
            l = 0;
            dif = 0;
            
            var = vars[setvar++ ];
            SCIP_CALL( SCIPsetSolVal(scip, sol, var, 0.0) ); 
            y[i][j] = 0;
            for (k = j; k < j + string_size; k ++) {
                if (symbol_value(alphabet, alpha_size, hResult[l++]) != sc[i][k]) {
                  dif ++;
                }
            }
            if (dif < menor_dif_sc) {
                ind = setvar -1;
                yInd = j;
            }
            menor_dif_sc = menor_dif_sc < dif ? menor_dif_sc : dif;
            
        }
        var = vars[ind];
        SCIP_CALL( SCIPsetSolVal(scip, sol, var, 1.0) ); 
        y[i][yInd] = 1;
        dc = menor_dif_sc > dc ? menor_dif_sc : dc; 
    }
    var = vars[setvar++ ];
    SCIP_CALL( SCIPsetSolVal(scip, sol, var, dc) ); 
    

    //Checar df
    df = string_size;
    for (i = 0; i < sf_size; i ++) {
        maior_dif_sf = -1;
        for (j = 0; j <= (sum_sf[i+1]+string_size-1) - sum_sf[i] - string_size ; j ++) {
            l = 0;
            dif = 0;
            for (k = j; k < j + string_size; k ++) {
                if (symbol_value(alphabet, alpha_size, hResult[l++]) != sf[i][k]) {
                  dif ++;
                }
            }
            maior_dif_sf = dif > maior_dif_sf ? dif : maior_dif_sf;
        }
        df = maior_dif_sf < df ? maior_dif_sf : df; 
    }
    var = vars[setvar++ ];
    SCIP_CALL( SCIPsetSolVal(scip, sol, var, df) ); 
    
        #ifdef DEBUG
    //rest 1
    for(i = 0; i < string_size; i++) {
        sum = 0;
        for (j = 0; j < alpha_size; j++) {
            sum += x[i][j];
        }
        if (sum != alpha_size - 1) {
            printf("Problema na restricao (1) na posicao %d da string x.\n", i);
            valid = 0;
        }
    }
    //rest 2
    for (i = 0; i < sc_size; i++) {
        for (k = 0; k < sum_sc[i + 1] - sum_sc[i]; k++) {
            sum = 0;
            for (l = 0; l < string_size; l++) {
                ind = value_symbol(alphabet, alpha_size, sc[i][k + l]);
                if (ind != -1)
                    sum += x[l][ind];
            }
            if (sum > dc) {
                printf("Problema na restricao (2) na string %d posicao %d.\n", i, k);
                valid = 0;
            }
        }
    }
    //rest3
    for (i = 0; i < sf_size; i++) {
        for (k = 0; k < sum_sf[i + 1] - sum_sf[i]; k++) {
            sum = 0;
            for (l = 0; l < string_size; l++) {
                ind = value_symbol(alphabet, alpha_size, sf[i][k + l]);
                if (ind != -1)
                    sum += x[l][ind];
            }
            if (sum < df) {
                printf("Problema na restricao (3) na string %d posicao %d.\n", i, k);
                valid = 0;
            }
        }
    }
    
    //rest4
    for (i = 0; i < sc_size; i++) {
        sum = 0;
        for (k = 0; k < sum_sc[i + 1] - sum_sc[i]; k++) {
            sum += y[i][k];
        }
        if (sum < 1) {
            printf("Problema na restricao (4) na string %d.\n", i);
            valid = 0;
        }
    }
    
    if (dc > kc) {
        printf("Problema na restricao (5).\n");
        valid = 0;
    }
    
    if (df < kf) {
        printf("Problema na restricao (6).\n");
        valid = 0;
    }
    
    if (valid){
        printf("Solucao validada com sucesso!\n"); fflush(stdout);}
    else{
        printf("Falha na validacao da solucao!\n"); fflush(stdout);}

     
        #endif
        #ifdef DEBUG
    printf("RA dc:%d df: %d dc-df: %d\n", dc, df, dc - df);
        #endif
    /* verificar se a solucao eh viavel e armazena */
    SCIP_CALL( SCIPtrySol(scip, sol, FALSE, TRUE, FALSE, TRUE, &stored) );
    if(stored) {
        printf("found feasible rounded solution:\n"); fflush(stdout);
        SCIP_CALL( SCIPprintSol(scip, sol, NULL, FALSE) );

        *result = SCIP_FOUNDSOL;
    }
    else {
        printf("Solucao nao adicionada\n");fflush(stdout);
    }
    
    

//     for (; nlpcands>0; )
//    {
//     c = randomInteger(0,nlpcands-1);
//     oldsolval = lpcandssol[c];
//     assert( ! SCIPisFeasIntegral(scip, oldsolval) );
//     var = lpcands[c];
//     assert( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN ); /* var esta ativa no problema */

// #ifdef ALEATORIO
//     if (randomInteger(0,10)>9){
//        newsolval = SCIPfeasCeil(scip, oldsolval);
//     }
//     else{
//          newsolval = SCIPfeasFloor(scip, oldsolval);
//     }
// #else

//     /* transforma o indice da variavel em indices i,j */
//     i=SCIPvarGetProbindex(var)/m;
//     j=SCIPvarGetProbindex(var)%m;

//     /* se o item cabe na mochila j e item i ainda nao usado */
//     if(item[i].peso<folga[j] - EPSILON && !coberto[i]){
//        newsolval = SCIPfeasCeil(scip, oldsolval);
//        folga[j]-=item[i].peso;
//        coberto[i]=1;
//     }
//     else{
//          newsolval = SCIPfeasFloor(scip, oldsolval);
//     }
// #endif
//     PRINTF("my heuristic: var <%s>, x_%d_%d, val=%g newval=%g\n",  SCIPvarGetName(var), i, j, oldsolval, newsolval);
//     SCIPdebugMessage("my heuristic: var <%s>, val=%g newval=%g\n",  SCIPvarGetName(var), oldsolval, newsolval);
     
//     /* seta o novo valor da variavel na solucao */
//     SCIP_CALL( SCIPsetSolVal(scip, sol, var, newsolval) );
//     lpcands[c]=lpcands[--nlpcands];
//    }

//    /* verificar se a solucao eh viavel e armazena */
//    SCIP_CALL( SCIPtrySol(scip, sol, FALSE, TRUE, FALSE, TRUE, &stored) );
//    if( stored )
//     {
// #ifdef SCIP_DEBUG
//        SCIPdebugMessage("found feasible rounded solution:\n");
//        SCIP_CALL( SCIPprintSol(scip, sol, NULL, FALSE) );
// #endif
//        *result = SCIP_FOUNDSOL;
//}
    #endif
#endif

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the myheuristic primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurmyheuristic(
   SCIP*                scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create myheuristic primal heuristic data */
   heurdata = NULL;

   heur = NULL;

   /* include primal heuristic */
#if 0
   /* use SCIPincludeHeur() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP,
         heurCopymyheuristic, heurFreemyheuristic, heurInitmyheuristic, heurExitmyheuristic, heurInitsolmyheuristic, heurExitsolmyheuristic, heurExecmyheuristic,
         heurdata) );
#else
   /* use SCIPincludeHeurBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecmyheuristic, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopymyheuristic) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreemyheuristic) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitmyheuristic) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitmyheuristic) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolmyheuristic) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolmyheuristic) );
#endif

   /* add myheuristic primal heuristic parameters */
   /* TODO: (optional) add primal heuristic specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}