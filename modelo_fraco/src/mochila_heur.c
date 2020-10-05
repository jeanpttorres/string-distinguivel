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

#include "probdata_mochila.h"
#include "heur_myheuristic.h"

/* configuracao da heuristica */
#define HEUR_NAME             "myheuristic"
#define HEUR_DESC             "primal heuristic template"
#define HEUR_DISPCHAR         '!'
#define HEUR_PRIORITY         1 /* comeca pelas heuristicas de maior prioridade */
#define HEUR_FREQ             1 /* a cada 1 nivel da arvore de B&B */
#define HEUR_FREQOFS          0 /* comecando do nivel 0 */
#define HEUR_MAXDEPTH         -1 /* nivel max para chamar a heuristica. -1 = sem limites */
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE /* chamado depois que o LP resolvido */
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */

/* conversao da var x_ij no indice da variavel no vetor de variaveis do problema */
#define VAR_X(i,j) (i*m + j)
#ifdef DEBUG
   #define PRINTF(...) printf(__VA_ARGS__)
#else
   #define PRINTF(...) 
#endif

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
int randomInteger (int low, int high)
{
  int k;
  double d;

  d = (double) rand () / ((double) RAND_MAX + 1);
  k = d * (high - low + 1);
  return low + k;
}

int getLPSolution(SCIP* scip, double x[])
{
   SCIP_PROBDATA* probdata;
   int                   m;
   SCIP_VAR** vars;
   SCIP_VAR* var;
   double solval;
   int v;
   int nvars;

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);
  
   m=SCIPprobdataGetM(probdata);

   nvars = SCIPprobdataGetNVars(probdata);
   vars = SCIPprobdataGetVars(probdata);

   PRINTF("\nValue: %lf", SCIPgetLPObjval(scip));
   for( v=0; v< nvars; v++ )
     {
       var = vars[v];
       solval = SCIPgetVarSol(scip,var);
#ifdef DEBUG_HEUR       
       if( solval > EPSILON )
	 {
	   printf("\nx_%d_%d=%lf", v/m, v%m, solval);	   
	 }
#endif
       x[v]=solval;
     }
   return 1;
}


/*
 * Callback methods of primal heuristic
 */

/* TODO: Implement all necessary primal heuristic methods. The methods with an #if 0 ... #else #define ... are optional */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyMyheuristic)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeMyheuristic)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitMyheuristic)
{  /*lint --e{715}*/


   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitMyheuristic)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
SCIP_DECL_HEURINITSOL(heurInitsolMyheuristic)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolMyheuristic)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecMyheuristic)
{  /*lint --e{715}*/
   SCIP_PROBDATA* probdata;
   int                   n;
   int                   m;
   int*                  L;
   itemType*             item;
   SCIP_SOL*             sol;                /**< solution to round */
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandssol;
   int nlpcands;
   SCIP_VAR* var;
   int c;
   SCIP_Real oldsolval;
   SCIP_Real newsolval;
   SCIP_Bool stored;

   double* x;
   long *folga;
   int *coberto;
   int i, j;

  
   assert(result != NULL);
   assert(SCIPhasCurrentNodeLP(scip));

   *result = SCIP_DIDNOTRUN;

   /* continua somente se LP concluido */
   if ( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* continua somente se o valor do LP for menor que o limitante */
   if( SCIPisGE(scip, SCIPgetLPObjval(scip), SCIPgetCutoffbound(scip)) )
      return SCIP_OKAY;


   /* recupera as variaveis fracionarias que deveriam ser inteiras */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, NULL, &nlpcands, NULL, NULL) );

  
   /* se solucao do LP for inteira, pare   */
   if ( nlpcands == 0 )
     return SCIP_OKAY;

   /* cria uma estrutura para guardar uma  solucao */
   SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );

   /* copia a solucao do LP atual para a estrutura sol */
   SCIP_CALL( SCIPlinkLPSol(scip, sol) );

   /* recupera os dados da mochila */
   probdata=SCIPgetProbData(scip);
   assert(probdata != NULL);
  
   n=SCIPprobdataGetN(probdata);
   m=SCIPprobdataGetM(probdata);
   L=SCIPprobdataGetL(probdata);
   item=SCIPprobdataGetItem(probdata);

   x=(double*)malloc(sizeof(double)*n*m);
   folga=(long*)malloc(sizeof(long)*m);
   coberto=(int*)calloc(n, sizeof(int));
   getLPSolution(scip, x); // x nao esta sendo usado...

   for(i=0;i<m;i++){
     folga[i]=L[i];
     for(j=0;j<n;j++){       
       if(x[VAR_X(j,i)] < - EPSILON || x[VAR_X(j,i)] > 1 - EPSILON){
	   folga[i]-=item[j].peso;
	   coberto[j]=1;
	 }
     }
   }

#ifdef DEBUG_HEUR   
   for(i=0;i<m;i++){
     printf("\nFolga da mochila %d = %d", i, folga[i]);
   }

   for(i=0;i<n;i++){
     if(coberto[i]){
       printf("\nItem %d foi usado", i);
     }
   }
#endif
  
   for (; nlpcands>0; )
   {
     c = randomInteger(0,nlpcands-1);
     oldsolval = lpcandssol[c];
     assert( ! SCIPisFeasIntegral(scip, oldsolval) );
     var = lpcands[c];
     assert( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN ); /* var esta ativa no problema */

#ifdef ALEATORIO
     if (randomInteger(0,10)>9){
       newsolval = SCIPfeasCeil(scip, oldsolval);
     }
     else{
         newsolval = SCIPfeasFloor(scip, oldsolval);
     }
#else

     /* transforma o indice da variavel em indices i,j */
     i=SCIPvarGetProbindex(var)/m;
     j=SCIPvarGetProbindex(var)%m;

     /* se o item cabe na mochila j e item i ainda nao usado */
     if(item[i].peso<folga[j] - EPSILON && !coberto[i]){
       newsolval = SCIPfeasCeil(scip, oldsolval);
       folga[j]-=item[i].peso;
       coberto[i]=1;
     }
     else{
         newsolval = SCIPfeasFloor(scip, oldsolval);
     }
#endif
     PRINTF("my heuristic: var <%s>, x_%d_%d, val=%g newval=%g\n",  SCIPvarGetName(var), i, j, oldsolval, newsolval);
     SCIPdebugMessage("my heuristic: var <%s>, val=%g newval=%g\n",  SCIPvarGetName(var), oldsolval, newsolval);
     
     /* seta o novo valor da variavel na solucao */
     SCIP_CALL( SCIPsetSolVal(scip, sol, var, newsolval) );
     lpcands[c]=lpcands[--nlpcands];
   }

   /* verificar se a solucao eh viavel e armazena */
   SCIP_CALL( SCIPtrySol(scip, sol, FALSE, TRUE, FALSE, TRUE, &stored) );
   if( stored )
     {
#ifdef SCIP_DEBUG
       SCIPdebugMessage("found feasible rounded solution:\n");
       SCIP_CALL( SCIPprintSol(scip, sol, NULL, FALSE) );
#endif
       *result = SCIP_FOUNDSOL;
     }
   
   /* libera memoria alocada */
   free(x);
   free(folga);
   free(coberto);
   
   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the myheuristic primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurMyheuristic(
   SCIP*                 scip                /**< SCIP data structure */
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
         heurCopyMyheuristic, heurFreeMyheuristic, heurInitMyheuristic, heurExitMyheuristic, heurInitsolMyheuristic, heurExitsolMyheuristic, heurExecMyheuristic,
         heurdata) );
#else
   /* use SCIPincludeHeurBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecMyheuristic, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyMyheuristic) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeMyheuristic) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitMyheuristic) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitMyheuristic) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolMyheuristic) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolMyheuristic) );
#endif

   /* add myheuristic primal heuristic parameters */
   /* TODO: (optional) add primal heuristic specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
