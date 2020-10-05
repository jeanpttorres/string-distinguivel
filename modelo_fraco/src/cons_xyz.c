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

/**@file   cons_xyz.c
 * @brief  constraint handler for xyz constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "cons_xyz.h"
#include "probdata_mochila.h"


/* fundamental constraint handler properties */
#define CONSHDLR_NAME          "xyz"
#define CONSHDLR_DESC          "constraint handler template"
#define CONSHDLR_ENFOPRIORITY         -1 /**< priority of the constraint handler for constraint enforcing -- used after LP was solved . If you use negative priority, it implies that the variables in the solution are integers */
#define CONSHDLR_CHECKPRIORITY        -1 /**< priority of the constraint handler for checking feasibility -- used to check if a primal solution is feasible. idem */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                              *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

/* optional constraint handler properties */
/* TODO: remove properties which are never used because the corresponding routines are not supported */
#define CONSHDLR_SEPAPRIORITY         0 /**< priority of the constraint handler for separation */
#define CONSHDLR_SEPAFREQ            1 /**< frequency for separating cuts; zero means to separate only in the root node. -1 = never and 5 means every node in depth=multiple of 5, p.e. */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */

#define CONSHDLR_PROPFREQ            -1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_PROP_TIMING       SCIP_PROPTIMING_BEFORELP/**< propagation timing mask of the constraint handler*/

#define CONSHDLR_PRESOLTIMING    SCIP_PRESOLTIMING_MEDIUM /**< presolving timing of the constraint handler (fast, medium, or exhaustive) */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */




/* TODO: (optional) enable linear or nonlinear constraint upgrading */
#if 0
#include "scip/cons_linear.h"
#include "scip/cons_nonlinear.h"
#define LINCONSUPGD_PRIORITY          0 /**< priority of the constraint handler for upgrading of linear constraints */
#define NONLINCONSUPGD_PRIORITY       0 /**< priority of the constraint handler for upgrading of nonlinear constraints */
#endif


/*
 * Data structures
 */

/* TODO: fill in the necessary constraint data */

/** constraint data for xyz constraints */
struct SCIP_ConsData
{
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
};


/*
 * Local methods
 */



int deltaConstraint(SCIP* scip, SCIP_CONSHDLR* conshdlr, SCIP_SOL* sol, int tamanho, int* result)
{
   SCIP_PROBDATA* probdata;
   SCIP_Real solval;
   SCIP_VAR** vars, *var;
   SCIP_ROW* row; 
   int n, *L, i, nvars;

   printf("\nSeparaDif... construindo\n");
   assert(scip != NULL);
   
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);
   
   nvars = SCIPprobdataGetNVars(probdata);
   vars = SCIPprobdataGetVars(probdata);
   int kc = SCIPprobdataGetKc(probdata);

   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, conshdlr, "delta", kc,  -SCIPinfinity(scip) , FALSE, FALSE, TRUE) ); 

   SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

   /* inclui coeficientes na nova restricao : sum_{i in sol} x_i <= |S| - 1 */
   for(i=0;i<nvars;i++){
     var = vars[i];
     solval = SCIPgetSolVal(scip, sol, var);
     if(solval>0.0001){
       SCIP_CALL( SCIPaddVarToRow(scip, row, vars[i], 1.0) );
     }
   }

   SCIP_CALL( SCIPflushRowExtensions(scip, row) );

   if( SCIPisCutEfficacious(scip, sol, row) )
     {
       SCIP_Bool infeasible;
       SCIP_CALL( SCIPaddCut(scip, sol, row, FALSE, &infeasible) );
       if ( infeasible ){
	 *result = SCIP_CUTOFF;
	 printf("\n***** problema ficou inviavel com o corte (poda o no)\n");
       }
       else{
	 *result = SCIP_SEPARATED;
	 printf("\n***** corte foi separado\n");
       }
       SCIProwPrint(row, NULL);
     }
   else{
	 printf("\n***** efficacious...\n");     
   }
   SCIP_CALL( SCIPreleaseRow(scip, &row) );
   //   SCIP_CALL( SCIPwriteOrigProblem(scip, "mochila2.lp", "lp", FALSE) );
   
   return SCIP_OKAY;
}


/*
 * Linear constraint upgrading
 */

#ifdef LINCONSUPGD_PRIORITY
/** tries to upgrade a linear constraint into a xyz constraint */
static
SCIP_DECL_LINCONSUPGD(linconsUpgdXyz)
{  /*lint --e{715}*/
   SCIP_Bool upgrade;

   assert(upgdcons != NULL);

   /* check, if linear constraint can be upgraded to xyz constraint */
   upgrade = FALSE;
   /* TODO: put the constraint's properties here, in terms of the statistics given by nposbin, nnegbin, ... */

   if( upgrade )
   {
      SCIPdebugMessage("upgrading constraint <%s> to xyz constraint\n", SCIPconsGetName(cons));

      /* create the bin Xyz constraint (an automatically upgraded constraint is always unmodifiable) */
      assert(!SCIPconsIsModifiable(cons));
      SCIP_CALL( SCIPcreateConsXyz(scip, upgdcons, SCIPconsGetName(cons), nvars, vars, vals, lhs, rhs,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
            SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
            SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
            SCIPconsIsStickingAtNode(cons)) );
   }

   return SCIP_OKAY;
}
#endif


/*
 * Callback methods of constraint handler
 */

/* TODO: Implement all necessary constraint handler methods. The methods with #if 0 ... #else #define ... are optional */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define conshdlrCopyXyz NULL
#endif

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_CONSFREE(consFreeXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consFreeXyz NULL
#endif


/** initialization method of constraint handler (called after problem was transformed) */
#if 0
static
SCIP_DECL_CONSINIT(consInitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitXyz NULL
#endif


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_CONSEXIT(consExitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitXyz NULL
#endif


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#if 0
static
SCIP_DECL_CONSINITPRE(consInitpreXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitpreXyz NULL
#endif


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#if 0
static
SCIP_DECL_CONSEXITPRE(consExitpreXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitpreXyz NULL
#endif


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_CONSINITSOL(consInitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitsolXyz NULL
#endif


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_CONSEXITSOL(consExitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitsolXyz NULL
#endif


/** frees specific constraint data */
#if 0
static
SCIP_DECL_CONSDELETE(consDeleteXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeleteXyz NULL
#endif


/** transforms constraint data into data belonging to the transformed problem */
#if 0
static
SCIP_DECL_CONSTRANS(consTransXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consTransXyz NULL
#endif


/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
#if 0
static
SCIP_DECL_CONSINITLP(consInitlpXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitlpXyz NULL
#endif


/** separation method of constraint handler for LP solutions */
//#if 0
static
SCIP_DECL_CONSSEPALP(consSepalpXyz)
{  /*lint --e{715}*/
   printf("method of xyz constraint handler (conssepalp) implemented!\n");
   // SCIPABORT(); /*lint --e{527}*/
   SCIP_Bool found;
   int imaior, imenor, tamanho;
   
   *result = SCIP_FEASIBLE;

   found = violaPesoMedio(scip, NULL, &imaior, &imenor, &tamanho);
   // if a restricao violada was found, we generate a cut constraint 
   if( found ){

     deltaConstraint(scip, conshdlr,NULL, tamanho, result);
     //     *result = SCIP_INFEASIBLE; // SCIP_CONSADDED ou SEPARATED
   }
  
   return SCIP_OKAY;
}
//#else
//#define consSepalpXyz NULL
//#endif


/** separation method of constraint handler for arbitrary primal solutions */
//#if 0
static
SCIP_DECL_CONSSEPASOL(consSepasolXyz)
{  /*lint --e{715}*/
   printf("method of xyz constraint handler (conssepasol) implemented!\n");
   // SCIPABORT(); /*lint --e{527}*/
   SCIP_Bool found;
   int imaior, imenor, tamanho;
   
   *result = SCIP_FEASIBLE;

   found = violaPesoMedio(scip, sol, &imaior, &imenor, &tamanho);
   // if a restricao violada was found, we generate a cut constraint 
   if( found ){

     deltaConstraint(scip, conshdlr, sol, tamanho, result);
     //     *result = SCIP_INFEASIBLE; // SCIP_CONSADDED ou SEPARATED
   }

   
   return SCIP_OKAY;
}
//#else
//#define consSepasolXyz NULL
//#endif


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpXyz)
{  /*lint --e{715}*/
   printf("method of xyz constraint handler (consenfolp)... wait\n");
   //   SCIPABORT(); /*lint --e{527}*/
   SCIP_Bool found;
   int imaior, imenor, tamanho;
   
   *result = SCIP_FEASIBLE;

   found = violaPesoMedio(scip, NULL, &imaior, &imenor, &tamanho);
   // if a restricao violada was found, we generate a cut constraint 
   if( found ){
     
     deltaConstraint(scip, conshdlr, NULL, tamanho, result);
     //     *result = SCIP_INFEASIBLE; // SCIP_CONSADDED ou SEPARATED
   }
   
   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsXyz)
{  /*lint --e{715}*/
   printf("method of xyz constraint handler (consenfops)... wait \n");
   //   SCIPABORT(); /*lint --e{527}*/
   int imaior, imenor, tamanho;

   SCIP_Bool found;
   
   *result = SCIP_FEASIBLE;

   found = violaPesoMedio(scip, NULL, &imaior, &imenor, &tamanho);
   // if a restricao violada was found, we generate a cut constraint 
   if( found ){
     *result = SCIP_INFEASIBLE;
   }
   
   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckXyz)
{  /*lint --e{715}*/
   printf("method of xyz constraint handler (consCheckxyz)... wait\n");
   //   SCIPABORT(); /*lint --e{527}*/

   SCIP_Bool found;
   int imaior, imenor, tamanho;
   
   *result = SCIP_FEASIBLE;

   found = violaPesoMedio(scip, sol, &imaior, &imenor, &tamanho);
   // if a restricao violada was found, we generate a cut constraint 
   if( found ){
     *result = SCIP_INFEASIBLE;
   }
   
   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
#if 0
static
SCIP_DECL_CONSPROP(consPropXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPropXyz NULL
#endif


/** presolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSPRESOL(consPresolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPresolXyz NULL
#endif


/** propagation conflict resolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSRESPROP(consRespropXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consRespropXyz NULL
#endif


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockXyz)
{  /*lint --e{715}*/
   printf("method of xyz constraint handler not implemented yet.. conslock (not done)\n");
   //   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSACTIVE(consActiveXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consActiveXyz NULL
#endif


/** constraint deactivation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDEACTIVE(consDeactiveXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeactiveXyz NULL
#endif


/** constraint enabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSENABLE(consEnableXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consEnableXyz NULL
#endif


/** constraint disabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDISABLE(consDisableXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDisableXyz NULL
#endif

/** variable deletion of constraint handler */
#if 0
static
SCIP_DECL_CONSDELVARS(consDelvarsXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDelvarsXyz NULL
#endif


/** constraint display method of constraint handler */

static
SCIP_DECL_CONSPRINT(consPrintXyz)
{  /*lint --e{715}*/
   printf("method of xyz constraint handler (consprint) partial implemented\n");
   //   SCIPABORT(); /*lint --e{527}*/

   SCIPinfoMessage(scip, file, "constraint xyz para pesos medios...\n");

   return SCIP_OKAY;
}


/** constraint copying method of constraint handler */
#if 0
static
SCIP_DECL_CONSCOPY(consCopyXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consCopyXyz NULL
#endif


/** constraint parsing method of constraint handler */
#if 0
static
SCIP_DECL_CONSPARSE(consParseXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consParseXyz NULL
#endif


/** constraint method of constraint handler which returns the variables (if possible) */
#if 0
static
SCIP_DECL_CONSGETVARS(consGetVarsXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetVarsXyz NULL
#endif

/** constraint method of constraint handler which returns the number of variables (if possible) */
#if 0
static
SCIP_DECL_CONSGETNVARS(consGetNVarsXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetNVarsXyz NULL
#endif

/** constraint handler method to suggest dive bound changes during the generic diving algorithm */
#if 0
static
SCIP_DECL_CONSGETDIVEBDCHGS(consGetDiveBdChgsXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetDiveBdChgsXyz NULL
#endif


/*
 * constraint specific interface methods
 */

/** creates the handler for xyz constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrXyz(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* create xyz constraint handler data */
   conshdlrdata = NULL;
   /* TODO: (optional) create constraint handler specific data here */

   conshdlr = NULL;

   /* include constraint handler */
#if 0
   /* use SCIPincludeConshdlr() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS,
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_NEEDSCONS,
         CONSHDLR_PROP_TIMING, CONSHDLR_PRESOLTIMING,
         conshdlrCopyXyz,
         consFreeXyz, consInitXyz, consExitXyz,
         consInitpreXyz, consExitpreXyz, consInitsolXyz, consExitsolXyz,
         consDeleteXyz, consTransXyz, consInitlpXyz,
         consSepalpXyz, consSepasolXyz, consEnfolpXyz, consEnfopsXyz, consCheckXyz,
         consPropXyz, consPresolXyz, consRespropXyz, consLockXyz,
         consActiveXyz, consDeactiveXyz,
         consEnableXyz, consDisableXyz, consDelvarsXyz,
         consPrintXyz, consCopyXyz, consParseXyz,
         consGetVarsXyz, consGetNVarsXyz, consGetDiveBdChgsXyz, conshdlrdata) );
#else
   /* use SCIPincludeConshdlrBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpXyz, consEnfopsXyz, consCheckXyz, consLockXyz,
         conshdlrdata) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrActive(scip, conshdlr, consActiveXyz) );
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyXyz, consCopyXyz) );
   SCIP_CALL( SCIPsetConshdlrDeactive(scip, conshdlr, consDeactiveXyz) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteXyz) );
   SCIP_CALL( SCIPsetConshdlrDelvars(scip, conshdlr, consDelvarsXyz) );
   SCIP_CALL( SCIPsetConshdlrDisable(scip, conshdlr, consDisableXyz) );
   SCIP_CALL( SCIPsetConshdlrEnable(scip, conshdlr, consEnableXyz) );
   SCIP_CALL( SCIPsetConshdlrExit(scip, conshdlr, consExitXyz) );
   SCIP_CALL( SCIPsetConshdlrExitpre(scip, conshdlr, consExitpreXyz) );
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolXyz) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeXyz) );
   SCIP_CALL( SCIPsetConshdlrGetDiveBdChgs(scip, conshdlr, consGetDiveBdChgsXyz) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsXyz) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsXyz) );
   SCIP_CALL( SCIPsetConshdlrInit(scip, conshdlr, consInitXyz) );
   SCIP_CALL( SCIPsetConshdlrInitpre(scip, conshdlr, consInitpreXyz) );
   SCIP_CALL( SCIPsetConshdlrInitsol(scip, conshdlr, consInitsolXyz) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpXyz) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseXyz) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolXyz, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintXyz) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropXyz, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
         CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropXyz) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpXyz, consSepasolXyz, CONSHDLR_SEPAFREQ, CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransXyz) );
#endif

#ifdef LINCONSUPGD_PRIORITY
   if( SCIPfindConshdlr(scip,"linear") != NULL )
   {
      /* include the linear constraint upgrade in the linear constraint handler */
      SCIP_CALL( SCIPincludeLinconsUpgrade(scip, linconsUpgdXyz, LINCONSUPGD_PRIORITY, CONSHDLR_NAME) );
   }
#endif

   /* add xyz constraint handler parameters */
   /* TODO: (optional) add constraint handler specific parameters with SCIPaddTypeParam() here */

   printf("\nmanipulador incluido...\n");
   return SCIP_OKAY;
}

/** creates and captures a xyz constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsXyz(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_Real*            coefs,              /**< array with coefficients of constraint entries */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs,                /**< right hand side of constraint */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   /* TODO: (optional) modify the definition of the SCIPcreateConsXyz() call, if you don't need all the information */

   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   printf("method of xyz constraint (create) implemented \n");
   //   SCIPABORT(); /*lint --e{527} --e{715}*/

   /* find the xyz constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("xyz constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   consdata = NULL;
   /* TODO: create and store constraint specific data here */

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}

/** creates and captures a xyz constraint with all its constraint flags set to their
 *  default values
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicXyz(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_Real*            coefs,              /**< array with coefficients of constraint entries */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs                 /**< right hand side of constraint */
   )
{
   SCIP_CALL( SCIPcreateConsXyz(scip, cons, name, nvars, vars, coefs, lhs, rhs,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}
