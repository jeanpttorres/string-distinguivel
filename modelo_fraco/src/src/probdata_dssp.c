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

/**@file   probdata_dssp.c
 * @brief  Problem data for dssp problem
 * @author Timo Berthold
 * @author Stefan Heinz
 *
 * This file handles the main problem data used in that project. For more details see \ref PROBLEMDATA page.
 *
 * @page PROBLEMDATA Main problem data
 *
 * The problem data is accessible in all plugins. The function SCIPgetProbData() returns the pointer to that
 * structure. We use this data structure to store all the information of the dssp problem. Since this structure is
 * not visible in the other plugins, we implemented setter and getter functions to access this data. The problem data
 * structure SCIP_ProbData is shown below.
 *
 * \code
 *  ** @brief Problem data which is accessible in all places
 *  *
 *  *   This problem data is used to store the input of the dssp instance, all variables which are created, and all
 *  *   constraints.
 *  *
 * \endcode
 *
 * The function SCIPprobdataCreate(), which is called in the \ref reader_dssp.c "reader plugin" after the input file was
 * parsed, initializes the problem data structure and creates the problem in the SCIP environment. For this, it creates
 * for each item of the binpacking problem one set covering constraint and creates an initial set of variables for the
 * packings. Note that the set covering constraints have to have the <code>modifiable</code>-flag set to TRUE. This is
 * necessary to tell the solver that these constraints are not completed yet. This means, during the search new
 * variables/packings might be added.  The solver needs this information because certain reductions are not allowed.
 * See the body of the function SCIPprobdataCreate() for more details.
 *
 * A list of all interface methods can be found in probdata_dssp.h.
 **/

/*

Modificações:

2015-10-09: Inclusão de heuristica primal para povar a base inicial. A melhor solucao eh usada para setar uma limitante primal (UB) dentro do SCIP

*/

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/*
#define SCIP_DEBUG
*/
#include <string.h>
#include "probdata_dssp.h"

#include "scip/cons_linear.h"
#include "scip/cons_setppc.h"
#include "scip/scip.h"

/** @brief Problem data which is accessible in all places
 *
 * This problem data is used to store the input of the dssp, all variables which are created, and all
 * constrsaints.
 */
struct SCIP_ProbData
{
    SCIP_VAR**          vars;           /**< all exiting variables in the problem */
    SCIP_CONS**         conss;          /**< all constraints */
    int                 nvars;          /**< number of generated variables */
    int                 varssize;       /**< size of the variable array */
    int                 ncons;          /**< number of constraints */
    int                 sc_size;        /**< number of strings in Sc */
    int                 sf_size;        /**< number of strings in Sf */
    int                 alpha_size;     /**< size of the alphabet */
    int                 string_size;    /**< size of the strings */
    int                 kc;             /**< dc bound */
    int                 kf;             /**< df bound */
    char**              sc;             /**< array of strings Sc */
    char**              sf;             /**< array of strings Sf */
    int*                sum_sc;         /**< array of length sum of Sc strings */
    int*                sum_sf;         /**< array of length sum of Sf strings */
    char*               alphabet;       /**< array of alphabet characters */
    const char*         probname;       /**< filename of the instance */
};

/**@name Event handler properties
 *
 * @{
 */

#define EVENTHDLR_NAME         "addedvar"
#define EVENTHDLR_DESC         "event handler for catching added variables"

/**@} */

/**@name Callback methods of event handler
 *
 * @{
 */

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecAddedVar)
{  /*lint --e{715}*/
    assert(eventhdlr != NULL);
    assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
    assert(event != NULL);
    assert(SCIPeventGetType(event) == SCIP_EVENTTYPE_VARADDED);

    /*SCIPdebugMessage("exec method of event handler for added variable to probdata\n");*/

    /* add new variable to probdata */
    return SCIP_OKAY;
}

/**@} */


/**@name Local methods
 *
 * @{
 */

/** creates problem data */
static
SCIP_RETCODE probdataCreate(
    SCIP*               scip,           /**< SCIP data structure */
    SCIP_PROBDATA**     probdata,       /**< pointer to problem data */
    SCIP_VAR**          vars,           /**< all exist variables */
    SCIP_CONS**         conss,          /**< all pmr constraints */
    int                 nvars,          /**< number of variables */
    int                 ncons,          /**< number of constraints */
    int                 sc_size,        /**< number of strings in Sc */
    int                 sf_size,        /**< number of strings in Sf */
    int                 alpha_size,     /**< size of the alphabet */
    int                 string_size,    /**< size of the strings */
    int                 kc,             /**< dc bound */
    int                 kf,             /**< df bound */
    char**              sc,             /**< array of strings Sc */
    char**              sf,             /**< array of strings Sf */
    int*                sum_sc,         /**< array of length sum of Sc strings */
    int*                sum_sf,         /**< array of length sum of Sf strings */
    char*               alphabet,       /**< array of alphabet characters */
    const char*         probname
    )
{
    assert(scip != NULL);
    assert(probdata != NULL);

    /* allocate memory */
    SCIP_CALL(SCIPallocMemory(scip, probdata));

    if(nvars > 0)
    {
        /* copy variable array */
        (*probdata)->vars = vars;
        SCIP_CALL(SCIPduplicateMemoryArray(scip, &(*probdata)->vars, vars, nvars)); /* NEEDED for transformed problem*/
    }
    else {
        (*probdata)->vars = NULL;
    }

    /* duplicate arrays */
    SCIP_CALL(SCIPduplicateMemoryArray(scip, &(*probdata)->conss, conss, ncons)); /* NEEDED for transformed problem */

    (*probdata)->nvars = nvars;
    (*probdata)->varssize = nvars;
    (*probdata)->ncons = ncons;
    (*probdata)->sc_size = sc_size;
    (*probdata)->sf_size = sf_size;
    (*probdata)->alpha_size = alpha_size;
    (*probdata)->string_size = string_size;
    (*probdata)->kc = kc;
    (*probdata)->kf = kf;
    (*probdata)->sc = sc;
    (*probdata)->sf = sf;
    (*probdata)->sum_sc = sum_sc;
    (*probdata)->sum_sf = sum_sf;
    (*probdata)->alphabet = alphabet;
    (*probdata)->probname = probname;

    return SCIP_OKAY;
}

/** frees the memory of the given problem data */
static
SCIP_RETCODE probdataFree(
    SCIP*                 scip,               /**< SCIP data structure */
    SCIP_PROBDATA**       probdata            /**< pointer to problem data */
    )
{

    int i;

    assert(scip != NULL);
    assert(probdata != NULL);

    /* release all variables */
    for(i = 0; i < (*probdata)->nvars; ++i)
    {
        SCIP_CALL(SCIPreleaseVar(scip, &(*probdata)->vars[i]));
    }

    /* release all constraints */
    for(i = 0; i < (*probdata)->ncons; ++i)
    {
        SCIP_CALL(SCIPreleaseCons(scip, &(*probdata)->conss[i]));
    }

    /* free memory of arrays */
    SCIPfreeMemoryArray(scip, &(*probdata)->vars);
    SCIPfreeMemoryArray(scip, &(*probdata)->conss);
    for (i = 0; i < (*probdata)->sc_size; i++) {
        SCIPfreeMemoryArray(scip, &(*probdata)->sc[i]);
    }
    SCIPfreeMemoryArray(scip, &(*probdata)->sc);

    for (i = 0; i < (*probdata)->sf_size; i++) {
        SCIPfreeMemoryArray(scip, &(*probdata)->sf[i]);
    }
    SCIPfreeMemoryArray(scip, &(*probdata)->sum_sc);

    SCIPfreeMemoryArray(scip, &(*probdata)->sum_sf);

    SCIPfreeMemoryArray(scip, &(*probdata)->alphabet);

    /* free probdata */
    SCIPfreeMemory(scip, probdata);

    return SCIP_OKAY;
}


/**@} */


/**@name Callback methods of problem data
 *
 * @{
 */

/** frees user data of original problem (called when the original problem is freed) */
static
SCIP_DECL_PROBDELORIG(probdelorigdssp)
{
    /*SCIPdebugMessage("free original problem data\n");*/

    SCIP_CALL(probdataFree(scip, probdata));

    return SCIP_OKAY;
}

/** creates user data of transformed problem by transforming the original user problem data
 *  (called after problem was transformed) @@@ ?*/
static
SCIP_DECL_PROBTRANS(probtransdssp)
{
    /* create transform probdata */
    SCIP_CALL(probdataCreate(scip, targetdata, sourcedata->vars, sourcedata->conss, sourcedata->nvars, sourcedata->ncons, sourcedata->sc_size, sourcedata->sf_size, sourcedata->alpha_size, sourcedata->string_size, sourcedata->kc, sourcedata->kf, sourcedata->sc, sourcedata->sf, sourcedata->sum_sc, sourcedata->sum_sf, sourcedata->alphabet, sourcedata->probname));

    /* transform all constraints */
    SCIP_CALL(SCIPtransformConss(scip, (*targetdata)->ncons, (*targetdata)->conss, (*targetdata)->conss));

    /* transform all variables */
    SCIP_CALL(SCIPtransformVars(scip, (*targetdata)->nvars, (*targetdata)->vars, (*targetdata)->vars));

    return SCIP_OKAY;
}

/** frees user data of transformed problem (called when the transformed problem is freed) */
static
SCIP_DECL_PROBDELTRANS(probdeltransdssp)
{
    /*SCIPdebugMessage("free transformed problem data\n");*/

    SCIP_CALL(probdataFree(scip, probdata));

    return SCIP_OKAY;
}

/** solving process initialization method of transformed data (called before the branch and bound process begins) */
static
SCIP_DECL_PROBINITSOL(probinitsoldssp)
{
    SCIP_EVENTHDLR* eventhdlr;

    assert(probdata != NULL);

    /* catch variable added event */
    eventhdlr = SCIPfindEventhdlr(scip, "addedvar");
    assert(eventhdlr != NULL);

    SCIP_CALL(SCIPcatchEvent(scip, SCIP_EVENTTYPE_VARADDED, eventhdlr, NULL, NULL));

    return SCIP_OKAY;
}

/** solving process deinitialization method of transformed data (called before the branch and bound data is freed) */
static
SCIP_DECL_PROBEXITSOL(probexitsoldssp)
{
    SCIP_EVENTHDLR* eventhdlr;

    assert(probdata != NULL);

    /* drop variable added event */
    eventhdlr = SCIPfindEventhdlr(scip, "addedvar");
    assert(eventhdlr != NULL);

    SCIP_CALL(SCIPdropEvent(scip, SCIP_EVENTTYPE_VARADDED, eventhdlr, NULL, -1));


    return SCIP_OKAY;
}

/**@} */


/**@name Interface methods
 *
 * @{
 */

/** sets up the problem data */
SCIP_RETCODE SCIPprobdataCreate(
    SCIP*               scip,           /**< SCIP data structure */
    const char*         probname,       /**< problem name */
    int                 sc_size,        /**< number of strings in Sc */
    int                 sf_size,        /**< number of strings in Sf */
    int                 alpha_size,     /**< size of the alphabet */
    int                 string_size,    /**< size of the strings */
    int                 kc,             /**< dc bound */
    int                 kf,             /**< df bound */
    char**              sc,             /**< array of strings Sc */
    char**              sf,             /**< array of strings Sf */
    int*                sum_sc,         /**< array of length sum of Sc strings */
    int*                sum_sf,         /**< array of length sum of Sf strings */
    char*               alphabet        /**< array of alphabet characters */
    )
{
    int i, j, k, l, ncons, nvars;
    char name[SCIP_MAXSTRLEN];
    SCIP_PROBDATA* probdata;
    SCIP_CONS** conss;
    SCIP_VAR** vars, *var;

    assert(scip != NULL);

    /* create event handler if it does not exist yet */
    if(SCIPfindEventhdlr(scip, EVENTHDLR_NAME) == NULL)
    {
        SCIP_CALL(SCIPincludeEventhdlrBasic(scip, NULL, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecAddedVar, NULL));
    }

    /* create problem in SCIP and add non-NULL callbacks via setter functions */
    SCIP_CALL(SCIPcreateProbBasic(scip, probname));

    SCIP_CALL(SCIPsetProbDelorig(scip, probdelorigdssp));
    SCIP_CALL(SCIPsetProbTrans(scip, probtransdssp));
    SCIP_CALL(SCIPsetProbDeltrans(scip, probdeltransdssp));
    SCIP_CALL(SCIPsetProbInitsol(scip, probinitsoldssp));
    SCIP_CALL(SCIPsetProbExitsol(scip, probexitsoldssp));

    /* set objective sense */
    SCIP_CALL(SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE));

    /* tell SCIP that the objective will be always integral */
    SCIP_CALL(SCIPsetObjIntegral(scip));

    /* Number of constraints */
    SCIP_CALL(SCIPallocBufferArray(scip, &conss, CONS_1 + CONS_2 + CONS_3 + CONS_4));

    /* Number of variables */
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
    
    // /* (4) cosntraint - disjunctive constraint */
    // if (definedVarsSize > 0) {
    //     ncons = i;
    //     for (; i < ncons + CONS_5; i++) {
    //         SCIP_CALL(SCIPcreateConsBasicLinear (scip, &conss[i], "(5)", 0, NULL, NULL, 1.0, SCIPinfinity(scip)));
    //         SCIP_CALL(SCIPaddCons(scip, conss[i]));
    //     }
    // }

    nvars = 0;
    ncons = i;
    
    //int bcpa = 0;

    /**
    * add variables xij
    */
    for(i = 0; i < string_size; i++) {
        for (j = 0; j < alpha_size; j++) {
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x_%d_%c", i, symbol_value(alphabet, alpha_size, j));
            SCIPdebugMessage("create variable %s\n", name);

            /* create a basic variable object */
            SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY));
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

    /* create a basic variable object */
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
    SCIPdebugMessage("create variable %s\n", name);

    /* create a basic variable object */
    SCIP_CALL(SCIPcreateVarBasic(scip, &var, name, (double) kf, (double) string_size, -1.0, SCIP_VARTYPE_INTEGER));
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

    /* create problem data */
    SCIP_CALL(probdataCreate(scip, &probdata, vars, conss, nvars, ncons, sc_size, sf_size, alpha_size, string_size, kc, kf, sc, sf, sum_sc, sum_sf, alphabet, probname));

    SCIP_CALL(SCIPwriteOrigProblem(scip, "dssp.lp", "lp", FALSE)); /* grava na saida padrao ou em file */

    /* set user problem data */
    SCIP_CALL(SCIPsetProbData(scip, probdata));

    /* free local buffer arrays */
    SCIPfreeBufferArray(scip, &conss);
    SCIPfreeBufferArray(scip, &vars);

    return SCIP_OKAY;
}

/** returns size of Sc */
int SCIPprobdataGetSc_size(
    SCIP_PROBDATA*        probdata            /**< problem data */
    )
{
    return probdata->sc_size;
}

/** returns Size of Sf */
int SCIPprobdataGetSf_size(
    SCIP_PROBDATA*        probdata            /**< problem data */
    )
{
    return probdata->sf_size;
}

/** returns size of alphabet */
int SCIPprobdataGetAlpha_size(
    SCIP_PROBDATA*        probdata            /**< problem data */
    )
{
    return probdata->alpha_size;
}

/** returns size of string */
int SCIPprobdataGetString_size(
    SCIP_PROBDATA*        probdata            /**< problem data */
    )
{
    return probdata->string_size;
}

/** returns dc bound */
int SCIPprobdataGetKc(
  SCIP_PROBDATA*        probdata            /**< problem data */
  )
{
  return probdata->kc;
}

/** returns df bound */
int SCIPprobdataGetKf(
    SCIP_PROBDATA*        probdata            /**< problem data */
    )
{
    return probdata->kf;
}

/** returns Sc */
char** SCIPprobdataGetSc(
    SCIP_PROBDATA*        probdata            /**< problem data */
    )
{
    return probdata->sc;
}

/** returns Sf */
char** SCIPprobdataGetSf(
    SCIP_PROBDATA*        probdata            /**< problem data */
    )
{
    return probdata->sf;
}

/** returns sum_sc */
int* SCIPprobdataGetSum_sc(
    SCIP_PROBDATA*        probdata            /**< problem data */
    )
{
    return probdata->sum_sc;
}

/** returns sum_sf */
int* SCIPprobdataGetSum_sf(
    SCIP_PROBDATA*        probdata            /**< problem data */
    )
{
    return probdata->sum_sf;
}

/** returns alphabet */
char* SCIPprobdataGetAlphabet(
    SCIP_PROBDATA*        probdata            /**< problem data */
    )
{
    return probdata->alphabet;
}

/** returns array of all variables itemed in the way they got generated */
SCIP_VAR** SCIPprobdataGetVars(
    SCIP_PROBDATA*        probdata            /**< problem data */
    )
{
    return probdata->vars;
}

/** returns number of variables */
int SCIPprobdataGetNVars(
    SCIP_PROBDATA*        probdata            /**< problem data */
    )
{
    return probdata->nvars;
}

/** returns array of set partitioning constrains */
SCIP_CONS** SCIPprobdataGetConss(
    SCIP_PROBDATA*        probdata            /**< problem data */
    )
{
    return probdata->conss;
}

/** returns array of set partitioning constrains */
int SCIPprobdataGetNcons(
    SCIP_PROBDATA*        probdata            /**< problem data */
    )
{
    return probdata->ncons;
}

/** returns Probname of the instance */
const char* SCIPprobdataGetProbname(
    SCIP_PROBDATA*        probdata            /**< problem data */
    )
{
    return probdata->probname;
}

/**
 * @function value_symbol
 * @param char symbol
 * @return int value which represents the numeric value of the symbol in the
 * given alphabet
*/
int value_symbol (char *alphabet, int alpha_size, char symbol) {
    for (int i = 0; i < alpha_size; i++) {
        if (alphabet[i] == symbol)
            return i;
    }

    printf("Character %c wasn't found in the given alphabet\n", symbol);
    return -1;
}

char symbol_value (char *alphabet, int alpha_size, int value) {
    if (value < 0 || value >= alpha_size)
        printf("Position %d invalid\n", value);
    return alphabet[value];
}



/**@} */
