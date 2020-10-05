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

/**@file   probdata_dssp.h
 * @brief  Problem data for dssp problem
 * @author Timo Berthold
 * @author Stefan Heinz
 *
 * This file handles the main problem data used in that project. For more details see \ref PROBLEMDATA page.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PROBDATA_DSSP__
#define __SCIP_PROBDATA_DSSP__

#include "scip/scip.h"

#define EPSILON 0.000001

/**@name Constraints quantity
 * 
 * @{
 */

#define N_CONS_1  string_size * sc_size         /**< Conjunto de restrições (1) */
#define N_CONS_2  (string_size - 1) * (sum_sc[sc_size] - sc_size)                     /**< Conjunto de restrições (2) */
#define N_CONS_3  (sc_size - 1) * alpha_size * string_size     /**< Conjunto de restrições (3) */
#define N_CONS_4  sc_size             /**< Conjunto de restrições (4) */
#define N_CONS_5  sum_sf[sf_size]
#define N_CONS_6  alpha_size * sc_size
/**@} */

/**@name Variables quantity
 * 
 * @{
 */

#define N_VARS_X  alpha_size * string_size * sum_sc[sc_size] * sc_size     /**< Variáveis xij */
#define N_VARS_D  2          /**< Variáveis dc e df */

/**@} */

/**@name Constraints position
 * 
 * @{
 */
 

#define N_ROWS_1(i, j)        (i) * string_size + (j)                                /**< Conjunto de restrições (1) */
#define N_ROWS_2(i, j, k)     CONS_1 + (string_size-1)*(sum_sc[(i)] - (i)) + (j)*(sum_sc[(i) + 1] - sum_sc[(i)] - 1) + (k)           /**< Conjunto de restrições (2) */
#define N_ROWS_3(i, j, k)     CONS_1 + CONS_2 + (i-1)*(alpha_size*string_size) + (j)*alpha_size + (k)  /**< Conjunto de restrições (3) */
#define N_ROWS_4(i)           CONS_1 + CONS_2 + CONS_3 + (i)     /**< Conjunto de restrições (4) */
#define N_ROWS_5(i, j)        CONS_1 + CONS_2 + CONS_3 + CONS_4 + sum_sf[(i)] + (j)
#define N_ROWS_6(i, l)        CONS_1 + CONS_2 + CONS_3 + CONS_4 + (sc_size * alpha_size * i) + (l)

/**@} */




/**@name Constraints quantity
 *
 * @{
 */

#define CONS_1  string_size         /**< Conjunto de restrições (1) */
#define CONS_2  sum_sc[sc_size]     /**< Conjunto de restrições (2) */
#define CONS_3  sum_sf[sf_size]     /**< Conjunto de restrições (3) */
#define CONS_4  sc_size             /**< Conjunto de restrições (4) */
#define CONS_6  alpha_size * sc_size


/**@} */

/**@name Variables quantity
 *
 * @{
 */

#define VARS_X  alpha_size * string_size     /**< Variáveis xij */
#define VARS_Y  sum_sc[sc_size]              /**< Váriaveis yki */
#define VARS_D  2                            /**< Variáveis dc e df */


/**@} */

/**@name Constraints position
 *
 * @{
 */


#define ROWS_1(i)       i                                /**< Conjunto de restrições (1) */
#define ROWS_2(i, k)    CONS_1 + sum_sc[i] + k           /**< Conjunto de restrições (2) */
#define ROWS_3(i, k)    CONS_1 + CONS_2 + sum_sf[i] + k  /**< Conjunto de restrições (3) */
#define ROWS_4(i)       CONS_1 + CONS_2 + CONS_3 + i     /**< Conjunto de restrições (4) */

/**@} */

/** sets up the problem data */
extern
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
    );


/** returns size of sc */
extern
int SCIPprobdataGetSc_size(
    SCIP_PROBDATA*      probdata            /**< problem data */
    );

/** returns size of sf */
extern
int SCIPprobdataGetSf_size(
    SCIP_PROBDATA*      probdata            /**< problem data */
    );

/** returns size of alphabet */
extern
int SCIPprobdataGetAlpha_size(
    SCIP_PROBDATA*      probdata            /**< problem data */
    );

/** returns size of string */
extern
int SCIPprobdataGetString_size(
    SCIP_PROBDATA*      probdata            /**< problem data */
    );

/** returns kc */
extern
int SCIPprobdataGetKc(
    SCIP_PROBDATA*      probdata            /**< problem data */
    );

/** returns kf */
extern
int SCIPprobdataGetKf(
    SCIP_PROBDATA*      probdata            /**< problem data */
    );

/** returns array of strings Sc */
extern
char** SCIPprobdataGetSc(
    SCIP_PROBDATA*      probdata            /**< problem data */
    );

/** returns array of strings Sf */
extern
char** SCIPprobdataGetSf(
    SCIP_PROBDATA*      probdata            /**< problem data */
    );

/** returns sum_sc */
extern
int* SCIPprobdataGetSum_sc(
    SCIP_PROBDATA*      probdata            /**< problem data */
    );

/** returns sum_sf */
extern
int* SCIPprobdataGetSum_sf(
    SCIP_PROBDATA*      probdata            /**< problem data */
    );

/** returns alphabet */
extern
char* SCIPprobdataGetAlphabet(
    SCIP_PROBDATA*      probdata            /**< problem data */
    );

/** returns array of variables */
extern
SCIP_VAR** SCIPprobdataGetVars(
    SCIP_PROBDATA*      probdata            /**< problem data */
    );

/** returns number of variables */
extern
int SCIPprobdataGetNVars(
    SCIP_PROBDATA*      probdata            /**< problem data */
    );

/** returns number of constrains */
extern
int SCIPprobdataGetNcons(
    SCIP_PROBDATA*      probdata            /**< problem data */
    );

/** returns array of constrains */
extern
SCIP_CONS** SCIPprobdataGetConss(
    SCIP_PROBDATA*      probdata            /**< problem data */
    );

/** returns the problem name */
extern
const char* SCIPprobdataGetProbname(
    SCIP_PROBDATA*      probdata            /**< problem data */
    );

/** returns the position of the character in the alphabet */
extern
int value_symbol (
    char*               alphabet,           /**< alphabet */
    int                 alpha_size,         /**< size of the alphabet */
    char                symbol              /**< symbol to look in the alphabet */
    );

/** returns the character at the value position in the alphabet */
extern
char symbol_value (
    char*               alphabet,           /**< alphabet */
    int                 alpha_size,         /**< size of the alphabet */
    int                 value               /**< symbol position to look in the alphabet */
    );

#endif
