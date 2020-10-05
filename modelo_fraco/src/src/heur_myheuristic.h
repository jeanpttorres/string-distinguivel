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

/**@file   heur_myheuristic.h
 * @ingroup PRIMALHEURISTICS
 * @brief  myheuristic primal heuristic
 * @author Tobias Achterberg
 *
 * template file for primal heuristic plugins
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_myheuristic_H__
#define __SCIP_HEUR_myheuristic_H__

#include "scip/scip.h"


#ifdef __cplusplus
extern "C" {
#endif

/** creates the myheuristic primal heuristic and includes it in SCIP */
EXTERN
SCIP_RETCODE SCIPincludeHeurmyheuristic(
   	SCIP*                 	scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
