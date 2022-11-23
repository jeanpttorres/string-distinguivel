#ifndef __PARAMETERS_DSSP_H__
#define __PARAMETERS_DSSP_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif
typedef enum{
             //   int heur; /* Default = 0 (none). 1 = RA, 2 = BCPA, 3 = VBPL, 4 = BLPL, 5 = SBPL*/
             NONE=0, RA, BCPA, VBPL, BLPL, SBPL, TOTALHEURTYPE
} heurT;
typedef enum{
             CSP=0, FSP=1, DSP=2
} problemT;

typedef enum{
             FIRST_IMPROV=1, BEST_IMPROV=2, LP_IMPROV=3
} localSearchT;

typedef enum{
             NORMAL_SHAKING=1, LP_SHAKING=2
} shakingT;

extern const char*heurName[];
extern const char*problemName[];

typedef struct{
   // global settings
   int time_limit; /* limit of execution time (in sec). Default = 600 (-1: unlimited) */
   int display_freq; /* frequency to display information about B&B enumeration. Default = -1 (-1: never) */
   int nodes_limit; /* limit of nodes to B&B procedure. Default = -1: unlimited (1: onlyrootnode) */
   int scipdefault; /* should use scip default settings? Default = 0.*/

   // presolving
   int presolve_maxrounds; /* Default = 0: never */
   int presolve_maxrestarts; /* Default = 0: never */

   // branching rules
   int pscost_priority; /* priority to branching using pscost. Default = 1000000 */

   // propagation
   int propag; /* Default = 0: disabled */
   int propag_maxrounds; /* Default = 0: never */
   int propag_maxroundsroot; /* Default = 0: never */

   // node selecting

   // problem 
   int substring; /* CSSP, FSSP, DSSSP? Default = 0. 0 =  STRING, 1 = SUB*/
   int model_novo; /* Default = 0 */
   int problem; /* 0=CSP, 1=FSP, 2=DSP. Default = DSP */

   // heuristics
   int heur; /* Default = 0 (none). 1 = RA, 2 = BCPA, 3 = VBPL, 4 = BLPL, 5 = SBPL*/
   int local_search; /* Type of local search. Default = 0 (none) . First = 1, Best = 2, LPImprov = 3 */
   int shaking; /* Type of shaking. Default = 0 (none) . Normal = 1, LP_SHAK = 2 */
   int vnstime; /* Max time for VNS loop. Default = 5 seconds */
   int ngsize; /* Size of neighbourhood intervals. Default = 3 ([0,1/3), [1/3,2/3) and [2/3,1]) */
   int vnssecond; /* Default = 0. ??? */
   int semy; /* should not consider variable y to lpshaking? Default = 1*/
   int bcpa_rounding;
   int ra_novo;
   int ra_reduce;
   int ra_merge;
   int bcpatime_limit;

   // parameter stamp
   char* parameter_stamp;

} parametersT;

int setParameters(int argc, char** argv, parametersT* Param);

extern parametersT param;
extern char* instance_filename;
extern const char* output_path;
extern char input_path[SCIP_MAXSTRLEN];


#ifdef __cplusplus
}
#endif

#endif


