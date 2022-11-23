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

/**
 * SBPL version: OK. Changes were made in lpimprovement and lpshaking.
 * BLPL version: BestImprovement has be included! Shaking was changed.  
 * TODO: parameters settings
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
//#include "reader_csp.h"
#include "parameters_dssp.h"
//#include "heur_myheuristic.h"
#include "utils.h"
#include "ra.h"
#include "bcpa.h"
#include "vns.h"

#define MAXINT 10e6

const char* heurName[]={"none", "RA", "BCPA", "VBPL", "BLPL", "SBPL"};

const char* problemName[]={"CSP", "FSP", "DSP"};

parametersT param;
// output path, path and filename of the instance
const char* output_path;
char input_path[SCIP_MAXSTRLEN];
char* instance_filename;
const char* current_path = ".";
clock_t antes, agora;
double tempoHeur, LB, UB;
double firstUB, heurValue;

//typedef struct{
//  char debug;
//} compilingOptionsT;

void splitfullfilename(char* fullfilename);
void configOutputName(char* name, char* probname, char* program);
SCIP_RETCODE printStatistic(SCIP* scip, char* outputname, int dc, int df);//, char* program);
SCIP_RETCODE configScip(SCIP** scip);
int readInstance(SCIP* scip, char* filename);
void printSol(SCIP* scip, char *filename, double tempo);
//void printSol(SCIP* scip, int* profit, int* custo, char* outputname, compilingOptionsT* compilingOptions);
void printHeurSol(SCIP* scip, char* filename, int* x, int dc, int df, double tempo);
void printCompilerSettings(FILE* file);
/** print best primal solution in file with extension ".sol"
 *  and in graphwiz format .gr
 */
int contFrac(SCIP* scip);

void splitfullfilename(char* fullfilename)
{
  // split fullfilename in path and filename
  instance_filename = strrchr(fullfilename, '/');
  if(instance_filename==NULL){
     instance_filename = fullfilename;
     input_path[0]='\0';
  }
  else{
     instance_filename++; // discard /  
     strncpy(input_path, fullfilename, instance_filename - fullfilename);
     input_path[instance_filename-fullfilename]='\0';
  }
}
//
void configOutputName(char* name, char* probname, char* program)
{
  char* program_filename;
 
 // remove path, if there exists on program name
 program_filename = strrchr(program, '/');
 if(program_filename==NULL){
   program_filename = program;
 }
 else{
   program_filename++; // discard /
 }

 (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s/%s-%s", output_path, instance_filename, program_filename);
 // append parameter stamp
 strcat(name, "-");
 strcat(name, param.parameter_stamp);
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

    strcpy(name, filename);
    strcat(name, ".lp.sol");

    file = fopen(name, "w");
    
    if(!file)
    {
        printf("\nProblem to create solution file: %s\n", name);
        return;
    }

    name[0] = '\0';

    fprintf(file, "Time: %lfs\n", tempo);
    fprintf(file, "Objective Value(dc - df) = %.0lf\n", SCIPgetPrimalbound(scip)); //SCIPsolGetOrigObj(bestSolution));
    fprintf(file, "Target String(x) = ");

    for(i = 0; i < string_size; i++) {
        for (j = 0; j < alpha_size; j++) {
            var = vars[i * alpha_size + j];
            solval = SCIPgetSolVal(scip, bestSolution, var);
            if(solval < EPSILON)
            {
                fprintf(file, "%c", symbol_value(alphabet, alpha_size, j));
            }
            //            else{
            //                fprintf(file, "-");
            //            }
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

    printCompilerSettings(file);
    fclose(file);
}
void printCompilerSettings(FILE* file)
{
    struct tm * ct;
    const time_t t = time(NULL);

  // save compiling and executions settings
    
#if defined DEBUG || defined DEBUG_VNS || defined DEBUG_RA || defined DEBUG_BCPA
   fprintf(file, "DEBUG=1\n");
#else
   fprintf(file, "DEBUG=0\n");
#endif
   
   fprintf(file, "Parameters settings file=%s\n", param.parameter_stamp);
   fprintf(file, "Instance file=%s\n", instance_filename);
   ct = localtime(&t);
   fprintf(file, "Date=%d-%.2d-%.2d\nTime=%.2d:%.2d:%.2d\n", ct->tm_year+1900, ct->tm_mon, ct->tm_mday, ct->tm_hour, ct->tm_min, ct->tm_sec);
}

int readInstance(SCIP* scip, char* filename)
{
  SCIP_RETCODE retcode;
  SCIP_RESULT result;

  retcode = readerDssp(scip,filename, &result);

  switch( retcode )
  {
  case SCIP_NOFILE:
    SCIPinfoMessage(scip, NULL, "file <%s> not found\n", filename);
    return 0;
  case SCIP_PLUGINNOTFOUND:
    SCIPinfoMessage(scip, NULL, "no reader for input file <%s> available\n", filename);
    return 0;
  case SCIP_READERROR:
    SCIPinfoMessage(scip, NULL, "error reading file <%s>\n", filename);
    return 0;
  default:
    SCIP_CALL( retcode );
  }
  return result;
}
int contFrac(SCIP* scip)
{
  int i, totalFrac, nvars;
  SCIP_VAR **vars;
  SCIP_PROBDATA* probdata;
  double solval;
  
  assert(scip != NULL);
  probdata = SCIPgetProbData(scip);
  assert(probdata != NULL);

  nvars = SCIPprobdataGetNVars(probdata);
  vars = SCIPprobdataGetVars(probdata);

  totalFrac = 0;
  for(i = 0; i < nvars; i++) {
    solval = SCIPgetVarSol(scip, vars[i]);
    if(!SCIPisFeasIntegral(scip, solval)) {
      totalFrac++;
    }
  }
  return totalFrac;
}
SCIP_RETCODE printStatistic(SCIP* scip, char* outputname, int dc, int df)
{
   FILE* fout;
   char filename[SCIP_MAXSTRLEN];
   SCIP_Bool outputorigsol = FALSE;
   //   compilingOptionsT compilingOptions;
   SCIP_SOL* bestSolution;
   double tempo;
   SCIP_PROBDATA* probdata;
   int string_size, alpha_size, sc_size, sf_size, frac, nvars;
   
   bestSolution = SCIPgetBestSol(scip);
   /*******************
    * Solution Output *
    *******************/
   //   SCIP_CALL( SCIPgetBoolParam(scip, "misc/outputorigsol", &outputorigsol) );
   /* I found the name of those statistical information looking the scip source code at file scip.c (printPricerStatistics(), for instance)  */ 
   outputorigsol = 1;
   if ( outputorigsol )
   {
      SCIP_SOL* bestsol;

      if(param.display_freq>0){
        SCIPinfoMessage(scip, NULL, "\nprimal solution (original space):\n");
        SCIPinfoMessage(scip, NULL, "=================================\n\n");
      }
      
      bestsol = SCIPgetBestSol(scip);
      if(param.display_freq>0){
        if ( bestsol == NULL )
          SCIPinfoMessage(scip, NULL, "no solution available\n");
        else
        {
          SCIP_SOL* origsol;
          
          SCIP_CALL( SCIPcreateSolCopy(scip, &origsol, bestsol) );
          SCIP_CALL( SCIPretransformSol(scip, origsol) );
          SCIP_CALL( SCIPprintSol(scip, origsol, NULL, FALSE) );
          SCIP_CALL( SCIPfreeSol(scip, &origsol) );
        }
      }
   }
   else
   {
     if(param.display_freq>0){
       SCIPinfoMessage(scip, NULL, "\nprimal solution (transformed space):\n");
       SCIPinfoMessage(scip, NULL, "====================================\n\n");
       
       SCIP_CALL( SCIPprintBestSol(scip, NULL, FALSE) );
     }
   }


   /**************
    * Statistics *
    **************/
   if(param.display_freq>0){
     SCIPinfoMessage(scip, NULL, "\nStatistics\n");
     SCIPinfoMessage(scip, NULL, "==========\n\n");
     SCIP_CALL( SCIPprintStatistics(scip, NULL) );
   }
   /*******************************
    * Save run output information *
    *******************************/
   tempo = (agora-antes)/((double) CLOCKS_PER_SEC);

   /*****************************
    * Save best solution found  *
    *****************************/
   printSol(scip, outputname, tempo);
   
   (void) SCIPsnprintf(filename, SCIP_MAXSTRLEN, "%s.out", outputname);
   //   strcat(filename, ".out");

   fout = fopen(filename, "w");
   if(!fout){
     printf("\nProblem to create output file: %s", filename);
   }
   else{
     /* recupera dados do problema*/
     probdata = SCIPgetProbData(scip);
     assert(probdata != NULL);
     string_size = SCIPprobdataGetString_size(probdata);
     sf_size = SCIPprobdataGetSf_size(probdata);
     sc_size = SCIPprobdataGetSc_size(probdata);
     alpha_size = SCIPprobdataGetAlpha_size(probdata);
     nvars = SCIPprobdataGetNVars(probdata);

     /* contabiliza frac */
     frac = contFrac(scip);
     
     UB = SCIPgetPrimalbound(scip);
     LB = SCIPgetDualbound(scip);
     if(param.heur!=NONE){
       fprintf(fout,"%s;%d;%d;%d;%d;%d;%d;%lli;%lf;%lf;%lf;%lf;%lf;%lli;%d;%lf;%lf;%lli;%d;%s;%lf;%lf;%lf;%d;%d", instance_filename,alpha_size, sc_size,sf_size,string_size,frac,nvars,SCIPgetNRootLPIterations(scip), tempo, LB, UB, SCIPgetGap(scip), SCIPgetDualboundRoot(scip),SCIPgetNTotalNodes(scip), SCIPgetNNodesLeft(scip), SCIPgetSolvingTime(scip),SCIPgetTotalTime(scip),SCIPgetMemUsed(scip),SCIPgetStatus(scip), heurName[param.heur], heurValue, tempoHeur, firstUB, dc, df);
       printf("%s;%d;%d;%d;%d;%d;%d;%lli;%lf;%lf;%lf;%lf;%lf;%lli;%d;%lf;%lf;%lli;%d;%s;%lf;%lf;%lf;%d;%d", instance_filename,alpha_size, sc_size,sf_size,string_size,frac,nvars,SCIPgetNRootLPIterations(scip), tempo, LB, UB, SCIPgetGap(scip), SCIPgetDualboundRoot(scip),SCIPgetNTotalNodes(scip), SCIPgetNNodesLeft(scip), SCIPgetSolvingTime(scip),SCIPgetTotalTime(scip),SCIPgetMemUsed(scip),SCIPgetStatus(scip), heurName[param.heur], heurValue, tempoHeur, firstUB, dc, df);
     }
     else{
       fprintf(fout,"%s;%d;%d;%d;%d;%d;%d;%lli;%lf;%lf;%lf;%lf;%lf;%lli;%d;%lf;%lf;%lli;%d;%s;-;-;-;-;-", instance_filename,alpha_size, sc_size,sf_size,string_size,frac,nvars,SCIPgetNRootLPIterations(scip), tempo, LB, UB, SCIPgetGap(scip), SCIPgetDualboundRoot(scip),SCIPgetNTotalNodes(scip), SCIPgetNNodesLeft(scip), SCIPgetSolvingTime(scip),SCIPgetTotalTime(scip),SCIPgetMemUsed(scip),SCIPgetStatus(scip), heurName[param.heur]);
       printf("%s;%d;%d;%d;%d;%d;%d;%lli;%lf;%lf;%lf;%lf;%lf;%lli;%d;%lf;%lf;%lli;%d;%s;-;-;-;-;-", instance_filename,alpha_size, sc_size,sf_size,string_size,frac,nvars,SCIPgetNRootLPIterations(scip), tempo, LB, UB, SCIPgetGap(scip), SCIPgetDualboundRoot(scip),SCIPgetNTotalNodes(scip), SCIPgetNNodesLeft(scip), SCIPgetSolvingTime(scip),SCIPgetTotalTime(scip),SCIPgetMemUsed(scip),SCIPgetStatus(scip), heurName[param.heur]);
     }
     /*     if(bestSolution!=NULL){
       fprintf(fout, ";bestsol in %lld;%lf;%d;%s", SCIPsolGetNodenum(bestSolution),
         SCIPsolGetTime(bestSolution),
         SCIPsolGetDepth(bestSolution),
         SCIPsolGetHeur(bestSolution) != NULL
         ? SCIPheurGetName(SCIPsolGetHeur(bestSolution))
         : (SCIPsolGetRunnum(bestSolution) == 0 ? "initial" : "relaxation"));
         }
     */

     fprintf(fout, ";%s\n", param.parameter_stamp);
     printf(";%s\n", param.parameter_stamp);
     fclose(fout);
   }
   return SCIP_OKAY;
}

/** 
 * creates a SCIP instance with default plugins, and set SCIP parameters 
 */
SCIP_RETCODE configScip(
  SCIP** pscip
   )
{
  SCIP* scip = NULL;

  /*********
    * Setup *
    *********/

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );
   
   /* include crtp reader */
   //   SCIP_CALL( SCIPincludeReaderCrtp(scip) );
  
   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   if(!param.scipdefault){
     /* for column generation instances, disable restarts */
     SCIP_CALL( SCIPsetIntParam(scip,"presolving/maxrestarts",0) ); //param.presolve_maxrestarts) ); // 0

     /* turn off all separation algorithms */
     SCIP_CALL( SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE) ); // turn off

     /* disable heuristics */
     SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, TRUE) );  // turn off

// #if defined(BCPA) //|| defined(VNS)
//     SCIP_CALL( SCIPincludeHeurmyheuristic(scip) );
// #endif

     /* disable presolving */
     SCIP_CALL( SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE) );

     SCIP_CALL( SCIPsetIntParam(scip, "propagating/maxrounds", 0) );
     SCIP_CALL( SCIPsetIntParam(scip, "propagating/maxroundsroot", 0) );

     SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );
     SCIP_CALL( SCIPsetIntParam(scip, "separating/maxrounds", 0) );

     // Branching priority using pseudocost rule! Perhaps, remove this should be better!
     //   if(param.pscost_priority>0){
     //      SCIP_CALL( SCIPsetIntParam(scip, "branching/pscost/priority", param.pscost_priority) ); //1000000
     //   }
   }
   
   SCIP_CALL( SCIPsetIntParam(scip, "display/freq", param.display_freq) ); //50
  
   /* set time limit */
   SCIP_CALL( SCIPsetRealParam(scip, "limits/time", param.time_limit) ); // 1800

   /* only root node */
   //#ifdef ONLYROOT   
   //   param.nodes_limit = 1;
   //#endif
   SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", param.nodes_limit) );

   SCIP_CALL( SCIPsetRealParam(scip, "numerics/lpfeastol", 1e-6) ); // primal feasibility tolerance
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/dualfeastol", 1e-6) ); // dual feasibility tolerance

   *pscip = scip;
   return SCIP_OKAY;
}
/**
 * set default+user parameters
 **/
int setParameters(int argc, char** argv, parametersT* pparam)
{
  typedef struct{
    const char* description;
    const char* param_name;
    void* param_var;
    enum {INT, DOUBLE, STRING} type;
    int ilb;
    int iub;
    double dlb;
    double dub;
    int idefault;
    double ddefault;    
  } settingsT;

  enum {time_limit,display_freq,nodes_limit,presolve_maxrounds,presolve_maxrestarts,pscost_priority,propag,propag_maxrounds,propag_maxroundsroot, param_stamp, param_output_path, substring,pheur, local_search, vnstime, ngsize, vnssecond, semy, bcpa_rounding, ra_novo, ra_reduce, ra_merge, modelo_novo, problem, shaking, bcpatime_limit, scipdefault, total_parameters};

  settingsT parameters[]={
            {"time limit", "--time", &(param.time_limit), INT, -1, MAXINT, 0,0,600,0},
            {"display freq", "--display", &(param.display_freq), INT, -1, MAXINT, 0,0,1,0},
            {"nodes limit", "--nodes", &(param.nodes_limit), INT, -1, MAXINT, 0,0,1,0},
            {"presolve maxrounds", "--presolve-rounds", &(param.presolve_maxrounds), INT, -1, MAXINT, 0,0,0,0},
            {"presolve maxrestarts", "--presolve-restarts", &(param.presolve_maxrestarts), INT, -1, MAXINT, 0,0,0,0},
            {"pscost priority", "--pscost-prior", &(param.pscost_priority), INT, -1, MAXINT, 0,0,1000000,0},
            {"domain propagation", "--propag", &(param.propag), INT, 0, 1, 0,0,0,0},
            {"propag maxrounds", "--propag-rounds", &(param.propag_maxrounds), INT, -1, MAXINT, 0,0,0,0},
            {"propag maxrounds root", "--propag-rootrounds", &(param.propag_maxroundsroot), INT, -1, MAXINT, 0,0,0,0},
            {"param stamp", "--param_stamp", &(param.parameter_stamp), STRING, 0,0,0,0,0,0},
            {"output path", "--output_path", &(output_path), STRING, 0,0,0,0,0,0},
            {"substring?", "--substring", &(param.substring), INT, 0, 1, 0, 0, 0, 0},
            {"heur (0=none 1=ra 2=bcpa 3=blpl 4=sbpl)", "--heur", &(param.heur), INT, 0, 5, 0, 0, 0, 0},
            {"local search (0=none, 1=first, 2=best, 3=lp)", "--ls", &(param.local_search), INT, 0, 4, 0, 0, 0, 0},
            {"vns time", "--vnstime", &(param.vnstime), INT, 0, 300, 0, 0, 5, 0},
            {"ng size", "--ngsize", &(param.ngsize), INT, 0, 5, 0, 0, 3, 0},
            {"vns second", "--vnssec", &(param.vnssecond), INT, 0, 1, 0, 0, 0, 0},
            {"sem y?", "--semy", &(param.semy), INT, 0, 1, 0, 0, 1, 0},
            {"bcpa rounding?", "--bcpa_rounding", &(param.bcpa_rounding), INT, 0, 1, 0, 0, 1, 0},
            {"ra novo?", "--ra_novo", &(param.ra_novo), INT, 0, 1, 0, 0, 0, 0},
            {"ra reduce?", "--ra_reduce", &(param.ra_reduce), INT, 0, 1, 0, 0, 0, 0},
            {"ra merge?", "--ra_merge", &(param.ra_merge), INT, 0, 1, 0, 0, 0, 0},
            {"model_novo?", "--model_novo", &(param.model_novo), INT, 0, 1, 0, 0, 0, 0},
            {"problem (0=CSP, 1=FSP, 2=DSP)", "--problem", &(param.problem), INT, 0, 2, 0, 0, 2, 0},
            {"shaking (0=none, 1=normal, 2=lp)", "--shaking", &(param.shaking), INT, 0, 2, 0, 0, 0, 0},
            {"bcpa time limit", "--bcpa_time", &(param.bcpatime_limit), INT, -1, MAXINT, 0,0,300,0},
            {"scip default", "--scipdefault", &(param.scipdefault), INT, 0, 1, 0,0,0,0}
  };
  int i, j, ivalue, error;
  double dvalue;
  FILE *fin;

  // total_parameters = sizeof(parameters)/sizeof(parameters[0]);
  
  
  if (pparam==NULL)
    return 0;
  
  // check arguments
  if(argc<2){
    printf("\nSintaxe: program <instance-file> <parameters-setting>. \nUse program --options to show options.\n");
    return 0;
  }

  // show options
  if(!strcmp(argv[1],"--options")){
    //    showOptions();
    printf("\noption                : default -   range        : description");
    for(i=0;i<total_parameters;i++){
      switch(parameters[i].type){
      case INT:
        printf("\n%-22s: %7d - [%3d,%8d] : %s", parameters[i].param_name, parameters[i].idefault, parameters[i].ilb, parameters[i].iub,parameters[i].description);
        break;
      case DOUBLE:
        printf("\n%-22s: %7.1lf - [%3.1lf,%8.1lf] : %s", parameters[i].param_name, parameters[i].ddefault, parameters[i].dlb, parameters[i].dub,parameters[i].description);
        break;
      case STRING:
        printf("\n%-22s:       * - [*,*]          : %s", parameters[i].param_name, parameters[i].description);
      }
    }
    printf("\n");
    return 0;
  }

  // check the existance of instance file
  fin = fopen(argv[1], "r");
  if(!fin){
      printf("\nInstance file not found: %s\n", argv[1]);
      return 0;
  }
  fclose(fin);  

  // set default parameters value
  for(i=0;i<total_parameters;i++){
    if(parameters[i].type==INT)
      *((int*)(parameters[i].param_var)) = parameters[i].idefault;
    else if (parameters[i].type==DOUBLE)
      *((double*)(parameters[i].param_var)) = parameters[i].ddefault;
    else
      *((char**) (parameters[i].param_var)) = NULL;
  }
  output_path = current_path;

  // set user parameters value
  error = 0;
  for(i=2;i<argc && !error;i+=2){
    for(j=0;j<total_parameters && strcmp(argv[i],parameters[j].param_name);j++)
      ;
    if(j>=total_parameters || i==argc-1){
      printf("\nParameter (%s) invalid or uncompleted.", argv[i]);
      error = 1;
    }
    else{
      switch(parameters[j].type){
      case INT:
        ivalue = atoi(argv[i+1]);
        if(ivalue < parameters[j].ilb || ivalue > parameters[j].iub){
          printf("\nParameter (%s) value (%d) out of range [%d,%d].", argv[i], ivalue, parameters[j].ilb, parameters[j].iub);
          error = 1;
          //          break;
        }
        else{
          *((int*)(parameters[j].param_var)) = ivalue;          
        }
        break;
      case DOUBLE:
        dvalue = atof(argv[i+1]);
        if(dvalue < parameters[j].dlb || dvalue > parameters[j].dub){
          printf("\nParameter (%s) value (%lf) out of range [%lf,%lf].", argv[i], dvalue, parameters[j].dlb, parameters[j].dub);
          error = 1;
          //          break;
        }
        else{
          *((double*)(parameters[j].param_var)) = dvalue;          
        }
        break;
      case STRING:
        *((char**) (parameters[j].param_var)) = argv[i+1];
      }
    }
  }

  // print parameters
  printf("\n\n----------------------------\nParameters settings");
  for(i=0;i<total_parameters;i++){
    printf("\nparameter %s (%s): default value=", parameters[i].param_name, parameters[i].description);
    switch(parameters[i].type){
    case INT:
      printf("%d - value = %d", parameters[i].idefault, *( (int*)parameters[i].param_var));
      break;
    case DOUBLE:
      printf("%lf - value = %lf", parameters[i].ddefault, *( (double*)parameters[i].param_var));
      break;
    case STRING:
      printf("(null) - value = %s", *( (char**)parameters[i].param_var));
      break;
    }
  }
  printf("\nerror = %d\n", error);
  if(!error){
    FILE* fout;
    char foutname[SCIP_MAXSTRLEN];

    if(param.parameter_stamp != NULL){
      // complete the fullname of the parameters stamp file
      sprintf(foutname, "%s/%s", output_path, param.parameter_stamp);
    }
    else{
      // define the parameters stamp file using stamp default = date-time
      struct tm * ct;
      const time_t t = time(NULL);
      param.parameter_stamp = (char*) malloc(sizeof(char)*100);
      ct = localtime(&t);
      snprintf(param.parameter_stamp, 100, "d%d%.2d%.2dh%.2d%.2d%.2d", ct->tm_year+1900, ct->tm_mon, ct->tm_mday, ct->tm_hour, ct->tm_min, ct->tm_sec);

      sprintf(foutname, "%s/%s", output_path, param.parameter_stamp);
    }
    // check if the stamp already exists
    fout = fopen(foutname, "r");
    // TODO: opendir() should be done first to avoid open a directory!
    if(!fout){
      // save parameters in the stamp file
      printf("\nwriting parameters in %s", foutname);
      fout = fopen(foutname, "w+");
      for(i=0;i<total_parameters;i++){
        fprintf(fout, "%s ", parameters[i].param_name);
        switch(parameters[i].type){
        case INT:
          fprintf(fout,"%d\n", *( (int*)parameters[i].param_var));
          break;
        case DOUBLE:
          fprintf(fout,"%lf\n", *( (double*)parameters[i].param_var));
          break;
        case STRING:
          fprintf(fout,"%s\n", *( (char**)parameters[i].param_var));
          break;
        }
      }
      fclose(fout);
    }
    else{
      // check if the stamp is valid
      char param_name[100];
      char svalue[100];
      while(!feof(fout)){
        fscanf(fout,"%s", param_name);
        for(j=0;j<total_parameters && strcmp(param_name,parameters[j].param_name);j++)
          ;
        if(j>=total_parameters){
          printf("\nParameter (%s) invalid or uncompleted.", param_name);
          error = 1;
          break;
        }
        else{
          switch(parameters[j].type){
          case INT:
            fscanf(fout, "%d\n", &ivalue);
            if(ivalue < parameters[j].ilb || ivalue > parameters[j].iub){
              printf("\nParameter (%s) value (%d) out of range [%d,%d].", param_name, ivalue, parameters[j].ilb, parameters[j].iub);
              error = 1;
            }
            else if(ivalue != *((int*)(parameters[j].param_var))){
              printf("\nParameter (%s) value (%d) differs to saved value = %d.", param_name, ivalue, *((int*)(parameters[j].param_var)));
              error = 1;
            }
            break;
          case DOUBLE:
            fscanf(fout, "%lf\n", &dvalue);
            if(dvalue < parameters[j].dlb || dvalue > parameters[j].dub){
              printf("\nParameter (%s) value (%lf) out of range [%lf,%lf].", param_name, dvalue, parameters[j].dlb, parameters[j].dub);
              error = 1;
            }
            else if(fabs(dvalue - *((double*)(parameters[j].param_var))) > EPSILON){
              printf("\nParameter (%s) value (%lf) differs to saved value = %lf.", param_name, dvalue, *((double*)(parameters[j].param_var)));
              error = 1;
            }
            break;
          case STRING:
            fscanf(fout, "%s\n", svalue);
            if(strcmp(svalue,*((char**)(parameters[j].param_var)))){
              printf("\nParameter (%s) value (%s) differs to saved value = %s.", param_name, svalue, *((char**)(parameters[j].param_var)));
              error = 1;
            }
            break;
          }
        } // each parameter
      } // while !feof
      fclose(fout);
    }
  }
  return !error;
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
    double tempo, LB, UB, LB_ROOT, ra = 1000000;//, vns = SCIP_INVALID;
    clock_t antes, agora;
    //SCIP_VAR **vars, *var;
    SCIP* scip = NULL;
    SCIP_PROBDATA *probdata;
    int sf_size;
    FILE *file; 
    double vns=1000000, bcpa = 1000000, vns0 = 1000000;
    SCIP_SOL *bestSolution;
    SCIP_Real solval;   
    FILE *fp;
    int success;


    antes = clock();

    /*********
    * Setup *
    *********/
    
    /* initialize SCIP */
    SCIP_CALL( SCIPcreate(&scip) );

    /* include dssp reader */
    SCIP_CALL( SCIPincludeReaderDSSP(scip) );
    //    SCIP_CALL( SCIPincludeReadercsp(scip) );


    /* include default SCIP plugins */
    SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
    
     /* include Xyz constraint handler */
     //SCIP_CALL( SCIPincludeConshdlrXyz(scip) );   


    /* for column generation instances, disable restarts */
    SCIP_CALL( SCIPsetIntParam(scip,"presolving/maxrestarts",0) );

    /* turn off all separation algorithms */
    if(!param.scipdefault){
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
    }
    /* set time limit */
    SCIP_CALL( SCIPsetRealParam(scip, "limits/time", 1800));//TIME_LIMIT) );

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
    strcpy(solFileName, "");
    char date_time_aux[256], date_aux[256];
    time_t t_aux = time(NULL);
    struct tm tm = *localtime(&t_aux);
    sprintf(date_aux, "%d-%02d", tm.tm_year + 1900, tm.tm_mon + 1);
    sprintf(date_time_aux, "-%d-%02d-%02d_%02d-%02d-%02d", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    
    int check = mkdir("results", 0777);

    strcat(solFileName, "results/");
    strcat(solFileName, date_aux);

    check = mkdir(solFileName, 0777); 

    strcat(solFileName, "/");
    char dir_grupo1[512], dir_grupo2[512], dir_grupo3[512];
    sprintf(dir_grupo1, "%sdata/", solFileName);
    sprintf(dir_grupo2, "%sdata/", solFileName);
    sprintf(dir_grupo3, "%sdata/", solFileName);
    check = mkdir(dir_grupo1, 0777);
    check = mkdir(dir_grupo2, 0777);
    check = mkdir(dir_grupo3, 0777);

    //    strcat(dir_grupo1, "instancias-tcc/");
    //    strcat(dir_grupo2, "instancias-tcc/");
    //    strcat(dir_grupo3, "instancias-tcc/");
    //    check = mkdir(dir_grupo1, 0777);
    //    check = mkdir(dir_grupo2, 0777);
    //    check = mkdir(dir_grupo3, 0777);

    strcat(dir_grupo1, "grupo1/");
    strcat(dir_grupo2, "grupo2/");
    strcat(dir_grupo3, "grupo3/");
    check = mkdir(dir_grupo1, 0777);
    check = mkdir(dir_grupo2, 0777);
    check = mkdir(dir_grupo3, 0777);

    strcat(solFileName, argv[2]);
    
    /* grava melhor solucao */
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
#ifdef VBPL
    strcat(solFileName, "-VBPL");
#elif BLPL
    strcat(solFileName, "-BLPL");
#else
    strcat(solFileName, "-PBPL");
#endif

#ifdef TIMESSS
    strcat(solFileName, "-5S");
#elif TIMESS
    strcat(solFileName, "-30S");
#elif TIMES
    strcat(solFileName, "-60S");
#elif TIMEM
    strcat(solFileName, "-120S");
#elif TIMEL
    strcat(solFileName, "-180S");
#elif TIMELL
    strcat(solFileName, "-300S");
#else
    strcat(solFileName, "-3000S");
#endif

#ifdef NSIZES
    strcat(solFileName, "-SIZE2");
#elif NSIZEM
    strcat(solFileName, "-SIZE3");
#elif NSIZEL
    strcat(solFileName, "-SIZE4");
#else
    strcat(solFileName, "-SIZE4");
#endif

#ifdef BEST
    strcat(solFileName, "-BEST");
#elif FIRST
    strcat(solFileName, "-FIRST");
#endif


#ifdef OMEGAS
    strcat(solFileName, "-OMEGAS");
#elif OMEGAM
    strcat(solFileName, "-OMEGAM");
#else
    strcat(solFileName, "-OMEGAL");
#endif

    strcat(solFileName, date_time_aux);
    //printf("SOLFILENAME: %s\n", solFileName);

    //    printSol(scip, solFileName, tempo);
    char solFileNameResumo[512];
    strcpy(solFileNameResumo, solFileName);
    strcat(solFileNameResumo, ".resumo");


    file = fopen(solFileNameResumo, "w");
    if(!file){
      printf("\nProblem to create the file %s", solFileNameResumo);
      exit(1);
    }
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
    success = heur_RA(scip, solFileName, tempo, &ra);
    char name[SCIP_MAXSTRLEN];
    agora = clock();
    tempo2 = (agora - antes) / ((double) CLOCKS_PER_SEC);
    FILE *file2;
    (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "RA.csv");
        
    file2 = fopen(name, "a");
    if(success){
      fprintf(file2, "%s;%lf;%lf;\n", argv[2], tempo+tempo2, ra);
    }
    else{
      fprintf(file2, "%s;%lf;-;\n", argv[2], tempo+tempo2);
    }
    fclose(file2);
    #endif
#endif
    // int bcpa;

#ifdef BCPA
    #ifdef ONLYROOT
    antes = clock();
    success = heur_BCPA(scip, argv[2], &bcpa);
    char name[SCIP_MAXSTRLEN];
    agora = clock();
    tempo2 = (agora - antes) / ((double) CLOCKS_PER_SEC);
    FILE *file2;
    (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "BCPA.csv");
        
    file2 = fopen(name, "a");
    if(success){
      fprintf(file2, "%s;%lf;%lf;\n", argv[2], tempo+tempo2, bcpa);
    }
    else{
      fprintf(file2, "%s;%lf;-;\n", argv[2], tempo+tempo2);
    }
    fclose(file2);
    #endif
#endif

#ifdef VNS
    #ifdef ONLYROOT
    antes = clock();    
    success = vns_main(scip, argv[2], tempo, &vns, &vns0);

    agora = clock();
    tempo2 = (agora - antes) / ((double) CLOCKS_PER_SEC);
    char name[SCIP_MAXSTRLEN];

    FILE *file2;
    (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "VNS.csv");
        
    file2 = fopen(name, "a");
    if(success){
      fprintf(file2, "%s;%lf;%lf;%lf\n", argv[2], tempo+tempo2, vns, vns0);//bcpa);
    }
    else{
      fprintf(file2, "%s;%lf;-;-\n", argv[2], tempo+tempo2);//bcpa);
    }
    fclose(file2);
    #endif
#endif

#ifdef ONLYROOT
    printRelaxedSol(scip, solFileName, tempo);
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
    fprintf(fp, ";%ld;%lf;%lf;%lf;%lf;%d;%lf;%lf;VNS\n", nNodes, LB, UB, LB_ROOT, tempo,SCIPgetStatus(scip), vns, vns0);
    
    fclose(fp);

    /*    if(UB>vns+EPSILON){
      UB = vns;
      }*/
    printf(";%lf;%lf;%lf;VNS\n", tempo2, vns, vns0);
    fprintf(file, ";%lf;%lf;%lf;%lf;%d;%lf;%lf;VNS\n", LB, UB, LB_ROOT, tempo,SCIPgetStatus(scip), vns, vns0);    
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

void printHeurSol(SCIP* scip, char* filename, int* x, int dc, int df, double tempo)
{
  FILE* file;
  char name[SCIP_MAXSTRLEN], *alphabet;
  int i, j, string_size, alpha_size;
  SCIP_PROBDATA* probdata;
  
  probdata = SCIPgetProbData(scip);
  assert(probdata != NULL);

  string_size = SCIPprobdataGetString_size(probdata);
  alphabet = SCIPprobdataGetAlphabet(probdata);
  alpha_size = SCIPprobdataGetAlpha_size(probdata);
  
  //  (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s-%s.sol", filename, heurName[param.heur]);
  (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s.sol", filename);

  file = fopen(name, "w");
    
  if(!file)
  {
    printf("\nProblem to create solution file: %s\n", name);
    return;
  }
  
  fprintf(file, "Time: %lfs\n", tempo);
  fprintf(file, "Objective Value(dc - df) = %d\n", dc-df);
  fprintf(file, "dc = %d\n", dc);
  fprintf(file, "df = %d\n", df);
  fprintf(file, "Target String(x) = ");
    
  for(i = 0; i < string_size; i++) {
    fprintf(file, "%c", symbol_value(alphabet, alpha_size, x[i]));
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
      if (x[i] != j) {
        fprintf(file, "0 ");
      } else {
        fprintf(file, "1 ");
      }
    }
    fprintf(file, "\n");        
  }
  printCompilerSettings(file);
  fclose(file);
}

/**
 * main ()
 */
int main(int argc, char** argv)
{
   SCIP* scip = NULL;
   int success, *heur, dc, df;
   char outputname[SCIP_MAXSTRLEN];

   // set default+user parameters
   if(!setParameters(argc, argv, &param))
     return 0;

   // check invalid parameters settings
   if(param.heur==SBPL && (param.shaking != LP_SHAKING || param.local_search != LP_IMPROV)){
     printf("\nConfiguracao invalida: SBPL incompativel com normal shaking and first/best improvement\n");
     return 0;
   }
   if(param.heur==BLPL && (param.shaking != NORMAL_SHAKING || param.local_search != BEST_IMPROV)){
     printf("\nConfiguracao invalida: BLPL incompativel com LP shaking and first/LP improvement\n");
     return 0;
   }
   if(param.heur==VBPL && (param.shaking != NORMAL_SHAKING || param.local_search != FIRST_IMPROV)){
     printf("\nConfiguracao invalida: SBPL incompativel com LP shaking and best/LP improvement\n");
     return 0;
   }

   // initialize global vars

   // start clock
   antes = clock();

   // create and config SCIP using default+user parameters 
   configScip(&scip);

   // split fullfilename in path and filename
   splitfullfilename(argv[1]);
   
   // read instance
   if(!readInstance(scip, argv[1])){
     return 1;
   }

   // initial settings
   heurValue = MAXINT;
   firstUB = MAXINT;
   if(param.heur!=NONE){     
     heur = (int*) malloc(sizeof(int)*SCIPprobdataGetString_size(SCIPgetProbData(scip)));
   }

   /*******************
    * Problem Solving *
    *******************/
   /* solve problem */
   SCIPinfoMessage(scip, NULL, "\nsolve problem\n");
   SCIPinfoMessage(scip, NULL, "=============\n\n");

   SCIP_CALL( SCIPsolve(scip) );
   agora = clock();

   configOutputName(outputname, argv[1], argv[0]);

   // call heuristics
   printf("\nSolving heuristic: %s...\n", heurName[param.heur]);
   tempoHeur = clock();
   switch(param.heur){
   case RA:
     success = heur_RA(scip, outputname, (agora-antes)/CLOCKS_PER_SEC, &heurValue, heur, &dc, &df);
     break;
   case BCPA:
     success = heur_BCPA(scip, outputname, &heurValue, heur, &dc, &df);
     break;
   case VBPL:
   case BLPL:
   case SBPL:
     success = vns_main(scip, outputname, (agora-antes)/CLOCKS_PER_SEC, &heurValue, &firstUB, heur, &dc, &df);
     break;
   case NONE:
     success = 1;
   }
   tempoHeur = (clock() - tempoHeur) / ((double) CLOCKS_PER_SEC); 

   if(!success){
     printf("\nError...\n");
     return 1;
   }

   // save solution found by heuristics
   if(param.heur!=NONE){
     printHeurSol(scip, outputname, heur, dc, df, tempoHeur);
   }
   // print statistics
   printStatistic(scip, outputname, dc, df);

   // print cost and time for irace
   //   if(SCIPgetGap(scip)>0){
   //     printf("\n%lf %.1lf", 100*SCIPgetGap(scip), (double) param.time_limit);
   //   }
   //   else{
   //     printf("\n%lf %.1lf",0.0,  SCIPgetTotalTime(scip));
   //   }
   ////printf("\n%lf %.1lf", 100*SCIPgetGap(scip)+SCIPgetTotalTime(scip)/param.time_limit, SCIPgetTotalTime(scip));

   /********************
    * Deinitialization *
    ********************/
   //   SCIP_CALL( SCIPfree(&scip) );
   if(param.heur!=NONE){     
     free(heur);
   }
   BMScheckEmptyMemory();

   return 0;
}

