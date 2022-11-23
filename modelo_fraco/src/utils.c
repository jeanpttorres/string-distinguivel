#include "utils.h"
#include "parameters_dssp.h"

/* put your local methods here, and declare them static */
int randomInteger (int low, int high)
{
  int k;
  double d;

  d = (double) rand () / ((double) RAND_MAX + 1);
  k = d * (high - low + 1);
  return low + k;
}

int objective_value(int  *sum_sc, int *sum_sf, char **Sc, char **Sf, int sc_size, int sf_size, int kc, int kf, char *alpha, int alpha_size, int string_size, int *solution) {
    int dif = 0;
    int dc = 99999;
    int df = -1;
    int aux = 0;
    int menor_dif_sc, maior_dif_sf;
    int i, j, k, l;
    
    dc = -1 * string_size;

    if (sc_size > 0) {
        for (i = 0; i < sc_size; i ++) {
          menor_dif_sc = (sum_sc[i+1]+string_size-1) - sum_sc[i];
            
            for (j = 0; j <= (sum_sc[i+1]+string_size-1) - sum_sc[i] -string_size ; j ++) {
                l = 0;
                dif = 0;
                
                for (k = j; k < j + string_size; k ++) {
                    if (symbol_value(alpha, alpha_size, solution[l++]) != Sc[i][k]) {
                    dif ++;
                    }
                }
                menor_dif_sc = menor_dif_sc < dif ? menor_dif_sc : dif;
                
            }
            dc = menor_dif_sc > dc ? menor_dif_sc : dc; 
        }
    }
    else dc = 0;
    
    //Checar se obedece em sf
    df = string_size;
    if (sf_size > 0) {
        for (i = 0; i < sf_size; i ++) {
            maior_dif_sf = -1;
            for (j = 0; j <= (sum_sf[i+1]+string_size-1) - sum_sf[i] - string_size ; j ++) { l = 0;
                dif = 0;
                for (k = j; k < j + string_size; k ++) {
                    if (symbol_value(alpha, alpha_size, solution[l++]) != Sf[i][k]) {
                    dif ++;
                    }
                }        
                // maior_dif_sf = dif > maior_dif_sf ? dif : maior_dif_sf;
                df = dif < df ? dif : df; 
            }

        }
    }
    else df = 0;

    //    printf("\ndc=%d df=%d\n", dc, df);
    return dc - df;
}

int dc_df_objective_value(int  *sum_sc, int *sum_sf, char **Sc, char **Sf, int sc_size, int sf_size, int kc, int kf, char *alpha, int alpha_size, int string_size, int *solution, int* dc, int* df) {
    int dif = 0;
    int aux = 0;
    int menor_dif_sc, maior_dif_sf;
    int i, j, k, l;
    
    //int dc = 99999;
    //int df = -1;

    *dc = -1 * string_size;
    if (sc_size > 0) {

        for (i = 0; i < sc_size; i ++) {
            menor_dif_sc = (sum_sc[i+1]+string_size-1) - sum_sc[i];
            
            for (j = 0; j <= (sum_sc[i+1]+string_size-1) - sum_sc[i] -string_size ; j ++) {
                l = 0;
                dif = 0;
                
                for (k = j; k < j + string_size; k ++) {
                    if (symbol_value(alpha, alpha_size, solution[l++]) != Sc[i][k]) {
                    dif ++;
                    }
                }
                
                menor_dif_sc = menor_dif_sc < dif ? menor_dif_sc : dif;
                
            }
            *dc = menor_dif_sc > *dc ? menor_dif_sc : *dc; 
        }
    }
    else *dc = 0;
    
    //Checar se obedece em sf
    *df = string_size;

    if (sf_size > 0) {
        for (i = 0; i < sf_size; i ++) {
            maior_dif_sf = -1;
            for (j = 0; j <= (sum_sf[i+1]+string_size-1) - sum_sf[i] - string_size ; j ++) {
                l = 0;
                dif = 0;
                for (k = j; k < j + string_size; k ++) {
                    if (symbol_value(alpha, alpha_size, solution[l++]) != Sf[i][k]) {
                    dif ++;
                    }
                }        
            //            maior_dif_sf = dif > maior_dif_sf ? dif : maior_dif_sf;
            *df = dif < *df ? dif : *df; 
            }
        //        *df = maior_dif_sf < *df ? maior_dif_sf : *df; 
        }
    }
    else *df = 0;
    

    return *dc - *df;
}

int valorSolucao(SCIP *scip, char *filename, int *heur, int *pdc, int *pdf){
   int i, j, k, l, dif, dc, df, menor_dif_sc, maior_dif_sf;
    int  *sum_sc, *sum_sf, string_size, alpha_size, sc_size, sf_size;
    char *alphabet, **sc, **sf;

    SCIP_PROBDATA *probdata;


    assert(scip != NULL);

    probdata = SCIPgetProbData(scip);
    assert(probdata != NULL);


    string_size = SCIPprobdataGetString_size(probdata);
    alphabet = SCIPprobdataGetAlphabet(probdata);
    alpha_size = SCIPprobdataGetAlpha_size(probdata);
    sc = SCIPprobdataGetSc(probdata);
    sf = SCIPprobdataGetSf(probdata);


    sum_sc = SCIPprobdataGetSum_sc(probdata);
    sum_sf = SCIPprobdataGetSum_sf(probdata);
    sc_size = SCIPprobdataGetSc_size(probdata);
    sf_size = SCIPprobdataGetSf_size(probdata);
 
    
    
    //Checar dc
    dif = 0;
    if (sc_size > 0) {
      dc = -string_size;
      for (i = 0; i < sc_size; i ++) {
        menor_dif_sc = (sum_sc[i+1]+string_size-1) - sum_sc[i];
        for (j = 0; j <= (sum_sc[i+1]+string_size-1) - sum_sc[i] -string_size ; j ++) {
	  l = 0;
	  dif = 0;
          
	  for (k = j; k < j + string_size; k ++) {
	    if (symbol_value(alphabet, alpha_size, heur[l++]) != sc[i][k]) {
	      dif ++;
	    }
            }
	  menor_dif_sc = menor_dif_sc < dif ? menor_dif_sc : dif;
          
        }
        dc = menor_dif_sc > dc ? menor_dif_sc : dc; 
      }
    }
    else
      dc = 0;
    
    //Checar se obedece em sf
    if (sf_size > 0) {
      df = string_size;
      for (i = 0; i < sf_size; i ++) {
        maior_dif_sf = -1;
        for (j = 0; j <= (sum_sf[i+1]+string_size-1) - sum_sf[i] - string_size ; j ++) {
	  l = 0;
	  dif = 0;
	  for (k = j; k < j + string_size; k ++) {
	    if (symbol_value(alphabet, alpha_size, heur[l++]) != sf[i][k]) {
	      dif ++;
	    }
	  }        
	  // maior_dif_sf = dif > maior_dif_sf ? dif : maior_dif_sf;
	  df = dif < df ? dif : df; 
        }
        
      }
    }
    else
      df = 0;

    *pdc = dc;
    *pdf = df;
    
    return dc - df;
}


void printRelaxedSol(SCIP* scip, char *filename, double tempo)
{
    int i, j, string_size, alpha_size, *sum_sc, nvars;
    char name[SCIP_MAXSTRLEN], *alphabet;
    FILE *file;
    
    SCIP_PROBDATA *probdata;
    SCIP_VAR **vars, *var;
    SCIP_Real solval, currentSolVal = 0.0;
    int totalFrac, desl, sc_size, *v, kcounter, k, l;
    
    assert(scip != NULL);
       
    probdata = SCIPgetProbData(scip);
    assert(probdata != NULL);
    
    nvars = SCIPprobdataGetNVars(probdata);
    vars = SCIPprobdataGetVars(probdata);
    string_size = SCIPprobdataGetString_size(probdata);
    alphabet = SCIPprobdataGetAlphabet(probdata);
    alpha_size = SCIPprobdataGetAlpha_size(probdata);
    sum_sc = SCIPprobdataGetSum_sc(probdata);
    sc_size = SCIPprobdataGetSc_size(probdata);
    
    strcpy(name, filename);
    //printf("%s\n", name);
    strcat(name, ".lp.sol");
    file = fopen(name, "w");
    
    if(!file)
    {
        //printf("\nProblem to create solution file: %s", name);
        return;
    }
    
    name[0] = '\0';
    
    fprintf(file, "Time: %lfs\n", tempo);
        fprintf(file, "Objective Value(dc - df) = %.0lf\n", SCIPgetPrimalbound(scip));
    
    fprintf(file, "Target String(x) = ");
    int posSolVal = 0;
    for(i = 0; i < string_size; i++) {
        desl=0;
        currentSolVal = 1.1;
        posSolVal = 0;
        for (j = 0; j < alpha_size; j++) {
            var = vars[i * alpha_size + j];
	    solval = SCIPgetVarSol(scip,var);
            if (solval < currentSolVal) posSolVal = j;
            if(SCIPisFeasIntegral(scip, solval) && solval < 0.00001)
            {
                fprintf(file, "%c", symbol_value(alphabet, alpha_size, j));
		        desl=1;
            }
        }
	if(desl==0){
        //fprintf(file, "%c", symbol_value(alphabet, alpha_size, posSolVal));

	    fprintf(file, "@");
	}
	
    }
   
    fprintf(file, "\n");
    
    /* dc */
    var = vars[nvars - 2];
    solval = SCIPgetVarSol(scip,var);
    fprintf(file, "dc = %.0lf\n", solval);
    
    var = vars[nvars - 1];
    solval = SCIPgetVarSol(scip,var);
    fprintf(file, "df = %.0lf\n", solval);

    v=(int*)malloc(sizeof(int)*nvars);
    for(i=0;i<nvars;i++){
      v[i]=0;
    }
    totalFrac=0;    

    if(param.model_novo){
       //#ifdef NOVO
    fprintf(file, "Var fracionarias = ");
    for (i = 0, kcounter = 0; i < sc_size; i ++) { //printf("F1\n"); fflush(stdout);
	 for (j = 0; j < string_size; j ++) { //printf("   F2\n"); fflush(stdout);
	    for (k = 0; k < sum_sc[i+1] - sum_sc[i]; k++) {  //printf("      F3\n"); fflush(stdout);
	      for (l = 0; l < alpha_size; l++) { //printf("         F4\n"); fflush(stdout);
                var = vars[kcounter++];
                solval = SCIPgetVarSol(scip, var); 
		if(!SCIPisFeasIntegral(scip, solval)){
		  fprintf(file, "\nx_%d_%d_%c_%d=%lf", i,j+k,symbol_value(alphabet, alpha_size, l),j, solval);
		}
	      }	
            }                    
	 }
    }

    fprintf(file, "Var inteiras 1.0 = ");
    for (i = 0, kcounter = 0; i < sc_size; i ++) { //printf("F1\n"); fflush(stdout);
	 for (j = 0; j < string_size; j ++) { //printf("   F2\n"); fflush(stdout);
	    for (k = 0; k < sum_sc[i+1] - sum_sc[i]; k++) {  //printf("      F3\n"); fflush(stdout);
	      for (l = 0; l < alpha_size; l++) { //printf("         F4\n"); fflush(stdout);
                var = vars[kcounter++];
                solval = SCIPgetVarSol(scip, var); 
		if(SCIPisFeasIntegral(scip, solval) && solval >= 1.0 - 0.00001){
		  fprintf(file, "\nx_%d_%d_%c_%d=%lf", i,j+k,symbol_value(alphabet, alpha_size, l),j, solval);
		}
	      }	
            }                    
	 }
    }    
    }
    else{
    //#else  
    fprintf(file, "Var fracionarias = ");
    for(i = 0; i < string_size; i++) {
        for (j = 0; j < alpha_size; j++) {
            var = vars[i * alpha_size + j];
	    solval = SCIPgetVarSol(scip,var);
	    if(!SCIPisFeasIntegral(scip, solval))
	      {
		fprintf(file, "\nx_%d_%c=%lf", i,symbol_value(alphabet, alpha_size, j), solval);
	      }
        }
    }


    fprintf(file, "Var inteiras 0.0 = ");
    for(i = 0; i < string_size; i++) {
        for (j = 0; j < alpha_size; j++) {
            var = vars[i * alpha_size + j];
	    solval = SCIPgetVarSol(scip,var);
	    if(SCIPisFeasIntegral(scip, solval) && solval <= 0.00001)
	      {
		fprintf(file, "\nx_%d_%c=%lf", i,symbol_value(alphabet, alpha_size, j), solval);
	      }
        }
    }

    fprintf(file, "Var fracionarias = ");
    for(i = 0; i < sc_size; i++) {
        for (j = 0; j < sum_sc[i + 1] - sum_sc[i]; j++) {
            var = vars[VARS_X + sum_sc[i] + j];
	    solval = SCIPgetVarSol(scip,var);
	    if(!SCIPisFeasIntegral(scip, solval))
	      {
		fprintf(file, "\ny_%d_%d=%lf", i, j, solval);
	      }
        }
    }
    fprintf(file, "Var inteiras 1.0 = ");
    for(i = 0; i < sc_size; i++) {
        for (j = 0; j < sum_sc[i + 1] - sum_sc[i]; j++) {
            var = vars[VARS_X + sum_sc[i] + j];
	    solval = SCIPgetVarSol(scip,var);
	    if(SCIPisFeasIntegral(scip, solval) && solval >= 1.0 - 0.00001)
	      {
		fprintf(file, "\ny_%d_%d=%lf", i, j, solval);
	      }
        }
    }
    //#endif
    }
    

    totalFrac=0;int totalNFRAC=0;
    for(i=0;i<nvars; i++){
      solval = SCIPgetVarSol(scip,vars[i]);
      if(!SCIPisFeasIntegral(scip, solval))
	{
	  totalFrac++;
	}else totalNFRAC +=1;
    }

    fprintf(file, "\nTotal fracionarias=%d\n", totalFrac);
    fclose(file);

    printf(";%d;%d;", totalFrac, totalNFRAC);
    free(v);
}

void validSol(SCIP *scip, char *filename) {
    int i, j, k, l, sum;
    int **x, **y, dc, df, *sum_sc, *sum_sf, string_size, alpha_size, sc_size, sf_size;
    char name[SCIP_MAXSTRLEN], **sc, **sf, *alphabet;
    char sol_x[999], aux[2];
    bool valid = true;
    FILE *file;
    FILE *input;

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
#ifndef NOVO
    (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s-antigo.debug", filename);
#else
    (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s-antigo.debug", filename);
#endif
    file = fopen(name, "w");
    fprintf(file, "Validacao...\n");
    vars = SCIPprobdataGetVars(probdata);
    string_size = SCIPprobdataGetString_size(probdata);
    alphabet = SCIPprobdataGetAlphabet(probdata);
    alpha_size = SCIPprobdataGetAlpha_size(probdata);
    sc = SCIPprobdataGetSc(probdata);
    sf = SCIPprobdataGetSf(probdata);
    sum_sc = SCIPprobdataGetSum_sc(probdata);
    sum_sf = SCIPprobdataGetSum_sf(probdata);
    sc_size = SCIPprobdataGetSc_size(probdata);
    sf_size = SCIPprobdataGetSf_size(probdata);
    
    x = (int **) calloc(string_size, sizeof(int *));
    for (i = 0; i < string_size; i++) {
        x[i] = (int *) calloc(alpha_size, sizeof(int));
    }

    y = (int **) calloc(sc_size, sizeof(int *));
    for (i = 0; i < sc_size; i++) {
        y[i] = (int *) calloc(sum_sc[i], sizeof(int));
    }
    
    
    
    for(i = 0; i < string_size; i++) {
        sum = 0;
        for (j = 0; j < alpha_size; j++) {
            var = vars[i * alpha_size + j];
            solval = SCIPgetSolVal(scip, bestSolution, var);
            x[i][j] = (int) solval;
            sum += x[i][j];
        }
        if (sum != alpha_size - 1) {
            fprintf(file, "Problema na restricao (1) na posicao %d da string x.\n", i);
            valid = false;
        }
    }

    for (i = 0; i < sc_size; i++) {
        for (k = 0; k < sum_sc[i + 1] - sum_sc[i]; k++) {
            var = vars[VARS_X + sum_sc[i] + k];
            y[i][k] = (int) SCIPgetSolVal(scip, bestSolution, var);
        }
    }

    var = vars[VARS_X + VARS_Y];
    dc = (int) (SCIPgetSolVal(scip, bestSolution, var) + 0.5);

    var = vars[VARS_X + VARS_Y + 1];
    df = (int) (SCIPgetSolVal(scip, bestSolution, var) + 0.5);

    for (i = 0; i < sc_size; i++) {
        for (k = 0; k < sum_sc[i + 1] - sum_sc[i]; k++) {
            sum = 0;
            for (l = 0; l < string_size; l++) {
                int pos = value_symbol(alphabet, alpha_size, sc[i][k + l]);
                if (pos != -1)
                    sum += x[l][pos];
            }
            if (y[i][k] == 1 && sum > dc) {
                fprintf(file, "Problema na restricao (2) na string %d posicao %d.\n", i, k);
                valid = false;
            }
        }
    }

    for (i = 0; i < sf_size; i++) {
        for (k = 0; k < sum_sf[i + 1] - sum_sf[i]; k++) {
            sum = 0;
            for (l = 0; l < string_size; l++) {
                int pos = value_symbol(alphabet, alpha_size, sf[i][k + l]);
                if (pos != -1)
                    sum += x[l][pos];
            }
            if (sum < df) {
                fprintf(file, "Problema na restricao (3) na string %d posicao %d.\n", i, k);
                valid = false;
            }
        }
    }

    for (i = 0; i < sc_size; i++) {
        sum = 0;
        for (k = 0; k < sum_sc[i + 1] - sum_sc[i]; k++) {
            sum += y[i][k];
        }
        if (sum < 1) {
            fprintf(file, "Problema na restricao (4) na string %d.\n", i);
            valid = false;
        }
    }

    input = fopen(filename, "r");
    if (!input) exit(3);

    strcat(sol_x, "");
    for(i = 0; i < string_size; i++) {
        for (j = 0; j < alpha_size; j++) {
            var = vars[i * alpha_size + j];
            solval = SCIPgetSolVal(scip, bestSolution, var);
            if(solval == 0.0)
            {
                aux[0] = symbol_value(alphabet, alpha_size, j);
                aux[1] = '\0';
                strcat(sol_x, aux);
            }
        }
    }

    if (valid){
        printf("Solucao validada com sucesso!"); fflush(stdout);}
    else{
        printf("Falha na validacao da solucao!"); fflush(stdout);}

    for (i = 0; i < string_size; i++) {
        free(x[i]);
    }
    free(x);

    for (i = 0; i < sc_size; i++) {
        free(y[i]);
    }
    free(y);

    fclose(file);
    fclose(input);
}
