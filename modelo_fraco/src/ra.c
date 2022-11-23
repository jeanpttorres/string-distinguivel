#include<string.h>
#include "ra.h"
#include "probdata_dssp.h"
#include "parameters_dssp.h"

int heur_RA(SCIP* scip, char *filename, double tempo, double *value, int *heur, int *pdc, int *pdf)
{
   int i, j, k, l, string_size, alpha_size;//, *heur;
    int ra;
    char name[SCIP_MAXSTRLEN], *alphabet;
    FILE *file;
    int         sc_size;
    int*        sum_sc;
    double maior;

    SCIP_SOL *bestSolution;
    SCIP_PROBDATA *probdata;
    SCIP_VAR **vars, *var;
    SCIP_Real solval;

    assert(scip != NULL);
    
    bestSolution = SCIPgetBestSol(scip);
  
    probdata = SCIPgetProbData(scip);
    assert(probdata != NULL);

    // nvars = SCIPprobdataGetNVars(probdata);
    vars = SCIPprobdataGetVars(probdata);
    string_size = SCIPprobdataGetString_size(probdata);
    alphabet = SCIPprobdataGetAlphabet(probdata);
    alpha_size = SCIPprobdataGetAlpha_size(probdata);
    sc_size = SCIPprobdataGetSc_size(probdata) == 0 ? SCIPprobdataGetSf_size(probdata) : SCIPprobdataGetSc_size(probdata);
    sum_sc = SCIPprobdataGetSum_sc(probdata);

    strcpy(name, filename);
    strcat(name, "-ra.sol");

    file = fopen(name, "w");
    
    if(!file)
    {
        printf("\nProblem to create solution file: %s\n", name);
        return 0;
    }

    name[0] = '\0';

    //    heur = (int *)malloc(sizeof(int) * string_size);
    if(param.ra_novo){
       //#ifdef NOVO
       double *guarda_maior; int change=0, ind;
       guarda_maior = malloc (string_size * sizeof (double));
       for(i=0;i<string_size;i++)guarda_maior[i]=-1.0;
       for(i=0;i<string_size;i++)heur[i]=0;
       
       int kcounter = 0, kpos = 0;
       //printf("REG2\n"); fflush(stdout);

       if(param.ra_merge){
          //#ifdef MERGE    
          /**
             Encontra o Melhor alinhamento em cada String em SC
             Para estes alinhamentos, encontra a variável de melhor valor em cada posição para computar 
             o alinhamento final
          */
          
          // printf("\nMERGE %d\n", sc_size); 
          
          double menor = 2;
          int melhor_posicao, jj;
          
          int *melhor_alinhamento = (int*)malloc(sizeof(int)*sc_size);
          double *melhor_alinhamento_valor = (double*)malloc(sizeof(double)*string_size);
          kcounter = 0;
          double menor_soma = 0, soma=0;;
          
          for (i = 0; i < sc_size; i ++) {
             menor_soma = 0;//2*string_size;
             //        printf("%I: %d\n", i);fflush(stdout);
             for (j = 0; j < string_size; j ++) {
                soma = 0;
                for (k = 0; k < sum_sc[i+1] - sum_sc[i]; k++) {
                   menor = 0; //2;
                   for (l = 0; l < alpha_size; l++) {
                      var = vars[kcounter++];
                      if(j+k > (3/4.0)*string_size){
                         continue;
                      }
                      solval = SCIPgetSolVal(scip, bestSolution, var); // guardar numa variavel
                      
                      //                    printf("%lf\n", solval);fflush(stdout);
                      
                      if (solval > menor + 0.00001) { //< menor) {
                         // printf("MAIOR %d %.2lf\n", l, solval); fflush(stdout);
                         melhor_posicao = j;
                         menor = solval;
                         ind = kcounter;
                    }
                   }
                   soma += menor;
                }
                if (soma > menor_soma + 0.00001) {
                   menor_soma = soma; 
                   melhor_alinhamento[i] = j;
                   //printf("\nMelhor soma da string %d eh na posicao %d com valor %lf", i, melhor_alinhamento[i], menor_soma);
                }
             }
          }

          
          for (i = 0; i < sc_size; i ++) //printf("%d ", melhor_alinhamento[i]);

             for (i = 0; i < string_size; i ++ ) melhor_alinhamento_valor[i] = 0;//2;

          kcounter=0;
          for (i = 0; i < sc_size; i ++) {
             for (j = 0; j < string_size; j ++) {
                for (k = 0; k < sum_sc[i+1] - sum_sc[i]; k++) {
                   menor = 0;//2;
                   for (l = 0; l < alpha_size; l++) {
                      var = vars[kcounter++];
                      solval = SCIPgetSolVal(scip, bestSolution, var); // guardar numa variavel
                      
                      if (solval > menor) {
                         // printf("MAIOR %d %.2lf\n", l, solval); fflush(stdout);
                         menor=solval;
                         ind = l;
                      }
                   }
                   
                   if (j == melhor_alinhamento[i]) {
                      if (menor > melhor_alinhamento_valor[j]) {
                         heur[j] = l;
                         melhor_alinhamento_valor[j] = menor;
                      }
                   }
                }
                
             }
          }
          // printf("\n"); 
          
          free(melhor_alinhamento);
       }
       else{
          if(param.ra_reduce){
             //#elif REDUCE
             /**
                Ordena as variáveis
                Seleciona a de melhor valor
                Torna 1 as variáveis que se tornam inviáveis quando a variável acima é escolhida
             */
             
             double menor = 2;
             int melhor_posicao, jj;
             
             int *melhor_alinhamento = (int*)malloc(sizeof(int)*sc_size);
             kcounter = 0;
             int menor_soma = 0;
             
             double *solutions, *positions;
             
             for (int count = 0; count < string_size; count +=1) {
                for (i = 0; i < sc_size; i ++) {
                   for (j = 0; j < string_size; j ++) {
                      for (k = 0; k < sum_sc[i+1] - sum_sc[i]; k++) {
                         menor = 2;
                         for (l = 0; l < alpha_size; l++) {
                            var = vars[kcounter++];
                            solval = SCIPgetSolVal(scip, bestSolution, var);
                            
                            if (solval < menor) {
                               
                               melhor_posicao = j;
                               maior=solval;
                               ind = l;
                            }
                         }                                             
                      }
                   }
                }
             }
          }
          else{
             //#else
             
             int melhor_posicao, jj;
             
             double *melhor_alinhamento_valor = (double*)malloc(sizeof(double)*string_size);
             
             for (i = 0; i < string_size; i ++) melhor_alinhamento_valor[i] = 2;
             kcounter = 0;
             int menor_soma = 0; 
             double menor;
             for (i = 0; i < sc_size; i ++) {
                for (j = 0; j < string_size; j ++) {
                   for (k = 0; k < sum_sc[i+1] - sum_sc[i]; k++) {
                      menor = 2;
                      for (l = 0; l < alpha_size; l++) {
                         var = vars[kcounter++];
                         solval = SCIPgetSolVal(scip, bestSolution, var); // guardar numa variavel
                         
                         if (solval < menor) {
                            // printf("MAIOR %d %.2lf\n", l, solval); fflush(stdout);
                            melhor_posicao = j;
                            menor=solval;
                            ind = l;
                         }
                      }
                      
                      if (menor < melhor_alinhamento_valor[j]) {
                         heur[j] = l;
                         melhor_alinhamento_valor[j] = menor;
                      }
                   }
                }
             }
          //#endif
          }
       }
    }
    else{
       //#else
       for(i = 0; i < string_size; i++) {
          double menor = 2.0;
          int ind = -1;
          for (j = 0; j < alpha_size; j++) {
             var = vars[i * alpha_size + j];
             solval = SCIPgetSolVal(scip, bestSolution, var);
             //fprintf(file, "%.2lf ", solval);
             if (solval < menor) {
                menor = solval;
                ind = j;
             }
          }
          heur[i] = ind;
       }
       //#endif
    }
    // for (i = 0; i < string_size; i ++) {
    //     printf("%d ", heur[i]); fflush(stdout);
    // } printf("\n");
    
    ra =  valorSolucao(scip, filename, heur, pdc, pdf);
    fprintf(file, "Time: %lfs\n", tempo);
    fprintf(file, "Objective Value(dc - df) = %d\n", ra);
    fprintf(file, "Target String(x) = ");
    
    for(i = 0; i < string_size; i++) {
       fprintf(file, "%c", symbol_value(alphabet, alpha_size, heur[i]));
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
          if (heur[i] != j) {
             fprintf(file, "0 ");
          } else {
             fprintf(file, "1 ");
          }
       }
       fprintf(file, "\n");        
    }
    //printf("#3\n"); fflush(stdout);
    fclose(file);
    
    *value = ra;
    return 1;       
}
