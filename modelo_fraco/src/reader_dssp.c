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

/**@file   reader_dssp.c
 * @brief  DSSP problem reader file reader
 * @author Edna A. Hoshino
 * @author Emerson Silva
 * @author Jean Torres
 *
 * This file implements the reader/parser used to read the DSSP input data. For more details see \ref READER.
 *
 * @page READER Parsing the input format and creating the problem
 *
 * This READER process files with extension ".dssp".
 *
 * @section READERINCLUDE The SCIPincludeReaderDSSP() interface method
 *
 * The interface method <code>SCIPincludeReaderDSSP()</code> is called to add the reader plugin to \SCIP (see
 * cmain.c). This means \SCIP gets informed that this reader is available for reading input files. Therefore, the
 * function <code>SCIPincludeReader()</code> is called within this method which passes all necessary information of the
 * reader to SCIP. This information includes the name of the reader, a description, and the file extension for which the
 * file reader is in charge. In our case we selected the file extension "dssp". This means that all files which have
 * this file extension are passed to our reader for parsing. Besides these information the call
 * <code>SCIPincludeReader()</code> also passes for each callback of the reader a function pointers
 * (some of them might be NULL pointers). These function
 * pointers are used by \SCIP to run the reader. For more information about all available reader callbacks we refer to
 * the <a href="http://scip.zib.de/doc/html/READER.html">How to add file readers</a> tutorial. In the remaining section
 * we restrict ourself to the callback <code>READERREAD</code> which is the only one we implemented for the dssp
 * example. All other callbacks are not required for this example.
 *
 * @section READERREAD The READERREAD callback method
 *
 * The READERREAD callback is in charge of parsing a file and creating the problem. To see the list of arguments this
 * functions gets see the file type_reader.h in the source of \SCIP. The following arguments are of interest in our
 * case. First of all the \SCIP pointer, the file name, and the SCIP_RESULT pointer. The \SCIP pointer gives us the
 * current environment. The file name states the file which we should open and parse. Last but not least, the SCIP_RESULT
 * pointer is required to tell \SCIP if the parsing process was successfully or
 * not. Note that in type_reader.h you also find a list of allowable result values for the SCIP_RESULT pointer and the
 * <code>SCIP_RETCODE</code> which is the return value of this function.
 *
 * @subsection PARSING Parsing the problem
 *
 * The file can be opened and parsed with your favorite methods. In this case we are using the functionality provided by
 * \SCIP since this has some nice side effects. We are using the function SCIPfopen() which can besides standard
 * files also handle files which are packed. To find all files related to the parsing of a file, we refer to the file pub_misc.h
 * in the source of SCIP. Parsing the data out of the file is not that hard. Please look at the code and comments
 * therein for more details.
 *
 * @subsection CREATING Creating the problem
 *
 * After parsing the file the final task for the reader is to create the problem. In our case, we pass the collected data
 * to the \ref probdata_dssp.h "main problem data plugin". For this, we use the interface methods
 * SCIPprobdataCreate() which is provided by the
 * problem data plugin (see probdata_dssp.c). After that, the reader sets the result value for the SCIP_RESULT
 * pointer to <code>SCIP_SUCCESS</code> and returns with a proper <code>SCIP_RETCODE</code>.
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


/*#define SCIP_DEBUG*/

#include <assert.h>
#include <string.h>

#include "scip/cons_setppc.h"

#include "probdata_dssp.h"
#include "reader_dssp.h"
#include "parameters_dssp.h"

/**@name Reader properties
 *
 * @{
 */

#define READER_NAME             "dsspreader"
#define READER_DESC             "file reader for dssp data format"
#define READER_EXTENSION        "dssp"

/**@} */

int getNextLine(SCIP_FILE* file, int* lineno, char* buffer);

/** Get next not empty line of the filestream and return it into the buffer */
int getNextLine(SCIP_FILE* file, int* lineno, char* buffer)
{
    do{
        /* ignore empty line */
        if(SCIPfgets(buffer, SCIP_MAXSTRLEN, file) == NULL)
            return 0;
        (*lineno)++;
    } while (buffer[0] == '\r' || buffer[0] == '\n');

    return 1;
}

/**@name Callback methods
 *
 * @{
 */

/** instance problem reader method */
static
SCIP_DECL_READERREAD(readerReaddssp)
{
   return readerDssp(scip, filename, result);
}

int readerDssp(SCIP* scip, const char* filename, SCIP_RESULT* result)
   
{  /*lint --e{715}*/
    SCIP_FILE* file;
    SCIP_Bool error = FALSE;

    char buffer[SCIP_MAXSTRLEN];

    /* dssp data */
    int sc_size = 0, sf_size = 0, alpha_size = 0, string_size = 0, kc = 0, kf = 0, *sum_sc = NULL, *sum_sf = NULL;
    char **sc = NULL, **sf = NULL, *alphabet = NULL;
    int i, count, nread, lineno;

    *result = SCIP_DIDNOTRUN;

    /* open file */
    file = SCIPfopen(filename, "r");
    if(file == NULL)
    {
        /*SCIPerrorMessage("cannot open file <%s> for reading\n", filename);*/
        SCIPprintSysError(filename);
        return SCIP_NOFILE;
    }

    lineno = 0;

    /* read problem dimension */
    if(!SCIPfeof(file))
    {
        /* get next line */
        if(!getNextLine(file, &lineno, buffer))
            return SCIP_READERROR;

        /* Le tamanho do alfabeto */
        nread = sscanf(buffer, "%d\n", &alpha_size);

        if(nread < 1)
        {
            SCIPwarningMessage(scip, "invalid input line %d in file <%s>: <%s>\n", lineno, filename, buffer);
            return SCIP_READERROR;
        }

        /* get next line */
        if(!getNextLine(file, &lineno, buffer))
            return SCIP_READERROR;

        /* Le tamanho dos conjuntos */
        nread = sscanf(buffer, "%d %d\n", &sc_size, &sf_size);

        if(param.problem==CSP){
           //#ifdef CSP
           sf_size = 0;
           //#endif
        }
        if(param.problem==FSP){
           //#ifdef FSP
           sc_size = 0;
           //#endif
        }

        // printf("\n\nSC_SIZE %d \n", sc_size); getchar();
        if(nread < 2)
        {
            SCIPwarningMessage(scip, "invalid input line %d in file <%s>: <%s>\n", lineno, filename, buffer);
            return SCIP_READERROR;
        }

        /* get next line */
        if(!getNextLine(file, &lineno, buffer))
            return SCIP_READERROR;

        /* Le tamanho da string */
        nread = sscanf(buffer, "%d\n", &string_size);
        // kc = string_size;

        if(nread < 1)
        {
            SCIPwarningMessage(scip, "invalid input line %d in file <%s>: <%s>\n", lineno, filename, buffer);
            return SCIP_READERROR;
        }

        /* get next line */
        if(!getNextLine(file, &lineno, buffer))
            return SCIP_READERROR;

        /* Le constantes kc e kf */
        nread = sscanf(buffer, "%d %d\n", &kc, &kf);

        if(param.substring){
           //#ifdef SUB
           if (string_size == kc) {
              string_size = (int) (string_size*0.75);
              kc = (int) (kc*0.75);
           }
        //#endif
        }
        
        if(nread < 2)
        {
            SCIPwarningMessage(scip, "invalid input line %d in file <%s>: <%s>\n", lineno, filename, buffer);
            return SCIP_READERROR;
        }

        /* allocate buffer memory for storing the strings in Sc */
        SCIP_CALL(SCIPallocMemoryArray(scip, &sc, sc_size));
        for (i = 0; i < sc_size; i++) {
            SCIP_CALL(SCIPallocMemoryArray(scip, &(sc[i]), SCIP_MAXSTRLEN));
        }

        /* allocate buffer memory for storing the strings in Sf */
        SCIP_CALL(SCIPallocMemoryArray(scip, &sf, sf_size));
        for (i = 0; i < sf_size; i++) {
            SCIP_CALL(SCIPallocMemoryArray(scip, &(sf[i]), SCIP_MAXSTRLEN));
        }

        /* allocate buffer memory for storing the strings size in Sc */
        SCIP_CALL(SCIPallocMemoryArray(scip, &sum_sc, sc_size + 1));
        sum_sc[0] = 0;

        /* allocate buffer memory for storing the strings size in Sf */
        SCIP_CALL(SCIPallocMemoryArray(scip, &sum_sf, sf_size + 1));
        sum_sf[0] = 0;

        /* allocate buffer memory for storing the alphabet characters */
        SCIP_CALL(SCIPallocMemoryArray(scip, &alphabet, alpha_size));

        /* read alphabet */
        for (i = 0; i < alpha_size && !SCIPfeof(file) && !error; i++)
        {
            /* get next line */
            if(!getNextLine(file, &lineno, buffer)){
                error = TRUE;
                break;
            }

            nread = sscanf(buffer, "%c\n", &alphabet[i]);

            if(nread < 1)
            {
                SCIPwarningMessage(scip, "invalid input line %d in file <%s>: <%s>\n", lineno, filename, buffer);
                error = TRUE;
                break;
            }
        }

        /* parse Sc */
        for (i = 0; i < sc_size && !SCIPfeof(file) && !error; i++)
        {
            /* get next line */
            if(!getNextLine(file, &lineno, buffer)){
                error = TRUE;
                break;
            }

            /* parse Sc strings */
            nread = sscanf(buffer, "%s\n", sc[i]);
            if(nread < 1)
            {
                SCIPwarningMessage(scip, "invalid input line %d in file <%s>: <%s>\n", lineno, filename, buffer);
                error = TRUE;
                break;
            }

            // Soma no vetor sum_sc
            sum_sc[i + 1] = sum_sc[i] + strlen(sc[i]) - string_size + 1;
        }

        /* parse Sf */
        count = 1;
        for (i = 0; i < sf_size && !SCIPfeof(file) && !error; i++)
        {
            /* get next line */
            if(!getNextLine(file, &lineno, buffer)){
                error = TRUE;
                break;
            }

            /* parse Sf strings */
            nread = sscanf(buffer, "%s\n", sf[i]);
            if(nread < 1)
            {
                SCIPwarningMessage(scip, "invalid input line %d in file <%s>: <%s>\n", lineno, filename, buffer);
                error = TRUE;
                break;
            }
            count++;

            // Soma no vetor sum_sf
            sum_sf[i + 1] = sum_sf[i] + strlen(sf[i]) - string_size + 1;
        }

        (void)SCIPfclose(file);

        if(count < sf_size)
        {
            SCIPwarningMessage(scip, "error in parsing line <%d> in file <%s>: <%s>\n", lineno, filename, buffer);
            return SCIP_READERROR;
        }

        /* create a new problem in SCIP */
        SCIP_CALL(SCIPprobdataCreate(scip, filename, sc_size, sf_size, alpha_size, string_size, kc, kf, sc, sf, sum_sc, sum_sf, alphabet));
    }

    if(error)
        return SCIP_READERROR;

    *result = SCIP_SUCCESS;

    return SCIP_OKAY;
}

/**@} */


/**@name Interface methods
 *
 * @{
 */

/** includes the dssp file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderDSSP(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
    SCIP_READERDATA* readerdata;
    SCIP_READER* reader;

    /* create dssp reader data */
    readerdata = NULL;

    SCIP_CALL(SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, readerdata));
    assert(reader != NULL);

    SCIP_CALL(SCIPsetReaderRead(scip, reader, readerReaddssp));

    return SCIP_OKAY;
}

/**@} */
