#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2014 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic Licence.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#@file    Makefile
#@brief   Makefile for SCIP project
#@author  Timo Berthold
#@author  Stefan Heinz

#-----------------------------------------------------------------------------
# path
#-----------------------------------------------------------------------------

#Lexa
#SCIPDIR         =       /usr/local/scipoptsuite-3.2.1/scip-3.2.1

#Edna
#SCIPDIR         =       /usr/local/scip/scipoptsuite-3.2.1/scip-3.2.1

#Servidor
#SCIPDIR         =     /usr/local/scip
SCIPDIR          =     ~/Programming/scipoptsuite-7.0.1/
#PC
#SCIPDIR         =       /home/ubuntu/environment/scipoptsuite-3.2.1/scip-3.2.1
# ./bin/dssp-NOTRACE-GLOBAL-STRING-RA-DSP-BASIC -f data/instancias-tcc/grupo1/rand-4-15-15-500-1.dssp

#-----------------------------------------------------------------------------
# OPTIONS
#-----------------------------------------------------------------------------

TRACE         =       NOTRACE
# DEBUG, DEBUG_VNS
#EXECUTE       =       GLOBAL
# ONLYROOT or GLOBAL
#PROBLEM       =       SUB
# SUB or STRING
#HEUR       =       NONE
# RA , BCPA or NONE
#TYPE       = 	   DSP
#DSP, CSP, FSP
#VNS_TYPE = NOTYPE
#BLPL, SBPL, VBPL
#RA_TYPE = ROUNDING
#REDUCE, MERGE 
#MODEL = VELHO
#NOVO
#TIMEH = TIMEX
#SSS=5, SS=30, S=60, M=120, L=180, LL=300
#NSIZE = NSIZEM
#S=2, M=3, L=4, LL=5
#OMEGA = OMEGAM
#S=2, M=3, L=4, LL=5
#SBPL_SEMY
#IMPROVEMENT=BEST
#BEST, FIRST





# make EXECUTE=ONLYROOT PROBLEM=SUB HEUR=RA
# scip ONLYROOT STRING RA ../instancias/rand0-2-5-5-250-1.dssp 

#-----------------------------------------------------------------------------
# include default project Makefile from SCIP
#-----------------------------------------------------------------------------

include $(SCIPDIR)/make/make.project

#-----------------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------------

MAINNAME	=	dssp-$(TRACE)
#-$(EXECUTE)-$(PROBLEM)-$(HEUR)-$(TYPE)-$(VNS_TYPE)-$(RA_TYPE)-$(MODEL)-$(TIMEH)-$(NSIZE)-$(OMEGA)-$(IMPROVEMENT)
CMAINOBJ	=	cmain.o \
			probdata_dssp.o \
			reader_dssp.o \
			utils.o \
			ra.o \
			bcpa.o \
			vns.o
#			reader_csp.o  \
#			heur_myheuristic.o

CXXMAINOBJ	=

MAINSRC		=	$(addprefix $(SRCDIR)/,$(CMAINOBJ:.o=.c))
MAINSRC		+=	$(addprefix $(SRCDIR)/,$(CXXMAINOBJ:.o=.cpp))
MAINDEP		=	$(SRCDIR)/depend.cmain.$(OPT)

MAIN		=	$(MAINNAME).$(BASE).$(LPS)$(EXEEXTENSION)
MAINFILE	=	$(BINDIR)/$(MAIN)
MAINSHORTLINK	=	$(BINDIR)/$(MAINNAME)
MAINOBJFILES	=	$(addprefix $(OBJDIR)/,$(CMAINOBJ))
MAINOBJFILES	+=	$(addprefix $(OBJDIR)/,$(CXXMAINOBJ))

#-----------------------------------------------------------------------------
# External libraries
#-----------------------------------------------------------------------------

FLAGS		+=
LDFLAGS		+=

#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

ifeq ($(VERBOSE),false)
.SILENT:	$(MAINFILE) $(MAINOBJFILES) $(MAINSHORTLINK)
endif

.PHONY: all
all:            $(SCIPDIR) $(MAINFILE) $(MAINSHORTLINK)

.PHONY: lint
lint:		$(MAINSRC)
		-rm -f lint.out
		$(SHELL) -ec 'for i in $^; \
			do \
			echo $$i; \
			$(LINT) $(SCIPDIR)/lint/scip.lnt +os\(lint.out\) -u -zero \
			$(FLAGS) -UNDEBUG -UWITH_READLINE -UROUNDING_FE $$i; \
			done'
.PHONY: scip
scip:
		@$(MAKE) -C $(SCIPDIR) libs $^

.PHONY: doc
doc:
		@-(cd doc && ln -fs ../$(SCIPDIR)/doc/scip.css);
		@-(cd doc && ln -fs ../$(SCIPDIR)/doc/pictures/scippy.png);
		@-(cd doc && ln -fs ../$(SCIPDIR)/doc/pictures/miniscippy.png);
		@-(cd doc && ln -fs ../$(SCIPDIR)/doc/scipfooter.html footer.html);
		cd doc; $(DOXY) $(MAINNAME).dxy

$(MAINSHORTLINK):	$(MAINFILE)
		@rm -f $@
		cd $(dir $@) && ln -s $(notdir $(MAINFILE)) $(notdir $@)

$(OBJDIR):
		@-mkdir -p $(OBJDIR)

$(BINDIR):
		@-mkdir -p $(BINDIR)

.PHONY: clean
clean:
#		rm -r bin
		rm -r obj
# $(OBJDIR)
# ifneq ($(OBJDIR),)
# 		@-(rm -f $(OBJDIR)/*);
# 		@echo "-> remove main objective files"
# endif
# 		@-rm -f $(MAINFILE) $(MAINLINK) $(MAINSHORTLINK)
# 		@echo "-> remove binary"

.PHONY: test
test:           $(MAINFILE)
		@-(cd check && ln -fs ../$(SCIPDIR)/check/check.sh);
		@-(cd check && ln -fs ../$(SCIPDIR)/check/evalcheck.sh);
		@-(cd check && ln -fs ../$(SCIPDIR)/check/evalcheck_cluster.sh);
		@-(cd check && ln -fs ../$(SCIPDIR)/check/check.awk);
		@-(cd check && ln -fs ../$(SCIPDIR)/check/getlastprob.awk);
		@-(cd check && ln -fs ../$(SCIPDIR)/check/configuration_set.sh);
		@-(cd check && ln -fs ../$(SCIPDIR)/check/configuration_logfiles.sh);
		@-(cd check && ln -fs ../$(SCIPDIR)/check/configuration_tmpfile_setup_scip.sh);
		@-(cd check && ln -fs ../$(SCIPDIR)/check/run.sh);
		cd check; \
		$(SHELL) ./check.sh $(TEST) $(MAINFILE) $(SETTINGS) $(notdir $(MAINFILE)) $(TIME) $(NODES) $(MEM) $(THREADS) $(FEASTOL) $(DISPFREQ) $(CONTINUE) $(LOCK) "example" $(LPS) $(VALGRIND) $(CLIENTTMPDIR) $(OPTCOMMAND);

.PHONY: tags
tags:
		rm -f TAGS; ctags -e src/*.c src/*.h $(SCIPDIR)/src/scip/*.c $(SCIPDIR)/src/scip/*.h;

.PHONY: depend
depend:		$(SCIPDIR)
		$(SHELL) -ec '$(DCC) $(FLAGS) $(DFLAGS) $(MAINSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/\([0-9A-Za-z\_]*\).c|$$\(OBJDIR\)/\2.o: $(SRCDIR)/\2.c|g'\'' \
		>$(MAINDEP)'

-include	$(MAINDEP)

$(MAINFILE):	$(BINDIR) $(OBJDIR) $(SCIPLIBFILE) $(LPILIBFILE) $(NLPILIBFILE) $(MAINOBJFILES)
		@echo "-> linking $@"
		$(LINKCXX) $(MAINOBJFILES) \
		$(LINKCXX_L)$(SCIPDIR)/lib $(LINKCXX_l)$(SCIPLIB)$(LINKLIBSUFFIX) \
                $(LINKCXX_l)$(LPILIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(NLPILIB)$(LINKLIBSUFFIX) \
                $(OFLAGS) $(LPSLDFLAGS) \
		$(LDFLAGS) $(LINKCXX_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.c
		@echo "-> compiling $@"
		echo "$(CC) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CFLAGS) -g -D$(TRACE)  -c -std=c99 -Wno-sign-compare -Wno-unused-result $< $(CC_o)$@"
# -D$(EXECUTE) -D$(PROBLEM) -D$(HEUR) -D$(TYPE) -D$(VNS_TYPE) -D$(RA_TYPE) -D$(MODEL) -D$(TIMEH) -D$(NSIZE) -D$(OMEGA) -D$(IMPROVEMENT) -c -std=c99 -Wno-sign-compare -Wno-unused-result $< $(CC_o)$@"
		$(CC) $(FLAGS) $(BINOFLAGS) $(CFLAGS) -g -c -std=c99 -pedantic  -D$(TRACE)  -Wno-sign-compare -Wno-unused-result $< $(CC_o)$@
#-D$(EXECUTE) -D$(PROBLEM) -D$(HEUR) -D$(TYPE) -D$(VNS_TYPE) -D$(RA_TYPE) -D$(MODEL) -D$(TIMEH) -D$(NSIZE) -D$(OMEGA) -D$(IMPROVEMENT) -Wno-sign-compare -Wno-unused-result $< $(CC_o)$@
#		$(CC) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CFLAGS) -g -c -std=c99 -pedantic  -D$(TRACE) -Wno-sign-compare -Wno-unused-result $< $(CC_o)$@
#-D$(EXECUTE) -D$(PROBLEM) -D$(HEUR) -D$(TYPE) -D$(VNS_TYPE) -D$(RA_TYPE) -D$(MODEL) -D$(TIMEH) -D$(NSIZE) -D$(OMEGA) -D$(IMPROVEMENT) -Wno-sign-compare -Wno-unused-result $< $(CC_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp
		@echo "-> compiling $@"
		$(CXX) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CXXFLAGS) -g -c -std=c99 -Wdeclaration-after-statement -Wno-sign-compare -Wno-unused-result $< $(CXX_o)$@

#---- EOF --------------------------------------------------------------------
