##################################################
##  A pgi makefile for ASP on rufus		##
##  Matt Alvarado, 02/13/2012		        ## 
##################################################
SHELL = /bin/sh
##
## === Define files and directories
##

# Source tree defaults to the directory holding this Makefile.
# Build outputs (objects, .mod files, library) default to the same tree
# but can be redirected with ASP_BUILD_DIR (run make from that directory), e.g.
#   cd build && make -f .../SRC_TUV/ASP/Makefile all ASP_BUILD_DIR=$PWD ASP_FC=ifx
ASP_DIR ?= $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
ASP_BUILD_DIR ?= $(ASP_DIR)
ASP_LIB = $(ASP_BUILD_DIR)/libAspLibrary.a
ASP_INCLUDES_DIR = $(ASP_DIR)/includes
AR = ar
AR_FLAGS = -rv

PRNTDIR = $(ASP_DIR)/OUTPUT/
INCDIR  = $(ASP_DIR)/INCLUDE/
EXEDIR  = $(ASP_DIR)
MODULES_OUTPUT_DIR = $(ASP_BUILD_DIR)/ModuleFiles/
# Create the object and module directories (they are not kept in git)
$(shell mkdir -p $(ASP_BUILD_DIR)/ModuleFiles $(ASP_BUILD_DIR)/ObjectFiles)

# ASP source: the common_ASP submodule at src/ by default
ASP_SRC_DIR ?= $(ASP_DIR)/src
SOURCEDIR = $(ASP_SRC_DIR)/
CHEMDIR = $(ASP_SRC_DIR)/Chemistry/
COAGDIR = $(ASP_SRC_DIR)/Coagulation/
COMPDIR = $(ASP_SRC_DIR)/ComputationalSubroutines/
CONDDIR = $(ASP_SRC_DIR)/Condensation/
EULERDIR = $(ASP_SRC_DIR)/EulerianEnvironment/
INFDIR = $(ASP_SRC_DIR)/InfrastructuralCode/
PARTDIR = $(ASP_SRC_DIR)/LagrangianParticles/
OUTPUTDIR = $(ASP_SRC_DIR)/OutputRoutines/
UNIFACDIR = $(ASP_SRC_DIR)/UNIFAC/
OBJDIR = $(ASP_BUILD_DIR)/ObjectFiles/

FILES	= $(SOURCEDIR)ModelParameters.f90 \
	  $(SOURCEDIR)ReadInputDeck.f90  \
          $(SOURCEDIR)StepASP4SAM.f90 \
	  $(SOURCEDIR)Time.f90 \
	  $(SOURCEDIR)EulerianBoxDriver.f90

MODULES = $(CHEMDIR)Chemistry.f90 \
          $(COAGDIR)Coagulation.f90 \
          $(CONDDIR)Condensation.f90 \
	  $(EULERDIR)GridPointFields.f90 \
          $(INFDIR)InfrastructuralCode.f90 \
	  $(PARTDIR)LagrangianAerosols.f90 \
          $(OUTPUTDIR)OutputRoutines.f90 \
	  $(UNIFACDIR)unifac.f90 \
          $(UNIFACDIR)soa_unidriv.f90 \
          $(UNIFACDIR)soa_unidriva.f90 \
	  $(COMPDIR)CDFNor.f90 \
          $(COMPDIR)dgaus8.f \
	  $(COMPDIR)Gibbs-DuhemEvaluationFunction.f90 \
	  $(COMPDIR)indexx.f \
          $(COMPDIR)lsodes.f \
          $(COMPDIR)PPND16.f
	 
CHEMINC = $(CHEMDIR)AqPhaseChemistryInits.h  $(CHEMDIR)AqPhaseReactionIntegration.h  \
	  $(CHEMDIR)ChemicalPropertyRoutines.h  $(CHEMDIR)ChemistryDerivativeStructures.h \
	  $(CHEMDIR)ChemistryParametersInitRoutines.h  $(CHEMDIR)EvolveGasChemistry.h \
	   $(CHEMDIR)GasPhaseChemistryInits.h  $(CHEMDIR)OrgPhaseChemInits.h

COAGINC = $(COAGDIR)CoagulationKernels.h

CONDINC = $(CONDDIR)AqOrgCondensationIntegrator.h \
          $(CONDDIR)OrgPhaseReactions.h \
          $(CONDDIR)ParticleAcidity.h \
          $(CONDDIR)binsolu.h \
	  $(CONDDIR)CondensationInitialization.h  \
          $(CONDDIR)CondensationIntegrator.h \
	  $(CONDDIR)CondensationRelatedFunctions.h  \
	  $(CONDDIR)HydrophobicCondensationFunctions.h \
	  $(CONDDIR)OrgCondensationIntegrator.h

EULERINC = $(EULERDIR)GasPhaseChemistry.h \
           $(EULERDIR)Initialize.h  \
	   $(EULERDIR)InterpolateGrid.h \
           $(EULERDIR)Pressure.h  \
           $(EULERDIR)SetGrid.h \
	   $(EULERDIR)Temperature.h \
           $(EULERDIR)Water.h

INFINC = $(INFDIR)Counters.h  \
         $(INFDIR)InputOutputCommands.h \
	 $(INFDIR)StringCommands.h

PARTINC = $(PARTDIR)InitializationAndSamplingRoutines.h  \
          $(PARTDIR)ParticleAttributes.h \
	  $(PARTDIR)ParticleStructure.h  \
          $(PARTDIR)SCTailorMakeParticle.h \
	  $(PARTDIR)SortAerosol.h \
          $(PARTDIR)Thermodynamics.h


UNIFACINC = $(UNIFACDIR)unifacparam.h \
            $(UNIFACDIR)unifacparamMattAq.h

ALLINC = $(CHEMINC) $(COAGINC) $(CONDINC) $(EULERINC) $(INFINC) $(PARTINC) $(UNIFACINC)

FILESINC = $(ALLINC)

FILESDUMP = asp 
	
##
## === Compile, link and execute:
## 

# Compile command -- single processor

# Compiler: ifort by default; set ASP_FC=ifx for the Intel oneAPI compiler.
# (FC has a built-in make default, so a separate name is used.)
ASP_FC ?= ifort
ifeq ($(ASP_FC),ifx)
  ASP_FC_OPTS ?= -DIFORT -fp-model=precise
else
  ASP_FC_OPTS ?= -DIFORT -fp-stack-check -fp-model source
endif
FC = $(ASP_FC) $(ASP_FC_OPTS)

LDFLAGS := -Wl,-zmuldefs

FFLAGS :=  -fpp -FI -O -r8 -prec-div -g -traceback \
	   -I$(MODULES_OUTPUT_DIR) -module $(MODULES_OUTPUT_DIR)

#FFLAGS :=  -fpp -FI -r8 -prec_div -g -O2 -fpe0 -debug extended -traceback \
          -I$(MODULES_OUTPUT_DIR) -module $(MODULES_OUTPUT_DIR)
F90FLAGS := $(FFLAGS) -FR -DIFORT


OBJECTS :=
OBJECTS := $(OBJECTS) \
           $(OBJDIR)CDFNor.o $(OBJDIR)indexx.o $(OBJDIR)lsodes.o \
           $(OBJDIR)PPND16.o \
	   $(OBJDIR)dgaus8.o $(OBJDIR)DMiLay.o $(OBJDIR)MIEV0.o $(OBJDIR)MieSupport.o
OBJECTS := $(OBJECTS) \
	   $(OBJDIR)ModelParameters.o \
	   $(OBJDIR)InfrastructuralCode.o \
	   $(OBJDIR)GridPointFields.o \
	   $(OBJDIR)Time.o \
	   $(OBJDIR)Chemistry.o \
	   $(OBJDIR)Gibbs-DuhemEvaluationFunction.o \
           $(OBJDIR)LagrangianAerosols.o \
	   $(OBJDIR)Condensation.o $(OBJDIR)soa_unidriv.o \
           $(OBJDIR)soa_unidriva.o $(OBJDIR)unifac.o \
	   $(OBJDIR)OutputRoutines.o \
	   $(OBJDIR)Coagulation.o \
	   $(OBJDIR)ReadInputDeck.o \
	   $(OBJDIR)StepASP4SAM.o \
           $(OBJDIR)EulerianBoxDriver.o
OBJECTS := $(OBJECTS) \
           $(OBJDIR)opkda2.o \
           $(OBJDIR)opkda1.o \
           $(OBJDIR)opkdmain.o

EXEOBJECTS := $(OBJECTS) $(OBJDIR)DevelopmentDriver.o	

# === Modules: 0 grid define; 1 essential; d type define; s all
MODASP := $(OBJDIR)DMiLay.o $(OBJDIR)MIEV0.o $(OBJDIR)MieSupport.o $(OBJDIR)CDFNor.o $(OBJDIR)indexx.o $(OBJDIR)lsodes.o \
           $(OBJDIR)PPND16.o \
	   $(OBJDIR)dgaus8.o \
	   $(OBJDIR)ModelParameters.o \
	   $(OBJDIR)InfrastructuralCode.o \
	   $(OBJDIR)GridPointFields.o \
	   $(OBJDIR)Time.o \
	   $(OBJDIR)Chemistry.o \
	   $(OBJDIR)Gibbs-DuhemEvaluationFunction.o \
           $(OBJDIR)LagrangianAerosols.o \
	   $(OBJDIR)Condensation.o $(OBJDIR)soa_unidriv.o \
           $(OBJDIR)soa_unidriva.o $(OBJDIR)unifac.o \
	   $(OBJDIR)OutputRoutines.o \
	   $(OBJDIR)Coagulation.o \
	   $(OBJDIR)ReadInputDeck.o \
	   $(OBJDIR)StepASP4SAM.o \
           $(OBJDIR)EulerianBoxDriver.o

MODASP := $(MODASP) \
          $(OBJDIR)opkda2.o \
          $(OBJDIR)opkda1.o \
          $(OBJDIR)opkdmain.o

MODASPEXEC := $(MODASP) $(OBJDIR)DevelopmentDriver.o

MODS := $(MODASP)

EXECMODS := $(MODASPEXEC)

PROG = asp_v2.1.1_linux_ifort_dbl

##
## === Dependency
##

#ALL: $(EXEDIR)$(PROG)

#$(EXEDIR)$(PROG): $(EXEOBJECTS)
#	$(FC) $(FFLAGS) -o $@ $(EXEOBJECTS) $(EXECMODS) $(LDFLAGS)




OBJECTS : $(MODS) $(FILESINC)

$(OBJDIR)opkda1.o: $(COMPDIR)opkda1.f
	$(FC) -c $(FFLAGS) $<
	  \mv opkda1.o $(OBJDIR)
$(OBJDIR)opkda2.o: $(COMPDIR)opkda2.f
	$(FC) -c $(FFLAGS) $<
	  \mv opkda2.o $(OBJDIR)
$(OBJDIR)opkdmain.o: $(COMPDIR)opkdmain.f $(OBJDIR)opkda1.o $(OBJDIR)opkda2.o
	$(FC) -c $(FFLAGS) $<
	  \mv opkdmain.o $(OBJDIR)

$(OBJDIR)CDFNor.o: $(COMPDIR)CDFNor.f90	
	  $(FC) -c $(F90FLAGS) $<
	  \mv CDFNor.o $(OBJDIR)
$(OBJDIR)indexx.o: $(COMPDIR)indexx.f
	  $(FC) -c $(FFLAGS) $<
	  \mv indexx.o $(OBJDIR)
$(OBJDIR)lsodes.o: $(COMPDIR)lsodes.f
	  $(FC) -c $(FFLAGS) $<
	  \mv lsodes.o $(OBJDIR)
$(OBJDIR)PPND16.o: $(COMPDIR)PPND16.f
	  $(FC) -c $(FFLAGS) $<
	  \mv PPND16.o $(OBJDIR)
$(OBJDIR)dgaus8.o: $(COMPDIR)dgaus8.f
	  $(FC) -c $(FFLAGS) $<
	  \mv dgaus8.o $(OBJDIR)
$(OBJDIR)MIEV0.o: $(SOURCEDIR)MIEV0.f
	  $(FC) -c $(FFLAGS) $<
	  \mv MIEV0.o $(OBJDIR)
$(OBJDIR)MieSupport.o: $(SOURCEDIR)MieSupport.f
	  $(FC) -c $(FFLAGS) $<
	  \mv MieSupport.o $(OBJDIR)
$(OBJDIR)DMiLay.o: $(SOURCEDIR)DMiLay.f
	  $(FC) -c $(FFLAGS) $<
	  \mv DMiLay.o $(OBJDIR)
$(OBJDIR)unifac.o: $(UNIFACDIR)unifac.f
	  $(FC) -c $(FFLAGS) $<
	  \mv unifac.o $(OBJDIR)
$(OBJDIR)soa_unidriv.o: $(UNIFACDIR)soa_unidriv.f \
	$(OBJDIR)unifac.o $(UNIFACDIR)unifacparam.h
	  $(FC) -c $(FFLAGS) $<
	  \mv soa_unidriv.o $(OBJDIR)
$(OBJDIR)soa_unidriva.o: $(UNIFACDIR)soa_unidriva.f \
	$(OBJDIR)unifac.o $(UNIFACDIR)unifacparamMattAq.h
	  $(FC) -c $(FFLAGS) $<
	  \mv soa_unidriva.o $(OBJDIR)
$(OBJDIR)ModelParameters.o: $(SOURCEDIR)ModelParameters.f90 
	  $(FC) -c $(F90FLAGS) $<
	  \mv ModelParameters.o $(OBJDIR)
$(OBJDIR)Gibbs-DuhemEvaluationFunction.o: \
	$(COMPDIR)Gibbs-DuhemEvaluationFunction.f90
	  $(FC) -c $(F90FLAGS) $<
	  \mv Gibbs-DuhemEvaluationFunction.o $(OBJDIR)
$(OBJDIR)InfrastructuralCode.o: $(INFDIR)InfrastructuralCode.f90 \
	$(INFINC)
	  $(FC) -c $(F90FLAGS) $<
	  \mv InfrastructuralCode.o $(OBJDIR)
$(OBJDIR)GridPointFields.o: $(EULERDIR)GridPointFields.f90 \
	$(OBJDIR)ModelParameters.o \
	$(OBJDIR)InfrastructuralCode.o $(EULERINC) $(INFINC)
	  $(FC) -c $(F90FLAGS) $<
	  \mv GridPointFields.o $(OBJDIR)       	
$(OBJDIR)Time.o: $(SOURCEDIR)Time.f90 $(OBJDIR)ModelParameters.o       	
	  $(FC) -c $(F90FLAGS) $<
	  \mv Time.o $(OBJDIR)   
$(OBJDIR)Chemistry.o: $(CHEMDIR)Chemistry.f90 \
       $(OBJDIR)Time.o 	$(OBJDIR)ModelParameters.o \
       $(OBJDIR)InfrastructuralCode.o \
       $(OBJDIR)GridPointFields.o $(CHEMINC) $(EULERINC) $(INFINC)
	  $(FC) -c $(F90FLAGS) $<
	  \mv Chemistry.o $(OBJDIR)   
$(OBJDIR)ReadInputDeck.o: $(SOURCEDIR)ReadInputDeck.f90 \
	$(OBJDIR)GridPointFields.o $(OBJDIR)Chemistry.o \
	$(OBJDIR)ModelParameters.o\
	$(OBJDIR)InfrastructuralCode.o $(CHEMINC) \
	$(EULERINC) $(INFINC)
	  $(FC) -c $(F90FLAGS) $<
	  \mv ReadInputDeck.o $(OBJDIR)   
$(OBJDIR)LagrangianAerosols.o: $(PARTDIR)LagrangianAerosols.f90 \
	$(OBJDIR)GridPointFields.o   $(OBJDIR)Chemistry.o   \
	$(OBJDIR)ModelParameters.o \
	$(OBJDIR)InfrastructuralCode.o   $(OBJDIR)Time.o \
	$(CHEMINC) $(EULERINC) $(INFINC) $(AEROINC) $(PARTINC)
	  $(FC) -c $(F90FLAGS) $< 
	  \mv LagrangianAerosols.o $(OBJDIR)  
$(OBJDIR)Condensation.o: $(CONDDIR)Condensation.f90  \
	$(OBJDIR)GridPointFields.o $(OBJDIR)Chemistry.o \
	$(OBJDIR)ModelParameters.o \
	$(OBJDIR)InfrastructuralCode.o \
	$(OBJDIR)LagrangianAerosols.o $(CHEMINC) \
	$(EULERINC) $(INFINC) $(AEROINC) $(CONDINC)
	  $(FC) -c $(F90FLAGS) $<
	  \mv Condensation.o $(OBJDIR)   	
$(OBJDIR)OutputRoutines.o: $(OUTPUTDIR)OutputRoutines.f90 \
	$(OBJDIR)GridPointFields.o   $(OBJDIR)Chemistry.o \
	$(OBJDIR)ModelParameters.o \
	$(OBJDIR)InfrastructuralCode.o   $(OBJDIR)LagrangianAerosols.o \
	$(OBJDIR)Condensation.o $(CHEMINC) \
	$(EULERINC) $(INFINC) $(AEROINC) $(CONDINC)
	  $(FC) -c $(F90FLAGS) $<  
	  \mv OutputRoutines.o $(OBJDIR) 
$(OBJDIR)Coagulation.o: $(COAGDIR)Coagulation.f90  \
	$(OBJDIR)GridPointFields.o   $(OBJDIR)Chemistry.o \
     	$(OBJDIR)ModelParameters.o \
	$(OBJDIR)InfrastructuralCode.o $(OBJDIR)LagrangianAerosols.o \
	$(OBJDIR)OutputRoutines.o $(CHEMINC) \
	$(EULERINC) $(INFINC) $(AEROINC) $(CONDINC) $(COAGINC) 
	  $(FC) -c $(F90FLAGS) $<  
	  \mv Coagulation.o $(OBJDIR) 
$(OBJDIR)StepASP4SAM.o: $(SOURCEDIR)StepASP4SAM.f90 \
	$(OBJDIR)GridPointFields.o   $(OBJDIR)Condensation.o \
       	$(OBJDIR)ModelParameters.o \
	$(OBJDIR)Time.o $(OBJDIR)LagrangianAerosols.o $(ALLINC)
	  $(FC) -c $(F90FLAGS) $<  
	  \mv StepASP4SAM.o $(OBJDIR) 
$(OBJDIR)EulerianBoxDriver.o: $(SOURCEDIR)EulerianBoxDriver.f90 \
	$(OBJDIR)GridPointFields.o   $(OBJDIR)Condensation.o \
       	$(OBJDIR)ModelParameters.o \
	$(OBJDIR)Time.o $(OBJDIR)LagrangianAerosols.o $(ALLINC)
	  $(FC) -c $(F90FLAGS) $<  
	  \mv EulerianBoxDriver.o $(OBJDIR) 
$(OBJDIR)DevelopmentDriver.o: $(SOURCEDIR)DevelopmentDriver.f90 \
	$(OBJECTS) $(ALLINC)
	  $(FC) -c $(F90FLAGS) $< 
	  \mv DevelopmentDriver.o $(OBJDIR)  
#	  \mv *.mod $(OBJDIR) 
#
# === Rules
#

.SUFFIXES: .mod .F .F90 .f90 .f

.F90.mod:
	$(FC) -c $(F90FLAGS) $<

.F.mod:
	$(FC) -c $(FFLAGS) $<

.F90.o:
	$(FC) -c $(F90FLAGS) $<

.F.o:
	$(FC) -c $(FFLAGS) $<

.f90.mod:
	$(FC) -c $(F90FLAGS) $<

.f.mod:
	$(FC) -c $(FFLAGS) $<

.f90.o:
	$(FC) -c $(F90FLAGS) $<

.f.o:
	$(FC) -c $(FFLAGS) $<

##
## === Functions
##

clean:
	rm -f  $(OBJDIR)*.o $(MODULES_OUTPUT_DIR)*.mod $(OBJDIR)*.il $(OBJDIR)*.d $(OBJDIR)*.pc* ./*.mod

# TODO: Original version does not have the F90 flags but instead FFLAGS; will
# 	this matter?  STILT-ASP will not call this make anyways...
# FIXME: Nothing re: standalone asp has been tested or attempted here
run_standalone_asp:
	$(MAKE) all
	$(FC) $(F90FLAGS) -o $(ASP_STANDALONE) $(EXEOBJECTS) $(EXEMODS)

# Will the second one work on ifort?  Probably not...
# Also, add asp_wrapper.o
#$(ASP_LIB): $(OBJECTS) $(MODS) asp_wrapper.o
#	$(AR) $(AR_FLAGS) $(ASP_LIB) asp_wrapper.o $(OBJECTS) $(MODS)

all: $(OBJECTS) $(MODS)
	$(AR) $(AR_FLAGS) $(ASP_LIB) $(OBJECTS) $(MODS)

install:
	echo ""; \
	if [ -a $(ASP_LIB) ]; then \
		echo "Moving ASP Library '$(ASP_LIB)' to $(LIBRARY_INSTALL_DIR)"; \
		mv $(ASP_LIB) $(LIBRARY_INSTALL_DIR); \
		echo "ASP library '$(ASP_LIB)' successfully installed in $(LIBRARY_INSTALL_DIR)'"; \
	else \
		echo "*** Error: Cannot install ASP library '$(ASP_LIB)' ***";\
		echo "*** Did you forget to run a 'make all' first? ***";\
	fi;

