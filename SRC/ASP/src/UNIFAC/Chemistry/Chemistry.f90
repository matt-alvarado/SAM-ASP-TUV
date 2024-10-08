!! ASP (c), 2004-2016, Matt Alvarado (malvarad@aer.com)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! Chemistry.f90
!! This is the main source file for chemical composition of the gas
!! and aerosol phases. It also contains the numerical integration 
!! routine for gas phase chemistry.
!! This is the source file linking many header files together.
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY							     !!
!!							                     !!
!! Month  Year   Name              Description				     !!
!! 07     2006   Matt Alvarado     Began Update History			     !!
!! 10/19  2006   Matt Alvarado     Changed ResetPhotoRates to do linear	     !!
!!			             interpolation.			     !!
!! 02/16  2012   Matt Alvarado     Removed Eulerian grids, making ASP        !!
!!                                 a one-box model or subroutine.            !!
!! 02/27  2012   Matt Alvarado     Removed Refractive Index Data, now        !!
!!                                  hardcoded in SUBROUTINE                  !!
!!                                  ShellRefIndAndRad(InParticle) in file    !!
!!                                  ParticleAttributes.h                     !!
!! 02/29  2012   Matt Alvarado     Added flag to AqPhaseChem.in and          !!
!!                                  SetAqPhaseChemistry to indicate what     !!
!!                                  refractive index should be used for an   !!
!!                                  electrolyte.                             !!
!! 01/04  2013   Matt Alvarado     Changed ReadPhotoRates to allow 50 rxns   !!
!! 01/11  2013   Matt Alvarado     Changed ReadPhotoRates to allow 51 rxns   !!
!! 01/24  2013   Matt Alvarado     Added OrgActivityFlag and OrgActivityIndex!!
!! 01/25  2013   Matt Alvarado     Added GasPhasePeroxy array                !!
!! 01/25  2013   Matt Alvarado     Added functions GetTotalPeroxy and        !!
!!                                       GetTotalAcylPeroxy                  !!
!! 01/25  2013   Matt Alvarado    Changed RecalculateReactionRates to pass   !!
!!                                  total peroxy concentration.              !!
!! 01/24  2013   Matt Alvarado     Added AqOrgActivityFlag,AqOrgActivityIndex!!
!!                                    Kappa, and GammaInf                    !!
!! 07/03  2013   Matt Alvarado     Added abilty to interpolate between 2     !!
!!                                 arbitrary points to ResetPhotoRates       !!
!! 09/13  2016   Matt Alvarado     Added array GasPhaseEmisFactor to store   !!
!!                                   fire emission factors in g/kg           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DEPENDENCIES		                    !!
!! 1. ModelParameters			    !!
!! 2. InfrastructuralCode		    !!
!! 3. GridPointFields			    !!
!! 4. Time				    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! LINKED HEADER FILES						    !!
!! 1. ChemistryDerivativeStructures.h				    !!
!! 2. ChemistryParametersInitRoutines.h				    !!
!! 3. ChemicalPropertyRoutines.h				    !!
!! 4. AqPhaseReactionIntegration.h				    !!
!!		NOTE: Does not integrate! (07/2006)		    !!
!! 5. GasPhaseChemistryInits.h					    !!
!! 6. AqPhaseChemistryInits.h					    !!
!! 7. EvolveGasChemistry.h					    !!
!! 8. OrgPhaseChemInits.h"					    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:	             !!
!! 1. SUBROUTINE SetPhaseIfs(iIFGAS, iIFAQ)
!! 2. SUBROUTINE SetChemistryParams()
!! 3. SUBROUTINE RecalculateReactionRates(IFGAS, IFAQ, IFSOLID)
!! 4. SUBROUTINE RecalculateHeteroRates(Radii, NumConc, NumBins, FirstTime)
!! 5. SUBROUTINE ResetHeteroRates
!! 6. SUBROUTINE ReadPhotolysisRates
!! 7. SUBROUTINE ResetPhotoRates
!! 8. SUBROUTINE RecalculatePhotoRates
!! 9. FUNCTION GetAqCharge (ChemIndex)
!!10. FUNCTION FindChem(Chem, Phase, SuppressError) RESULT (Index)
!!11. FUNCTION FindIon(InChem, SuppressError) RESULT (IonIndex)
!!12. FUNCTION FindAnion(InChem)
!!13. FUNCTION FindCation(InChem)
!!14. FUNCTION GetTotalPeroxy
!!14. FUNCTION GetAcylPeroxy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE Chemistry

    ! CMB: add these two
	PRIVATE :: parcelIndex, timestepIndex
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Specify which variables, arrays, subroutines, and functions	!!
	!! and may be accessed from outside the module			!!
	PUBLIC	::      EnthalpyOfVaporization,		&
		        HowManyGasChems,		&	!Variables
			HowManyAqChems,			&	!	  "
			HowManyAqCations,		&	!     "
			HowManyAqAnions,		&	!     "
			HowManyOrgChems,                &   !
			HowManyAqOrgChems,              &
			HowManyAqEqReactions,		&	!	  "
			HowManyEvolveGasChems,		&	!     "
			HowManyKMParams,		&	!     "
			GasPhaseReactionRates,		&
			GasPhaseReactionRateArraySize,  &	
			SetChemistryParams,		&	! Subroutines
			SetPhaseIfs,			&	!     "
			ReadPhotolysisRates,		&
			GasPhaseChemicalNames,		&	! Names for Error Checking and Messaging
			AqPhaseChemicalNames,		&	!     "
			AqCationNames,AqAnionNames,	&	!     "
			OrgPhaseChemicalNames,          &
			OrgMolecularMass,               &
			OrgDensity,			&
                        OrgActivityFlag,                &
                        OrgActivityIndex,               &
			AqOrgPhaseChemicalNames,        &
			AqOrgMolecularMass,             &
			AqOrgDensity,			&
			SolidSaltDensity,               &
			GasMolecularMass,		&! Chemical Mass Values
			AqMolecularMass,		&	!     "
			AqCationMass, AqAnionMass,	&	!     "
			AqCationCharge,			&! For Thermodynamic Calculations
			AqAnionCharge,			&!     "  (Charge)
			GetAqCharge,			&
			AqEquilibriaList,		&!(Eq Reactions Index List)
			KMRef,				&
			SolidMolecularMass,		&
			IFGAS, IFAQ, IFSOLID,		&! Logical Variables
			RecalculateReactionRates,	&! 
			RecalculateHeteroRates,         &
			RecalculatePhotoRates,          &
			FindChem,			&! Indexing Functions & Variables
			FindIon,			&	!	  "
			FindCation,FindAnion,		&
			GasPhase,AqPhase,               & 
			StepGasChemistry,		&
			MolecularDiffusionCoefficient,  &
			GasPhaseBackground,             &
                        GasPhaseEmisFactor,              & !g/kg dry matter burned, MJA, 2016-09-13
			AqOrgNumCarbon,                 &
			GasPhaseDepVel,                 &
                        CationIonicRefraction,          &
                        AnionIonicRefraction,           &
                        RefIndFlag,                     &
                        AqOrgActivityFlag,              &
                        AqOrgActivityIndex,             &
                        Kappa,                          &
                        ! CB: Add a Photolysis array for tracer
                        GammaInf,                       &
                        !PhotoRateTimes,                 &
                        !PhotoRateZenithAngles,          &
                        DegToRad,                        &
                        RadToDeg,                        &
                        piDouble,                       &
                        RecalculatePhotoRatesSZA,       &
                        getSzaIndices,                  &
                        ResetPhotoRatesSZA,  &
                        
                        ! CMB (AER): Add these from the locals; initialize
                        !            first one in initializeASP 
                        GasPhasePeroxy, &
                        GetTotalPeroxy, &
                        GetTotalAcylPeroxy, &
                        getTimestepIndex, setTimestepIndex, getParcelIndex, setParcelIndex
                        
                        
                        
                        !GasPhaseAcylPeroxy !, !, &
                        !PhotoRates  ! CMB: let's try allowing access here        
                        !GammaInf   ! CB out
					!     "
	!! These are indexed by the GAS PHASE INDEX
	CHARACTER (len=15), ALLOCATABLE :: GasPhaseChemicalNames(:)	   !
	REAL*8,		    ALLOCATABLE :: EnthalpyOfVaporization(:) 
! Peroxy Radical Flag, 2 if should be included in RO2_T and RCO3_T sums (acyl peroxy),
!                         1 if should be included in RO2_T sum only, 0 if not
	REAL*8,		    ALLOCATABLE :: GasPhasePeroxy(:)   
! The Latent Heat of Evaporation of Each Species (J/g)
	REAL*8,		    ALLOCATABLE :: GasMolecularMass(:)		   
! Molecular Mass (g / mole)
	REAL*8,		    ALLOCATABLE :: GasPhaseReactionRates(:,:)  
! This holds reaction rate information
	REAL*8,		    ALLOCATABLE :: GasPhaseBackground(:)
	REAL*8,		    ALLOCATABLE :: GasPhaseEmisFactor(:) !g/kg dry matter burned, MJA, 2016-09-13
	REAL*8,		    ALLOCATABLE :: GasPhaseDepVel(:)		
	REAL*8				:: PhotoRates(96,51)

        ! CB add (this will have to be made configurable in the future!)
        REAL*8, PARAMETER :: piDouble = 3.1415926535897932384626433
        REAL*8, PARAMETER :: DegToRad = piDouble/180.0
        REAL*8, PARAMETER :: RadToDeg = 1.0/DegToRad
        REAL*8 :: PhotoRateTimes(96), PhotoRateZenithAngles(96)
        !REAL*8 :: piDouble, DegToRad, RadToDeg
        ! end CB add

	INTEGER :: GasConcInputTypeFlag

!!! --
!!! -- DECODING KEY FOR GasAerEquilibriaList (a,b)
!!! -- 
!!! --   (i,1) : Gas Phase Index of the Species
!!! --   (i,2) : Aqueous Phase Index of the Species
!!! --   (i,3) : Associated Aq Phase Equilibrium Reaction
!!! --   (i,4) : Equilibrium Coefficient at 298 K (either true or infinite)
!!! --		imply that it condenses directly into a dissociated form)
!!! --   (i,5): - Delta Ho0 / R T0   for Equilibrium Coefficient 
!!! --            Temperature Dependence
!!! --   (i,6): - Delta Cp0 / R      for Equilibrium Coefficient 
!!! --            Temperature Dependence
	REAL*8,		    ALLOCATABLE :: GasAerEquilibriaList(:,:)
	INTEGER	:: GasPhaseReactionRateArraySize(2)

	!! These are indexed by the AQUEOUS PHASE INDEX
	CHARACTER (len=15), ALLOCATABLE :: AqPhaseChemicalNames(:) 
	REAL*8,		    ALLOCATABLE :: SolidSaltDensity(:)  
	Integer,	    ALLOCATABLE ::  RefIndFlag(:)
! Electrolyte Density (mole / cm3) for aq. species
	REAL*8,		    ALLOCATABLE :: AqPhaseReactionRateInfo(:,:) 
! This holds reaction rate information
	REAL*8,		    ALLOCATABLE :: AqPhaseReactionRates(:,:)   
! This holds reaction rates v. temperature
	REAL*8,		    ALLOCATABLE :: AqMolecularMass(:)	
! Molecular Mass (g / mole)
	INTEGER	:: AqPhaseReactionRateArraySize(2)

	INTEGER, ALLOCATABLE :: AqChemicalCharge(:)			
	! The charge allocations are reportd as integral values

	CHARACTER (len=15),	 ALLOCATABLE :: AqCationNames(:)	! 
	CHARACTER (len=15),	 ALLOCATABLE :: AqAnionNames(:)		! 
	INTEGER, ALLOCATABLE :: AqCationCharge(:)		
	! The charge allocations are reportd as + integral values
	INTEGER, ALLOCATABLE :: AqAnionCharge(:)	
	! The charge allocations are reportd as + integral values
	REAL*8,  ALLOCATABLE :: AqIonGrid(:,:,:)		
	! Maps Ion Pair to Electrolyte and Provides Equilibrium Data
	REAL*8,  ALLOCATABLE :: AqCationMass(:)			
	! The charge allocations are reportd as + integral values
	REAL*8,  ALLOCATABLE :: AqAnionMass(:)			
	! The charge allocations are reportd as + integral values
	REAL*8,  ALLOCATABLE :: CationIonicRefraction(:)	
		!Ionic refraction for each ion	
	REAL*8,  ALLOCATABLE :: AnionIonicRefraction(:)		
	!Ionic refraction for each ion	
	
!!! -- ION-ELECTROLYTE MAP and EQUILIBRIUM DATA -- !!!
!!! --
!!! -- DECODING KEY FOR AqEquilibriaList (a,b)
!!! -- 
!!! --   (i,1) : Electrolyte  : this is -1 if water is the only dissociator 
!!! --           (e.g., H2O <=> H+ & OH-)
!!! --   (i,2) : Cation 
!!! --   (i,3) : Anion
!!! --   (i,4) : Cation Stoicheometric Coefficient
!!! --   (i,5) : Anion  Stoicheometric Coefficient
!!! --   (i,6) : Flag for solid compound
!!! --   (i,7) : Blank 
!!! --   (i,8) : Number of H2O Molecules Needed on LHS of Dissociation
!!! --   (i,9) : Equilibrium Coefficient at 298 K (either true or infinite)
!!! --   (i,10): - Delta Ho0 / R T0   for Equilibrium Coefficient Temp Dep
!!! --   (i,11): - Delta Cp0 / R      for Equilibrium Coefficient Temp Dep
!!! --   (i,12): First Density Parameter for dissociated ion solution
!!! --   (i,13): Second Density Parameter for dissociated ion solution
!!! --   (i,14): Surface Tension Surface Excess at Saturation (mol / m2)
!!! --   (i,15): Surface Tension Adsorption coefficient (dimentionless)
!!!
!!! ----------- Indexing of activities for various reaction types:
!!! --- for simple reactions where everything fits KM, (I,17) = (I,18) = 0 and
!!! ---  (I,16): Index into KMRef for thermodynamic reaction data
!!!
!!! --- for dissociating electrolytes (like bisulfate) (I,18) = 0 and
!!! ---  (I,16): Index into AqEquilibriaList of numerator reaction
!!! ---  (I,17): Index into KMRef of denominator reaction 
!!! ---          (which is this actual reaction)
!!!
!!! --- for complex reactions in which the activity is a 
!!! --    composite of several other reactions
!!! --   (i,16): Numerator Activity for Equilibrium Reaction 
!!! --           (stored as reaction+exponent/10,000)
!!! --   (i,17): Denominator Activity for Equilibrium Reaction 
!!! --           (stored as reaction+exponent/10,000)
!!! --   (i,18): Another Numerator Activity 
!!! --- for levitocite (NH4)3H(SO4)3, follow Kim et al., 1993 (I,18) = 0 and
!!! --   (i,16): Index into AqEquilibriaList for Ammonium Sulfate
!!! --   (i,17): Index into AqEquilibriaList for Sulfuric Acid 
!!!
!!! --   (i,19): A flag: if this is a positive number then this reaction 
!!! --           will be assigned the mixed coefficient of 
!!!              another reaction (for bisulfate mixing and Levitocite)
!!! --   (i,20): DRH at 298 K
!!! --   (i,21): Temperature dependence of DRH
!!! --   (i,22): Third dissociated species index (for levitocite)
!!! --   (i,23): Third dissociated species stoichiometry
!!! --   (i,24): Third dissociated species index (for hydrates, should be water
!!! --   (i,25): Third dissociated species stoichiometry

! The list of electrolytes used in the model, 
! and also density and surface tension parameters
	REAL*8,  ALLOCATABLE :: AqEquilibriaList(:,:)	

!!! -- DECODING KEY FOR KMRef (a,b)
!!! -- 
!!! --   (i,1) : Electrolyte
!!! --   (i,2) : Cation 
!!! --   (i,3) : Anion
!!! --   (i,4) : Cation Stoicheometric Coefficient
!!! --   (i,5) : Anion  Stoicheometric Coefficient
!!! --   (i,6) : Kusik-Meissner Coefficient
!!! --   (i,7) : Kusik-Meissner Temperature Dependence
!!! --   (i,8) : Number of H2O Molecules Needed on LHS of Dissociation

	REAL*8,  ALLOCATABLE :: KMRef(:,:)				
	INTEGER :: HowManyKMParams

	!! These are used in the Aqueous Phase Reaction Rate Table
	INTEGER, PARAMETER :: NumbAqPhaseRateTemps = 50
	INTEGER, PARAMETER :: MinAqPhaseRateTemp   = 250
	INTEGER, PARAMETER :: MaxAqPhaseRateTemp   = 315

	!! These are indexed by the ORGANIC PHASE INDEX
	CHARACTER (len=15), ALLOCATABLE :: OrgPhaseChemicalNames(:) 
	REAL*8,		ALLOCATABLE :: OrgMolecularMass(:)		   
        ! Molecular Mass (g / mole)
	REAL*8,		ALLOCATABLE :: OrgDensity(:)		   
        ! Density (g / cm3)
	REAL*8,		ALLOCATABLE :: OrgActivityFlag(:)
        REAL*8,		ALLOCATABLE :: OrgActivityIndex(:)


	!! These are indexed by the HYDROPHILIC ORGANIC COMPOUND INDEX
	CHARACTER (len=15), ALLOCATABLE :: AqOrgPhaseChemicalNames(:) 
	REAL*8,		ALLOCATABLE :: AqOrgMolecularMass(:)		   
        ! Molecular Mass (g / mole)
	REAL*8,		ALLOCATABLE :: AqOrgDensity(:)		   
        ! Density (g / cm3)
	INTEGER,	ALLOCATABLE :: AqOrgNumCarbon(:)
	REAL*8,		ALLOCATABLE :: Kappa(:)
	REAL*8,		ALLOCATABLE :: GammaInf(:)
	REAL*8,		ALLOCATABLE :: AqOrgActivityFlag(:)
        REAL*8,		ALLOCATABLE :: AqOrgActivityIndex(:)		   
	
        !! Each array is sized by the number of evolving chemicals,
	!! but only a certain subset of these will evolve (placed first
	!! in each of the arrays).
	INTEGER :: HowManyGasChems,   HowManyEvolveGasChems, &
	           HowManyAqChems,    HowManyEvolveAqChems,  &
		   HowManyAqCations,  HowManyAqAnions, HowManyAqEqReactions, &
	     HowManyCondensingChems, HowManyOrgChems, HowManyEvolveOrgChems, &
		   HowManyAqOrgChems, HowManyEvolveAqOrgChems

	INTEGER :: AqWaterDissociationReaction

	!! The structures used in the linked lists are defined here
	INCLUDE "ChemistryDerivativeStructures.h"

	!! The actual ODE Set and Jacobian 
        !! are set as arrays of these linked lists
	TYPE (DiffEqList), ALLOCATABLE :: GasPhaseODESet(:),AqPhaseODESet(:)
	TYPE (DiffEqList), ALLOCATABLE :: GasPhaseJacobian(:,:), &
                                           AqPhaseJacobian(:,:)
	TYPE (DiffEqList), ALLOCATABLE :: GasPhaseSensitivityMatrix(:,:),  &
                                           AqPhaseSensitivityMatrix(:,:)
	
	INTEGER	:: HowManyNonZeroGasJacTerms
	INTEGER, ALLOCATABLE	:: GasIa(:), GasJa(:)
	INTEGER	:: HowManyNonZeroAqJacTerms
	INTEGER, ALLOCATABLE    :: AqIa(:), AqJa(:)

	!! Phases are indexed for some functions
	INTEGER, PARAMETER :: GasPhase   = 0
	INTEGER, PARAMETER :: AqPhase	 = 1

	!! LSODES Related Values
	INTEGER :: LRW, GasLIW, YLEN

	!! The Phases to Include
	LOGICAL :: IFGAS, IFAQ

! CB add
        integer :: parcelIndex, timeIndex
! end CB

	CONTAINS
! cmb add	
 integer function getParcelIndex() result(parcInd)
        !integer :: parcInd
        parcInd = parcelIndex
    end function
    subroutine setParcelIndex(parcInd)
        integer, intent(in) :: parcInd
        parcelIndex = parcInd
    end subroutine 
    
    integer function getTimestepIndex() result(timeInd)
        !integer :: timeInd
        timeInd = timestepIndex
    end function
    subroutine setTimestepIndex(timeInd)
        integer, intent(in) :: timeInd
        timestepIndex = timeInd
    end subroutine 
    ! end CMB
    
	SUBROUTINE SetPhaseIfs(iIFGAS, iIFAQ)

		!! Which Kinds of Chemistry Will Be Included?
		LOGICAL :: iIFGAS, iIFAQ

		!! Set Module Variables with those passed from Driver
		IFGAS   = iIFGAS
		IFAQ    = iIFAQ

	END SUBROUTINE SetPhaseIfs

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Publicly available routine to call several Chemical	!!
	!! Parameter initialization subroutines			!!
	!! This must be called after the grid point fields are  !!
	!! initialized.						!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE SetChemistryParams()

		IMPLICIT NONE

		!! Deal with Gas Phase Initialization
        !        print *, 'CMB: Setting GasPhaseChemistry'
                !call flush(6)
		CALL SetGasPhaseChemistry ()

		!! Deal with Aq Phase Initialization
        !        print *, 'CMB: Setting AqPhaseChemistry'
                !call flush(6)
		CALL SetAqPhaseChemistry ()

		!! Deal with hydrophobic org. phase initialization
        !        print *, 'CMB: Setting organic phase chemistry'
                !call flush(6)
		CALL SetOrgPhaseChemistry ()

		!! Deal with hydrophilic org. compound initialization
        !        print *, 'CMB: SEtting hydrophilic org phase chemistry'
                !call flush(6)
		CALL SetHydrophilicOrgPhaseChemistry ()

		!!!!!!!!!!!!!!!
		!! GAS PHASE !!
		IF (IFGAS) THEN

			!! Set the Gas Phase Chemistry ODE Set and Jacobian
            !            print *, 'CMB: setting gas phase ODE stuff'
                        !call flush(6)
			CALL SetGasPhaseODEandJacobian ()

			!! This Routine Sets the LSODES-related code
            !            print *, 'CMB: setting evolvegaschemconstants'
                        !call flush(6)
			CALL SetEvolveGasChemConstants ()

		END IF

		!!!!!!!!!!!!!!!!!!!
		!! AQUEOUS PHASE !!
		IF (IFAQ) then
            !            print *, 'CMB: Setting aqueous phase ODE set, etc'
            !            call flush(6)
			!! Set the Aqueous Phase Chemistry ODE Set and Jacobian
			CALL SetAqPhaseODEandJacobian ()
                endif

		!! Prepare the Dissociation Reactions
		!! (this always happens)
		!print *, 'CMB: Setting aq phase equil'
		!call flush(6)
		CALL SetAqPhaseEquilibriumReactions ()
		
		RETURN
	END SUBROUTINE SetChemistryParams
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! All of the routines that are used to initialize chemistry!!
	!! are contained within a .h file to be stored in the same  !!
	!! directory as this primary .f90 file.			    !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	INCLUDE "ChemistryParametersInitRoutines.h"
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! These routines calculate properties of particular chems, !!
	!! such as diffusivity and the like.			    !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	INCLUDE "ChemicalPropertyRoutines.h"

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Aqueous Phase Chemistry is Integrated Using Parameterized !!
	!! (against Temperature) Reaction Rates.  All routines are   !!
	!! stored here.                                              !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	INCLUDE "AqPhaseReactionIntegration.h"

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! The initialization routines for each phase are contained !!
	!! in their own files.					    !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	INCLUDE "GasPhaseChemistryInits.h"
	INCLUDE "AqPhaseChemistryInits.h"
	INCLUDE "EvolveGasChemistry.h"
	INCLUDE "OrgPhaseChemInits.h"

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! The reaction rates are recalculated sometime during		!!
	!! each timestep.  Gas phase reaction rates are stored at	!!
	!! each gridpoint and aqueous and solid are stored in each	!!
	!! of the particles.						!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE RecalculateReactionRates(IFGAS, IFAQ)

		USE Time, ONLY : CurrentTime, BeginningTime
		USE GridPointFields, ONLY : ResetGasPhaseReactionRates, &
                                            GridGasChem
		implicit none

		!! Input Variables
		LOGICAL :: IFGAS,IFAQ
  
                !!Internal Variables
                REAL*8 :: RO2_T, RCO3_T, C_H2O

                
		!!!!!!!!!!!!!!!
		!! Gas Phase !!
		IF (IFGAS) THEN
		    !write(678,*) 'In IFGAS() in RecalculateReactionRates: IFGAS = true, calling getters...'
                   RO2_T = GetTotalPeroxy() !molecules/cm3, Also resets RO2_T in chem array
                   RCO3_T = GetTotalAcylPeroxy() !molecules/cm3,Also resets RCO3_T in chem array
                   C_H2O = GridGasChem(1) !molecules/cm3
                   
                   !write(678,*) 'RO2_T = ', RO2_T
                   !write(678,*) 'RCO3_T = ', RCO3_T
                   !write(678,*) 'C_H2O = ', C_H2O
                   !write(678,*) 'In IFGAS: Calling ResetGasPhaseReactionRates'
                   CALL ResetGasPhaseReactionRates (GasPhaseReactionRates,GasPhaseReactionRateArraySize, RO2_T, RCO3_T, C_H2O)
		END IF	
		
		!!!!!!!!!!!!!!!!!!!
		!! Aqueous Phase !!
		IF (IFAQ) THEN
		    !print *, 'In IFAQ() in RecalculateReactionRates: IFAQ = true'
			IF (CurrentTime .EQ. BeginningTime) &
			CALL SetAqPhaseReactionRates (AqPhaseReactionRateInfo,AqPhaseReactionRateArraySize)
		END IF

	END SUBROUTINE RecalculateReactionRates

	SUBROUTINE RecalculateHeteroRates(Radii, NumConc, NumBins, FirstTime)
		!This subroutine just calls ResetHeteroRates

		USE Time, ONLY : CurrentTime, BeginningTime
		implicit none

		!! Input Variables
		REAL*8 :: Radii(:), NumConc(:)
		INTEGER  :: NumBins
		LOGICAL  :: FirstTime

		CALL ResetHeteroRates (GasPhaseReactionRates,GasPhaseReactionRateArraySize,Radii,NumConc,NumBins, FirstTime)
	
	END SUBROUTINE RecalculateHeteroRates

	SUBROUTINE ResetHeteroRates (ReactionRates,ArraySize,Radius,NumberConc,NumberBins, FirstTime)

		Use ModelParameters
		USE InfrastructuralCode, ONLY : Transcript, ERROR, Int2Str, & 
                                                Real2Str
		USE GridPointFields,	 ONLY :	GasPhaseChemicalRates, &
						GetTemp, &
						GetAirDensity
					  
		IMPLICIT NONE

		REAL*8,  INTENT(in) :: ReactionRates(:,:)
		REAL*8,  INTENT(in) :: Radius(:) !Radius in m
		REAL*8,  INTENT(in) :: NumberConc(:) 
                !Number Concentration in cm**(-3)
		
                INTEGER, INTENT(in) :: ArraySize(2)
		INTEGER, INTENT(in) :: NumberBins  
                !Number of particle types to sum over

		LOGICAL, INTENT(in) :: FirstTime

		INTEGER :: I,J,K,L,N
		REAL*8	:: Temp,Mgas,X(3),Uptake, MW, MWingrams, & 
                           Dg, R, c_bar, &
			   Surf_Conc, Term, Total, k_het, massofmolecule, A

		!! Loop over all reactions
		DO L = 1, ArraySize(1)
						
			!! Select only 1st order hetero rxns
			IF (INT(ReactionRates(L,1)) .EQ. 1 .AND. &
                            INT(ReactionRates(L,3)) .EQ. 3) THEN	
						
				!Assign uptake coefficient and molecular weight
				Uptake = ReactionRates(L,4)
				MW = ReactionRates(L,5)/1000 !kg/mol
				MWingrams = ReactionRates(L,5)
							
				!Assign Temperature
				TEMP = GetTemp()
						
		                !Calculate Diffusivity from MW
			       Dg = 0.375/(4.5*Angstrom)**2./GetAirDensity()* &
				SQRT(Rstar*Temp*AirMolecMass*ChemScale/	&
				Avogadro/MWingrams*(MWingrams+		&
				Avogadro*AirMolecMass)/2/Pi)
								
				Dg = Dg/1e4 !convert to m2/s from cm2/s	
	
				R = Rstar/1e7 !J/mol/K
			
				c_bar = SQRT(8.0*R*Temp/pi/MW) !m/s
								
				!Sum rate constants over all particle bins
				!See Jacob, 2000 Hetero. Chem. Review
				Total = 0.0
				DO N=1,NumberBins
				   Surf_Conc = 4.0*pi*(Radius(N)**2) &
                                                *(NumberConc(N)*1e6) !1/m
				   Term = Surf_Conc* ((Radius(N)/Dg + 4/Uptake/c_bar)**(-1.0)) !1/s
				   Total = Total + Term	
				   IF(FirstTime .AND. NumberBins .EQ. 1) THEN
					A = 100*Total/Surf_Conc !cm/s	
				   END IF
								
				END DO

				GasPhaseChemicalRates(L) = Total 						
				IF (FirstTime) then
                                    CALL TRANSCRIPT("Heterogeneous Reaction #"//TRIM(INT2STR(L))//&
                                                    " Rate Constant:"//TRIM(REAL2STR(GasPhaseChemicalRates(L)))//"")
                                endif

			!Second order heterogeneous reactions
			ELSE IF (INT(ReactionRates(L,1)) .EQ. 2 .AND. &
                                 INT(ReactionRates(L,3)) .EQ. 6) THEN	
				!Sum surface area over all particle bins
				Total = 0.0
				DO N=1,NumberBins
					Surf_Conc = 4.0*pi*(Radius(N)**2)*(NumberConc(N)*1e6) !1/m
					Total = Total + Surf_Conc	
				END DO
				
				GasPhaseChemicalRates(L) = (Total/100)*ReactionRates(L,4)						
				IF (FirstTime) then
                                    CALL TRANSCRIPT("2nd order Heterogeneous Reaction #"//TRIM(INT2STR(L))//&
                                                    " Rate Constant:"//TRIM(REAL2STR(GasPhaseChemicalRates(L)))//"")
                                endif
			END IF						
		END DO
	END SUBROUTINE ResetHeteroRates

	
	SUBROUTINE ReadPhotolysisRates()
	!!Matt Alvarado designed this subroutine to read in a text file 
        !!of photolysis
	!!rate data generated off-line by TUV v4.6 and store it in the 
        !! array PhotoRates(I,J)
	!! I = INT(UT*4 + 1), where UT is universal time in hours
	!! J = reaction rate index
	!! The file must contain 96 lines of data (every quarter hour)
	!! and 51 columns, with the first being the time in UT
        !! the second being solar zenith angle (sza)
	!! and the next 51 being the reaction in the order given in the header 
        !! of PhotoRates.in
	!! The lines must start at UT = 0.0
		
		USE InfrastructuralCode
		USE ModelParameters
		
		IMPLICIT NONE

		INTEGER :: FH, I, J
		CHARACTER (len=1024):: CurrLine, Time, Token, Sza


		FH = GetFileHandle()
            ! CB: separated into 2 lines for VIM
	    OPEN(UNIT=FH, FILE=TRIM(InputDeckSubDir)//"PhotoRates.in", &
                 STATUS='OLD')
		
		!The input file must have 96 lines and 51 reactions (53 column total, first one time, second solar zenith angle)
		
		DO I = 1, 96
			!! Get a new data line (non-comment) from the input deck
			CurrLine = GetLine(FH)
			!! Store the time
			CALL GetToken(CurrLine," ", Time)
			!! Store the sza
			CALL GetToken(CurrLine," ", Sza)

                        ! CB: Add these to the array
                        PhotoRateZenithAngles(i) = STR2REAL(Sza)
                        !print *, 'I = ', i
                        !print *, 'Sza = ', Sza(1:len(trim(Sza)))
                        !print *, 'PhotoRateZenithAngles = ', &
                        !         PhotoRateZenithAngles(i)
                        !print *, 'nextchar'
                        PhotoRateTimes(i) = STR2REAL(Time)
                        ! end CB add

			!Do for each photolysis reaction
			DO J = 1, 51
                                !write(*,*) I,J
				CALL GetToken(CurrLine," ", Token)
                                !write(*,*) Token
				Token = StripToken(Token)		
				PhotoRates(I,J) = STR2REAL(Token)
				!IF (J .EQ.25 ) WRITE(*,*) I, PhotoRates(I,J)
                                !CALL flush(6)
                                !STOP
			END DO
		END DO

	END SUBROUTINE ReadPhotolysisRates
        
        ! CB: I have added a routine to 
        SUBROUTINE getSzaIndices(sza, previous, next_indx)
        implicit none

        real*8, intent(in) :: sza
        integer, intent(out) :: previous, next_indx
        integer :: i
        real*8 :: current_zenith_angle, previous_angle
        real*8 :: next_angle

        previous = -1
        previous_angle = -360
        next_angle = 360
        next_indx = -1

        !print *, '-------- Input Zenith Angle  = ', &
        !         sza, ' --------------------'
        do i = 1, 96
            current_zenith_angle = PhotoRateZenithAngles(i)
            !print *, 'Current Zenith Angle: ', current_zenith_angle
            
            if (current_zenith_angle .le. sza) then
                if (current_zenith_angle .ge. previous_angle) then
                    !print *, '    resetting previous to ', i
                    previous_angle = current_zenith_angle
                    previous = i
                endif
            else !iif (current_zenith_angle .ge. sza) then
                if (current_zenith_angle .le. next_angle) then
                    !print *, '    resetting next to ', i
                    next_angle = current_zenith_angle
                    next_indx = i
                endif
!            else
!                print *, 'TODO: Equality case in Chemistry.', &
!                         'setSzaIndices(): sza of ', sza, ' is ', &
!                         'on a bin boundary...'
!                previous = i
!                next_indx = i
!            if ( (current_zenith_angle .lt. next_angle) .and. &
!                 (current_zenith_angle .gt. sza) ) then
!                next_angle = current_zenith_angle
!                next_indx = i
!            endif
!
!            if (current_zenith_angle < sza) then
!                print *, '    resetting previous to ', i
!                previous = i
!            else if (current_zenith_angle > sza) then
!                print *, '    resetting next to ', i
!                next = i
!            else
!                print *, 'TODO: Equality case in Chemistry.', &
!                         'setSzaIndices(): sza of ', sza, ' is ', &
!                         'on a bin boundary...'
!                previous = i
!                next = i
            endif
        enddo

        ! handle wraparound (this gets handled when we interpolate in
        ! cosine of sza later on I think)
        if (previous < 1) then
            previous = 96
        endif
        if (next_indx < 1) then
            next_indx = 1
        endif
    
        !stop ! let's fix this first
        end subroutine getSzaIndices

    SUBROUTINE ResetPhotoRatesSZA(ReactionRates, ArraySize, &
                                  FirstTime, sza, transmissivity)
        Use ModelParameters
        USE InfrastructuralCode, ONLY : Transcript, ERROR, & 
                                        Int2Str, Real2Str
        USE GridPointFields, ONLY : GasPhaseChemicalRates, &
                                    linear2dInterp

        IMPLICIT NONE

        REAL*8,  INTENT(in) :: ReactionRates(:,:)
        INTEGER, INTENT(in) :: ArraySize(2)
        !REAL*8,  INTENT(in) :: ut
        LOGICAL, INTENT(in) :: FirstTime
        
        INTEGER :: I,J,K,L, TimeIndex, RxnIndex
        REAL*8 :: OldRate, y1, y2, ut1, ut2, A, ut_temp
        REAL*8, PARAMETER :: Scale = 1.0

        ! cb add
        integer :: lowerZenithIndex, higherZenithIndex, ii
        real*8 :: DegHigher, DegLower, cosLower, cosHigher, szarad
        real*8 :: RateHigher, RateLower, InterpolatedRate, sza, transmissivity
        ! end cb add

        !! Loop over all reactions
        DO L = 1, ArraySize(1)
            !print *, ''
            !print *, 'L = ', L
            
            !do ii = 1, size(ReactionRates(L,:))
            !    print *, '    ReactionRates(', L, ', ', ii, ') = ', ReactionRates(L, ii)
            !enddo
                                    
            !! Only do photolysis reactions
            IF (INT(ReactionRates(L,1)) .EQ. 1 .AND. &
                INT(ReactionRates(L,3)) .EQ. 0.) THEN
                                                    
                OldRate = GasPhaseChemicalRates(L)
                !print *, 'OldRate/GasPhaseChemicalRates(', L, ') = ', oldrate                                    
                !Find reaction index number
                RxnIndex = INT(ReactionRates(L,4))
                !print *, 'rxnindex: ', rxnindex
!!Option 1 - if PhotoFlag = 0, find corresponding time in PhotoRates.in array
                IF (.not. PhotoFlag) THEN
! CB: This is not defined yet for mapping
!
                    print *, 'This option (PhotoFlag = false) ', &
                             'is currently not supported in ', &
                             'ResetPhotoRatesSZA.  Exiting...'
                    stop
!                    !Find lower time index number
!                    TimeIndex = INT(ut*4 + 1.0)
!
!                    ut1 = (TimeIndex - 1)/4.
!                    IF(TimeIndex .EQ. 96) THEN
!                        ut2 = 0.
!                    ELSE
!                        ut2 = (TimeIndex)/4.
!                    END IF
!
!                    !Early Value
!                    y1 = PhotoRates(TimeIndex, RxnIndex) 
!                    !Later Value
!                    IF(TimeIndex .EQ. 96) THEN
!                     y2 = PhotoRates(1, RxnIndex)
!                    ELSE
!                     y2 = PhotoRates(TimeIndex+1, RxnIndex)
!                    END IF
!
!                    A = (ut2 - ut)/(ut2 - ut1)
!                    !WRITE(*,*) A
!                                            
!                    GasPhaseChemicalRates(L) = Scale*(A*y1 + (1-A)*y2)

!!Option 2 - if PhotoFlag = 1, interpolate between starting and ending times of photolysis as given in the input deck
!! CB: This modification interpolates in the cosine of solar zenith
!angle
                ELSE

                    ! Find lower zenith angle index number
                    !print *, 'sza = ', sza
                    call getSzaIndices(sza, lowerZenithIndex, &
                                       higherZenithIndex)
                    !print *, 'Lower SZA index: ', lowerZenithIndex
                    !print *, 'Higher SZA index: ', higherZenithIndex
                    DegHigher=PhotoRateZenithAngles(higherZenithIndex)
                    DegLower=PhotoRateZenithAngles(lowerZenithIndex)

                    !print *, 'DegHigher = ', DegHigher
                    !print *, 'DegLower = ', DegLower
                    cosLower = dcos(DegLower*RadToDeg)
                    cosHigher = dcos(DegHigher*RadToDeg)
                    szarad = sza*DegToRad

                    RateLower = PhotoRates(lowerZenithIndex, RxnIndex)
                    RateHigher = PhotoRates(higherZenithIndex, &
                                            RxnIndex)

                    !print *, 'RateLower = ', RateLower
                    !print *, 'RateHigher = ', RateHigher

                    call linear2dInterp(cosLower, RateLower, &
                                        cosHigher, RateHigher, &
                                        dcos(szarad), InterpolatedRate)
                    !print *, 'Interpolated Rate = ', InterpolatedRate
                    !print *, 'Scale Factor = ', Scale
                    GasPhaseChemicalRates(L) = Scale*transmissivity*abs(InterpolatedRate)
                    !print *, 'GasPhaseChemicalRates(', l, ') = ', GasPhaseChemicalRates(L)
                    !!Find lower time index number
                    !ut1 = UTCStartTime
                    !TimeIndex = INT(ut1*4 + 1.0)
                    !!Early Value
                    !y1 = PhotoRates(TimeIndex, RxnIndex) 
                    !
                    !ut2 = PhotoEndTime
                    !TimeIndex = INT(ut2*4 + 1.0)
                    !!Later Value
                    !y2 = PhotoRates(TimeIndex, RxnIndex)
                    !
                    !!Calculate run end time in UTC
                    !ut2 = ut1+(RunLength/60.0)
                    !!Correct time for passing 24 hrs
                    !IF (ut .LT. ut1) THEN 
                    !      ut_temp = ut + 24.0
                    !ELSE
                    !      ut_temp = ut
                    !ENDIF
                    !A = (ut2 - ut_temp)/(ut2 - ut1)
                    !IF (L .EQ. 6) then
                    !    WRITE(*,*) ut, y1, y2, Scale*(A*y1 + (1-A)*y2)
                    !endif
                    !                        
                    !GasPhaseChemicalRates(L) = Scale*(A*y1 + (1-A)*y2)
                                                        
                ENDIF

                IF (FirstTime) then
                    CALL TRANSCRIPT("Photolysis Reaction #" & 
                                    //TRIM(INT2STR(L)) &
                                    //" Rate Constant:"// &
                                    TRIM(REAL2STR(&
                                    GasPhaseChemicalRates(L)))//"")
                endif
            END IF
        END DO

    end SUBROUTINE ResetPhotoRatesSZA

	SUBROUTINE ResetPhotoRates(ReactionRates,ArraySize, ut, FirstTime)
	!!This subroutine looks for the appropriate value of the 
        !! photolysis rate constant
	!!given the time in UT and the reaction rate index 
        !! (read from the gas chemical mechanism)	
	!! Changed on 10/19/06 to include linear interpolation - Matt Alvarado
		Use ModelParameters
		USE InfrastructuralCode, ONLY : Transcript, ERROR, & 
                                                Int2Str, Real2Str
		USE GridPointFields,	 ONLY : GasPhaseChemicalRates

		IMPLICIT NONE

		REAL*8,  INTENT(in) :: ReactionRates(:,:)
		INTEGER, INTENT(in) :: ArraySize(2)
		REAL*8,  INTENT(in) :: ut
		LOGICAL, INTENT(in) :: FirstTime
		
		INTEGER :: I,J,K,L, TimeIndex, RxnIndex
		REAL*8 :: OldRate, y1, y2, ut1, ut2, A, ut_temp
		REAL*8, PARAMETER :: Scale = 1.0


	        !! Loop over all reactions
		DO L = 1, ArraySize(1)
		    !print *, ''
		    !print *, 'L = ', L
					
			!! Only do photolysis reactions
			IF (INT(ReactionRates(L,1)) .EQ. 1 .AND. &
                            INT(ReactionRates(L,3)) .EQ. 0.) THEN	
								
				OldRate = GasPhaseChemicalRates(L)
				!print *, 'OldRate: ', oldrate				
				!Find reaction index number
				RxnIndex = INT(ReactionRates(L,4))
                !print *, 'RxnIndex: ', RxnIndex
!!Option 1 - if PhotoFlag = 0, find corresponding time in PhotoRates.in array
                IF (.not. PhotoFlag) THEN	
                							
				    !Find lower time index number
				    TimeIndex = INT(ut*4 + 1.0)

					ut1 = (TimeIndex - 1)/4.
					IF(TimeIndex .EQ. 96) THEN
	                                    ut2 = 0.
					ELSE
					    ut2 = (TimeIndex)/4.
					END IF

					!Early Value
					y1 = PhotoRates(TimeIndex, RxnIndex) 
					!Later Value
					IF(TimeIndex .EQ. 96) THEN
					 y2 = PhotoRates(1, RxnIndex)
					ELSE
					 y2 = PhotoRates(TimeIndex+1, RxnIndex)
					END IF

					A = (ut2 - ut)/(ut2 - ut1)
					!WRITE(*,*) A
								
					GasPhaseChemicalRates(L) = Scale*(A*y1 + (1-A)*y2)
!!Option 2 - if PhotoFlag = 1, interpolate between starting and ending times of photolysis as given in the input deck
                                ELSE	
				        !Find lower time index number
					ut1 = UTCStartTime	
			                TimeIndex = INT(ut1*4 + 1.0)
					!Early Value
					y1 = PhotoRates(TimeIndex, RxnIndex) 

					ut2 = PhotoEndTime
	                                TimeIndex = INT(ut2*4 + 1.0)
					!Later Value
					y2 = PhotoRates(TimeIndex, RxnIndex)

                                        !Calculate run end time in UTC
                                        ut2 = ut1+(RunLength/60.0)
                                        !Correct time for passing 24 hrs
                                        IF (ut .LT. ut1) THEN 
                                              ut_temp = ut + 24.0
                                        ELSE
                                              ut_temp = ut
                                        ENDIF
					A = (ut2 - ut_temp)/(ut2 - ut1)
					IF (L .EQ. 6) WRITE(*,*) ut, y1, y2, Scale*(A*y1 + (1-A)*y2)
								
					GasPhaseChemicalRates(L) = Scale*(A*y1 + (1-A)*y2)
								        
                                ENDIF

				IF (FirstTime) then
                                    CALL TRANSCRIPT("Photolysis Reaction #"//TRIM(INT2STR(L))//&
                                                    " Rate Constant:"//TRIM(REAL2STR(GasPhaseChemicalRates(L)))//"")
                                endif        
			END IF	
		END DO

	END SUBROUTINE ResetPhotoRates
	
	! CB: Make a version of this for SZA
	subroutine RecalculatePhotoRatesSZA(sza, transmissivity, FirstTime)
	    implicit none

	    real*8 :: sza, transmissivity
	    logical :: FirstTime

	    call ResetPhotoRatesSZA(GasPhaseReactionRates, &
	        GasPhaseReactionRateArraySize, &
	        FirstTime, sza, transmissivity)

	end subroutine RecalculatePhotoRatesSZA
        ! end CB add

	SUBROUTINE RecalculatePhotoRates(ut, FirstTime)
		
		!This subroutine just calls ResetPhotoRates, and must be called
		!after every time RecalculateGasPhaseReactionRates is called. 
		
		USE Time, ONLY : CurrentTime, BeginningTime
		implicit none

		!! Input Variables
		REAL*8   :: ut
		LOGICAL  :: FirstTime

		CALL ResetPhotoRates (GasPhaseReactionRates,GasPhaseReactionRateArraySize, ut, FirstTime)
	
	END SUBROUTINE RecalculatePhotoRates

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Get the charge of an aqueous species indexed between		 !!
	!! 1 and HowManyAqChems + HowManyAqCations + HowManyAqAnions !!
	INTEGER FUNCTION GetAqCharge (ChemIndex) !!!!!!!!!!!!!!!!!!!!!!

		IMPLICIT NONE
		
		INTEGER :: ChemIndex

		IF(ChemIndex .LE. HowManyAqChems) THEN
			GetAqCharge = 0
			RETURN
		
		ELSE IF (ChemIndex .LE. HowManyAqChems + HowManyAqCations) THEN
			GetAqCharge = AqCationCharge(ChemIndex - HowManyAqChems)
			RETURN
		ELSE 
			GetAqCharge = AqAnionCharge(ChemIndex - HowManyAqChems - HowManyAqCations)
			RETURN
		END IF

	END FUNCTION

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Find the index of a particular chemical	  !!
	!!   Phase = 0 for Gas, 1 for Aq, 2 for Org   !!
	!! A result value of 0 indicates an error.	  !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	FUNCTION FindChem(Chem, Phase, SuppressError) RESULT (Index)

		USE InfrastructuralCode, ONLY : ERROR, StripToken

		IMPLICIT NONE
		INTEGER				  :: i, Index, Phase
		CHARACTER*(*)		  :: Chem
		CHARACTER (len = 7)   :: PhaseStr
		LOGICAL, OPTIONAL     :: SuppressError
		LOGICAL				  :: LocalSuppressError

		IF (PRESENT(SuppressError)) THEN
			LocalSuppressError = SuppressError
		ELSE 
			LocalSuppressError = .FALSE.
		END IF

			Index = 0

			!! There are several possible phase
			SELECT CASE (Phase)

				CASE (GasPhase)
					DO i = 1, HowManyGasChems
						IF (Trim(Chem) .eq. Trim(GasPhaseChemicalNames(i))) THEN
							Index = i
							EXIT
						END IF
					END DO
				!! Aqueous Chems may be Ions as well.  Search all types.
				CASE (AqPhase)
					DO i = 1, HowManyAqChems
						IF (Trim(Chem) .eq. Trim(AqPhaseChemicalNames(i))) THEN
							Index = i
							EXIT
						END IF
					END DO
					DO i = 1, HowManyAqCations
						IF (Trim(Chem) .eq. Trim(AqCationNames(i))) THEN
							Index = i + HowManyAqChems
							EXIT
						END IF
					END DO
					DO i = 1, HowManyAqAnions
						IF (Trim(Chem) .eq. Trim(AqAnionNames(i))) THEN
							Index = i + HowManyAqChems + HowManyAqCations
							EXIT
						END IF
					END DO

				CASE (2) !Organic Phase
					DO i = 1, HowManyOrgChems
						IF (Trim(Chem) .eq. Trim(OrgPhaseChemicalNames(i))) THEN
							Index = i
							EXIT
						END IF
					END DO
				
				CASE (3) !Hydrophilic Organic Phase
					DO i = 1, HowManyAqOrgChems
						IF (Trim(Chem) .eq. Trim(AqOrgPhaseChemicalNames(i))) THEN
							Index = i
							EXIT
						END IF
					END DO
				
				CASE DEFAULT

					IF (.NOT.PRESENT(SuppressError) .OR. .NOT.SuppressError) THEN
						CALL ERROR("Function FindChem was passed a bad Phase Value with the"//			&
								   " Chemical name "//TRIM(Chem)//".  It then got confused and now"//	&
								   " is quitting...  Look into it.")
					ELSE
						Index = -1
					END IF
				END SELECT

			!! If the chemical isn't found, then barf.
			IF (Index .EQ. 0 .AND. .NOT.LocalSuppressError) THEN
				IF (Phase == 0) PhaseStr = "Gas"
				IF (Phase == 1) PhaseStr = "Aqueous"
				IF (Phase == 2) PhaseStr = "Organic"
				IF (Phase == 3) PhaseStr = "HydrophilicOrganic"
				CALL ERROR("Looked for Chemical "//Trim(Chem)//" in the "//PhaseStr//" Chemical Name Index and did not find it.")
			END IF

		RETURN
	END FUNCTION

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Find the index of a particular ION  !!
	!! Only Applied to Aqueous Phase Chems !!
	!! -- Returns a 2 element array with   !!
	!! -- 1 for Cation and -1 for anion in !!
	!! -- slot one and then then index in  !!
	!! -- slot two.						   !!
	!! -- Each is filled as 0 if no result !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	FUNCTION FindIon(InChem, SuppressError) RESULT (IonIndex)

		USE InfrastructuralCode, ONLY : ERROR, StripToken

		IMPLICIT NONE

		!! External Variables
		INTEGER		  :: IonIndex(2)
		CHARACTER*(*) :: InChem
		LOGICAL, OPTIONAL     :: SuppressError

		!! Internal Variables
		INTEGER	:: I
		CHARACTER (len=200) :: Chem
		LOGICAL             :: LocalSuppressError

		IF(.NOT.PRESENT(SuppressError)) THEN
			LocalSuppressError = .FALSE.
		ELSE
			LocalSuppressError = SuppressError
		END IF

		Chem = StripToken(InChem)

		!! Fill with default values
		IonIndex(1) = 0 ; IonIndex(2) = 0

		DO I = 1, HowManyAqAnions
			IF (Trim(Chem) .eq. Trim(AqAnionNames(i))) THEN
				IonIndex(1) = -1
				IonIndex(2) = I
				EXIT
			END IF
		END DO

		IF (IonIndex(1) .EQ. 0) THEN
			DO I = 1, HowManyAqCations
				IF (Trim(Chem) .eq. Trim(AqCationNames(i))) THEN
					IonIndex(1) = 1
					IonIndex(2) = I
					EXIT
				END IF
			END DO
		END IF

		!! If the chemical isn't found, then barf.
		IF (IonIndex(1) .EQ. 0 .AND. .NOT.LocalSuppressError) THEN
			CALL ERROR("Looked for Ion "//Trim(Chem)//" in the Aquoues Phase Ion Name Indices and did not find it.")
		END IF

	END FUNCTION FindIon

	INTEGER FUNCTION FindAnion(InChem)

		USE InfrastructuralCode, ONLY : ERROR, StripToken

		IMPLICIT NONE

		!! External Variables
		CHARACTER*(*) :: InChem

		!! Internal Variables
		INTEGER	:: I
		CHARACTER (len=200) :: Chem

		Chem = StripToken(InChem)
		FindAnion = 0

		DO i = 1, HowManyAqAnions
			IF (Trim(Chem) .eq. Trim(AqAnionNames(i))) THEN
				FindAnion = I
				EXIT
			END IF
		END DO

		!! If the chemical isn't found, then barf.
		IF (FindAnion .EQ. 0) THEN
			CALL ERROR("Looked for Anion "//Trim(Chem)//" in the Aquoues Phase Ion Name Indices and did not find it.")
		END IF

	END FUNCTION

	INTEGER FUNCTION FindCation(InChem)

		USE InfrastructuralCode, ONLY : ERROR, StripToken

		IMPLICIT NONE

		!! External Variables
		CHARACTER*(*) :: InChem

		!! Internal Variables
		INTEGER	:: I
		CHARACTER (len=200) :: Chem

		Chem = StripToken(InChem)
		FindCation = 0

		DO i = 1, HowManyAqAnions
			IF (Trim(Chem) .eq. Trim(AqCationNames(i))) THEN
				FindCation = I
				EXIT
			END IF
		END DO

		!! If the chemical isn't found, then barf.
		IF (FindCation .EQ. 0) THEN
			CALL ERROR("Looked for Cation "//Trim(Chem)//" in the Aquoues Phase Ion Name Indices and did not find it.")
		END IF

	END FUNCTION


        !! This function calculates the concentration (molecules/cm3)
        !! of total peroxy radicals (GasPhasePeroxy flags 1 or 2)
        !! and places it in RO2_T slot and returns as the value
        !! MJA, 01/25/2013
	REAL*8  FUNCTION GetTotalPeroxy()
               
                USE GridPointFields, ONLY : GridGasChem, ReplaceGridCellChemBurden

		IMPLICIT NONE           

                INTEGER :: q, ind
                REAL*8, ALLOCATABLE  :: Chem(:)
                
                ALLOCATE(CHEM(HowManyGasChems))
                
                GetTotalPeroxy = 0.0
		!Calculate peroxy sum
		!write(678,*)''
		DO q = 1, HowManyGasChems
		    !write(678,*) 'In GetTotalPeroxy(): q = ', q
		   IF(GasPhasePeroxy(q) .gt. 0.005) THEN
		       !write(678,*) '    In GetTotalPeroxy(): GasPhasePeroxy(', q, ') = gt 0.005 = ', &
		            !GasPhasePeroxy(q)
                       Chem(q) = GridGasChem(q)
                       !write(678,*)'    Chem(', q, ') = ', Chem(q)
                       GetTotalPeroxy = GetTotalPeroxy + Chem(q) !ppb
                       !write(678,*) '    GetTotalPeroxy is now = ', GetTotalPeroxy
                   ENDIF
	        END DO	
               
               !Place in RO2_T slot
               ind = FindChem("RO2_T", 0)
               !write(678,*)'In GetTotalPeroxy(): RO2_T index = ', ind
               !write(678,*)'In GetTotalPeroxy(): Assigning ', GetTotalPeroxy, &
               !     'to Chem(', ind, ')'
               Chem(ind) = GetTotalPeroxy
               !write(678,*)'In GetTotalPeroxy(): Going in to ReplaceGridCellChemBurden...'
               CALL ReplaceGridCellChemBurden (ind,Chem(ind)) 
               if (allocated(chem)) deallocate(chem)    ! CMB add
               !write(678,*)'In GetTotalPeroxy(): GridGasChem(', ind, ') = ', &
               !     GridGasChem(ind)

               !Return value
               RETURN
   
        END FUNCTION GetTotalPeroxy

        !! This function calculates the concentration (molecules/cm3)
        !! of total acyl peroxy radicals (GasPhasePeroxy flag 2)
        !! and places it in RCO3_T slot and returns as the value
        !! MJA, 01/25/2013
	REAL*8  FUNCTION GetTotalAcylPeroxy()

                USE GridPointFields, ONLY : GridGasChem, ReplaceGridCellChemBurden

		IMPLICIT NONE           

                INTEGER :: q, ind
                REAL*8, ALLOCATABLE  :: Chem(:)

                ALLOCATE(CHEM(HowManyGasChems))

                GetTotalAcylPeroxy = 0.0
		!Calculate peroxy sum
		!write(678,*)''
		DO q = 1, HowManyGasChems
		    !write(678,*) 'In GetTotalAcylPeroxy(): q = ', q
		   IF(GasPhasePeroxy(q) .gt. 1.005) THEN
		    !write(678,*) '    In GetTotalAcylPeroxy(): GasPhasePeroxy(', q, ') = gt 1.005 = ', &
            !        GasPhasePeroxy(q)
                       Chem(q) = GridGasChem(q)
                       !write(678,*)'    Chem(', q, ') = ', Chem(q)
                       GetTotalAcylPeroxy = GetTotalAcylPeroxy + Chem(q) !ppb
                       !write(678,*) '    GetTotalAcylPeroxy is now = ', GetTotalAcylPeroxy
                   ENDIF
	        END DO	

               !Place in RCO3_T slot
               ind = FindChem("RCO3_T", 0)
               !write(678,*)'In GetTotalAcylPeroxy(): RCO3_T index = ', ind
               !write(678,*)'In GetTotalAcylPeroxy(): Assigning ', GetTotalAcylPeroxy, &
               !     'to Chem(', ind, ')'
               Chem(ind) = GetTotalAcylPeroxy
               !write(678,*)'In GetTotalAcylPeroxy(): Going in to ReplaceGridCellChemBurden...'
               CALL ReplaceGridCellChemBurden (ind,Chem(ind)) 
               !write(678,*)'In GetTotalAcylPeroxy(): GridGasChem(', ind, ') = ', &
               !     GridGasChem(ind)
               
               
               if (allocated(Chem)) deallocate(chem)
                    
               !Return value
               RETURN
   
        END FUNCTION GetTotalAcylPeroxy
END MODULE Chemistry




