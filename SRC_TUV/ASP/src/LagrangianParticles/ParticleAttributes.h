!! ASP (c), 2004-2013, Matt Alvarado (malvarad@aer.com)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! ParticleAttributes.h
!! Contains functions for calculating particle properties 
!! (density, radius, mass
!! curvature correction, Knudsen number, Reynolds number, etc.)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY							     !!
!!									     !!
!! Month  Year   Name              Description				     !!
!! 07     2006   Matt Alvarado     Began Update History			     !!
!! 07/24  2006   Matt Alvarado     Removed "ISNAN" to fit pgf90		     !!
!! 09/07  2006   Matt Alvarado   Changed RecalculateRadius to give lower edge!!
!!				      when number conc. = 0.		     !!
!! 11/01  2006   Matt Alvarado     Created ShellRefIndAndRad
!! 05/18  2007   Matt Alvarado     Revised ShellRefIndAndRad
!!				   Created AerosolOptProp
!! 07/13  2007   Matt Alvarado     Added line to  ShellRefIndAndRad
!!	 			      to keep core rad below shell rad
!! 07/17  2007   Matt Alvarado     Added line to SurfaceTension
!!				      to keep value above that of pure Ethanol
!! 07/18  2007   Matt Alvarado     Removed calls to ParticleDensity from
!!				CurvatureCorrection and OrgCurvatureCorrection
!!			 Just use values in particle data structure instead
!!				   In ParticleDensity, if ion loop fails, force
!!				      density of ionic solution to 1.0 g/cm3
!! 07/26  2007   Matt Alvarado     Fixed AbsCoreRad calculation in 
!!                                 ShellRefIndAndRad
!! 09/07  2007   Matt Alvarado     Fixed shell refractive index calcultion
!! 09/17  2007   Matt Alvarado     Removed write statement from ParticleDensity
!! 10/01  2007   Matt Alvarado     Set ParticleDensity to skip bins 
!!                                 with 0 particles
!! 10/03  2007   Matt Alvarado     Set Denomsum limit in SurfaceTension to 
!!                                 1.0e40
!! 10/04  2007   Matt Alvarado     Protected ShellRefIndAndRad from near zero 
!!                                  number concentrations
!! 10/11  2007   Matt Alvarado     Set HypotheticalElectrolyteConcentrations 
!!                                 to warn of charge
!!                                 imbalance, but not stop.
!! 10/15  2007   Matt Alvarado     Fixed Assymetry parameter in 
!!                                 ShellRefIndAndRad
!! 08/30  2010   Matt Alvarado     Changed ShellRefIndandRad, 

!!                                 to allow 62 photolysis wavelengths
!! 02/16  2012   Matt Alvarado     Removed Eulerian grids, making ASP        !!
!!                                 a one-box model or subroutine.            !!
!! 02/27  2012   Matt Alvarado     Changed ShellRefIndandRad, 
!!                                  AerosolOpticalDepth, and AerosolOptProp 
!!                                  to allow 451 bins between 250 nm and 700 nm
!!                                 ShellRefIndandRad now used OPAC wavelength
!!                                  dependent complex refractive indices
!! 02/29  2012   Matt Alvarado     Changed ShellRefIndandRad to use different
!!                                  wavelength dependent refractive indices
!!                                  for different electrolytes based on flags
!!                                  in AqPhaseChems.in. Also added refractive
!!                                  index for dust.
!! 05/03  2012   Matt Alvarado     Changed AerosolOptProp to calculate
!!                                  backscattering coefficient.
!! 05/08  2012   Matt Alvarado     Added special routine for dry BC only particles
!!                                  to particle density.
!! 05/03  2012   Matt Alvarado     Fixed calculation of hemispherical
!!                                  backscattering coefficient in AerosolOptProp
!  08/17  2012   Matt Alvarado     Added capability for external mixtures to ShellRefIndAndRad
!! 11/08  2012   Matt Alvarado     Added Maxwell-Garnett mixing rule to ShellRefIndAndRad
!! 01/24  2013   Matt Alvarado     Updated OrgCurvatureCorrection to give 1.0
!!                                  if particle number = 0.0
!! 07/10  2013   Matt Alvarado     Made it so only RecalculateRadius resets 
!!                                  the particle%surfacetension by calling 
!!                                  SurfaceTension, and CurvatureCorrection
!!                                  uses the particle%surfacetension
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:		!!
!! 1. SUBROUTINE ResetAllOriginPressureLevels (OriginPressureLevel)
!! 2. FUNCTION ParticleKnudsen (InParticle)
!! 3. FUNCTION ParticleKnudsenForEnergy (InParticle)
!! 4. FUNCTION CunninghamSlipCorrectionFactor (InParticle)
!! 5. FUNCTION AerosolPH (InParticle)
!! 6. FUNCTION ParticleDensity (InParticle) !Sets radius and effective radius
!! 7. FUNCTION ParticleMass (InParticle)
!! 8. FUNCTION GetPatwardhandAndKumarWeightings (InParticle)
!! 9. FUNCTION HypotheticalElectrolyteConcentrations (InParticle)
!!10. FUNCTION TerminalVelocity (InParticle)
!!11. SUBROUTINE RecalculateRadius (InParticle) !Calls particle density
!!12. FUNCTION CurvatureCorrection (InParticle)
!!13. FUNCTION SurfaceTension (InParticle)
!!14. FUNCTION ReynoldsNumber (InParticle)
!!15. REAL*8 FUNCTION OrgCurvatureCorrection (InParticle)
!!16. REAL*8 FUNCTION MolalityAqCarbon(InParticle)
!!17. REAL*8 FUNCTION AqueousSolutionMass
!!18. SUBROUTINE ShellRefIndAndRad(InParticle)
!!19. REAL*8 FUNCTION AerosolOpticalDepth()
!!20. SUBROUTINE AerosolOptProp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
! Modified by CMB (AER): Gfortran complains when "Particle" is used both
!                        as a structure type declarator and a local variable
!                        name; I have replaced the variable names with
!                        "p_Particle"


SUBROUTINE ResetAllOriginPressureLevels (OriginPressureLevel)

	IMPLICIT NONE

	Type(Particle),POINTER :: p_Particle
	REAL    :: OriginPressureLevel

	p_Particle => Particles%First

	IF (ASSOCIATED(p_Particle)) p_Particle%OriginPressureLevel = OriginPressureLevel
		
	DO WHILE (ASSOCIATED(p_Particle%Next))
		p_Particle => p_Particle%Next
		IF (ASSOCIATED(p_Particle)) p_Particle%OriginPressureLevel = OriginPressureLevel
	END DO

	RETURN
END SUBROUTINE ResetAllOriginPressureLevels

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate the Knudsen Number of the Particle	!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION ParticleKnudsen (InParticle)

	USE GridPointFields, ONLY : MeanFreePathOfAir

	IMPLICIT NONE

	Type(Particle),POINTER :: InParticle

	ParticleKnudsen = MeanFreePathOfAir ()/ InParticle%EffectiveRadius
	
	!write(*,*) 'CRL: ParticleKnudsen = ',ParticleKnudsen
	!write(*,*) 'CRL: MeanFreePathOfAir = ',MeanFreePathOfAir ()
	!write(*,*) 'CRL: InParticle%EffectiveRadius= ',InParticle%EffectiveRadius
 	!write(*,*) ''

END FUNCTION ParticleKnudsen


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate the Knudsen Number Of Energy of the Particle !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION ParticleKnudsenForEnergy (InParticle)

	USE GridPointFields, ONLY : MolecularThermalDiffusivity, & 
                                    GetThermalVelocityOfAir

    IMPLICIT NONE

	Type(Particle),POINTER :: InParticle

	ParticleKnudsenForEnergy = 3 *MolecularThermalDiffusivity() &
			/ GetThermalVelocityOfAir()   &
			/ InParticle%EffectiveRadius

	RETURN
END FUNCTION ParticleKnudsenForEnergy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate the CUNNINGHAM SLIP CORRECTION FACTOR, which !!
!! has a number of possible forms.  This function selects !!
!! the appropriate one and returns the value.			  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION CunninghamSlipCorrectionFactor (InParticle)

    IMPLICIT NONE

	Type(Particle),POINTER :: InParticle
	REAL*8		   :: Knudsen

	Knudsen = ParticleKnudsen (InParticle)

	!! From Kasten 1968 (via Jacobson, 1997)
	CunninghamSlipCorrectionFactor = 1 + Knudsen * (1.249 + 0.42 * dexp(-0.87 / Knudsen))

	!write(*,*) 'CRL: CunninghamSlipCorrectionFactor = ',CunninghamSlipCorrectionFactor
	!write(*,*) 'CRL: Knudsen = ',Knudsen
	
	RETURN
END FUNCTION CunninghamSlipCorrectionFactor



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Returns the pH of an aerosol. !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION AerosolPH (InParticle)

	USE ModelParameters, ONLY : ProtonIndex

	IMPLICIT NONE

	!! External Variables
	TYPE(PARTICLE),POINTER :: InParticle

	!!Internal Variables
	REAL*8 :: ProtonMolarity


	IF (InParticle%NumberOfParticles .EQ. 0. .OR. InParticle%Dry) THEN
		AerosolPH = 14.

		RETURN
	END IF

	IF (InParticle%AqChems(ProtonIndex) .GT. 0) THEN
		!Molarity in mol/L solution
		ProtonMolarity = InParticle%AqChems(ProtonIndex)*InParticle%SolutionDensity*1000/AqueousSolutionMass(InParticle)
		
		AerosolPH = -1.*LOG10(ProtonMolarity)

	ELSE
		AerosolPH = 14.

	END IF


END FUNCTION AerosolPH



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate the PARTICLE DENSITY for a lagrangian particle !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION ParticleDensity (InParticle)
!MJA 05-008-2012 Added special routine for dry BC only particles

	USE ModelParameters,     ONLY : grams, cm, aqchemscale, &
                                        moles, pi, avogadro,    &
                                        SmallestAerosolPossible, & 
									micron ! cmb add

	USE Chemistry,		ONLY : IFSOLID, AqMolecularMass,       &
                                       HowManyAqChems,HowManyAqCations,&
				       HowManyAqAnions,                &
                                       HowManyAqEqReactions,           & 
                                       AqEquilibriaList,	       &
				       AqAnionCharge, AqCationCharge,  &
				       HowManyOrgChems,                &
                                       OrgMolecularMass, OrgDensity,   &
                                       SolidSaltDensity,               &
                                       HowManyAqOrgChems,              &
                                       AqOrgMolecularMass,             &
                                       AqOrgDensity

	USE InfrastructuralCode, ONLY : ERROR,REAL2STR, WARN, IsNaN, int2str

	! cmb add
	!use outputroutines, only : spillbeans
	! 
	
	IMPLICIT NONE

	!! External Variables
	TYPE(Particle),POINTER :: InParticle

	!! Internal Variables
	INTEGER :: I,K, SkipFlag, C
	REAL*8  ::Volume,Mass,ionMass,PSIij, TempElectrolytes(HowManyAqChems),&
	BinaryDensity, NumerSum, DenomSum, TempMolality, &
  	PatKumWeights(HowManyAqEqReactions,3), XX, YY, Residual, dR, &
	SolutionDensity, InsolMass, InsolVolume, InsolDensity, InorgSolnDensity

	REAL*8, PARAMETER :: ResidualThreshold = 0.0025 !**3.

        LOGICAL :: ForceSolnDens

	real*8 :: sum_electrolytes

        !Use this flag to force the aqueous solution density to be 1.0 g/cm3 and skip a lot of code
        ForceSolnDens = .TRUE.
        
        ! CB: Modify this equality check to have a designated error 
        !     tolerance        
        !IF (InParticle%NumberOfParticles .LE. 0. ) THEN
        if (inparticle%numberofparticles .le. 1.0E-40) then
			ParticleDensity = 1.0
			InParticle%EmbryoRadius    = 0.
			InParticle%EffectiveRadius = 0.
			InParticle%InsolubleRadius = 0.
			InParticle%ParticleDensity = 1.
			InParticle%SolutionDensity = 1.
			InParticle%InsolubleDensity = 1.
			InParticle%InorgSolnDensity = 1.
			RETURN
		Else
			! cb add
			!call system("/bin/rm -f fort.2222")
			!call spillbeans(tout_unit=2222)
			!write(*,*) 'Nonzero # of particles (', inparticle%numberofparticles, ')'
			!write(*,*) ''
			!write(*, *)"Number ",trim(INT2STR(I)),",  ",inparticle%numberofparticles,inparticle%sectional
			!write(*, *)"Effective Radius: ",trim(real2str(inparticle%EffectiveRadius/micron)),",  (",trim(real2str(inparticle%Edges(1)/micron))," to ",trim(real2str(inparticle%Edges(2)/micron)),")"
			!write(*, *)"Embryo Radius: ",trim(real2str(inparticle%EmbryoRadius/micron))
			!write(*, *)"Insoluble Radius: ",trim(real2str(inparticle%InsolubleRadius/micron))
			!write(*, *)"Absorbing Core Radius: ",trim(real2str(inparticle%AbsCoreRad/micron))
			!write(*, *)"Shell Real Ref. Ind. 550 nm",trim(real2str(inparticle%ShellRealRefracInd(1)))
			!write(*, *)"Shell Imag. Ref. Ind. 550 nm",trim(real2str(inparticle%ShellImagRefracInd(1)))
			!write(*, *)"Extinction Efficiency: 550 nm",trim(real2str(inparticle%ExtEff(1)))
			!write(*, *)"Scattering Efficiency: 550 nm",trim(real2str(inparticle%ScaEff(1)))
			!write(*, *)"BackScattering Efficiency: 550 nm",trim(real2str(inparticle%BackScaEff(1)))
			!write(*, *)"Assymetry Parameter: 550 nm",trim(real2str(inparticle%AssymParam(1)))
		
			!write(*, *)"Temperature ",trim(real2str(inparticle%Temperature))
			!write(*, *)"Ionic Strength ",trim(real2str(inparticle%IonicStr))
			!write(*, *)"Water Activity ",trim(real2str(inparticle%WaterActivity))
			!write(*, *)"SolutionDensity ",trim(real2str(inparticle%SolutionDensity))
			!!write(*, *)"InsolubleDensity ",trim(real2str(inparticle%InsolubleDensity))
			!write(*, *)"ParticleDensity ",trim(real2str(inparticle%ParticleDensity))
			!write(*, *)"Ph ", trim(real2str(AerosolPh(inparticle)))
			!write(*, *)"Curvature Correction ", trim(real2str(CurvatureCorrection(inparticle)))
			!write(*, *)"OrgCurvature Correction ", trim(real2str(OrgCurvatureCorrection(inparticle)))	
			!write(*, *)"Solution Surface Tension", trim(real2str(inparticle%SurfaceTension))
			
			!write(*, *)"Temperature ", trim(real2str(inparticle%temperature))
			
			if (inparticle%numberofparticles .GT. 0.0) then
				!write(*,*) '--> calculating HypotheticalElectrolyteConcentrations <--'
				!		SaltConc = HypotheticalElectrolyteConcentrations (inparticle)
			endif   
			!write(*,*) ''
			!write(*,*) 'AqChems:'
			!write(*,*) '--------'
			!do i = 1, size(inparticle%aqchems)
		!		write(*,*) 'i = ', i, ', inparticle%aqchems(i) = ', inparticle%aqchems(i)
		!	enddo
		!	write(*,*) ''
			!write(*,*) 'OrgChems:'
			!write(*,*) '---------'
			!do i = 1, size(inparticle%orgchems)
		!		write(*,*) 'i = ', i, ', inparticle%orgchems(i) = ', inparticle%orgchems(i)
		!	enddo
		!	write(*,*) ''
		!	write(*,*) 'AqOrgChems:'
		!	write(*,*) '-----------'
		!	do i = 1, size(inparticle%aqorgchems)
		!		write(*,*) 'i = ', i, ', inparticle%aqorgchems(i) = ', inparticle%aqorgchems(i)
		!	enddo
			! end cb add
		endif

	!MJA 05-08-2012 Special routine for dry (BC only) particles
        IF (InParticle%Dry) THEN
		ParticleDensity = OrgDensity(1) !BC density
		InsolVolume = InParticle%OrgChems(1)*OrgMolecularMass(1)/OrgDensity(1) 
                InParticle%EmbryoRadius    = 0.0
		InParticle%InsolubleRadius = (0.75 * InsolVolume / pi)**0.33333333333
                InParticle%EffectiveRadius = InParticle%InsolubleRadius

		InParticle%ParticleDensity = OrgDensity(1)
		InParticle%SolutionDensity = OrgDensity(1)
		InParticle%InsolubleDensity = OrgDensity(1)
		InParticle%InorgSolnDensity = OrgDensity(1)
		RETURN
	END IF

        Mass = 0 ; Volume = 0; ionMass = 0; NumerSum = 0; DenomSum = 0
        IF (ForceSolnDens) THEN
            InorgSolnDensity = 1.0
            DO I = 1, HowManyAqChems+HowManyAqAnions+HowManyAqCations
                Mass = Mass    + AqMolecularMass(I)*InParticle%AqChems(I)
            ENDDO
            Volume = Mass*InorgSolnDensity
        ELSE        


           !! Get electrolyte equivalents of solution containing dissociated ions,
           !! which assumes a two-parameter fit proscribed in Resch, 1995 
           !! (assuming that the density of water does not change with Temperature).
           !! This is a SOLUTION density, which includes the incorporated water.
           !! So we scale up by the mass of the water as well!
           TempElectrolytes = HypotheticalElectrolyteConcentrations (InParticle)
           ! cmb in: 
           !write(*,*) 'CMB: Output from HypotheticalElectrolyteConcentrations...'
           call flush(6)
           ! end cmb
           do i = 1, size(TempElectrolytes)
              !write(*,*) 'i = ', i, ', TempElectrolytes = ', TempElectrolytes(i)
              sum_electrolytes = sum_electrolytes + TempElectrolytes(i)
           enddo
           ! end cmb in
	
           !For some reason, This routine returns NaN if there isn't a write
           !statement here. MJA, 070507
           WRITE(1234,*) TempElectrolytes(1)
	   ! cb add
           close(1234)

           !! Get the weightings for the Patwardhan and Kumar mixing rule
           !WRITE(*,*)'Just Before PatKumWeights spilling beans'
           !call spillbeans()
           PatKumWeights = GetPatwardhandAndKumarWeightings (InParticle)

           !! Loop over each reaction and consider the appropriate equation
           DO I = 1, HowManyAqEqReactions
              !write(*,*) 'i (of HowManyAqEqReactions) = ', i

              !! Some quasi-electrolytes (such as HSO4-) are not quantifyable in this
              !! density parameterization and should have been pushed to upper level
              !! electrolytes anyways (such as H2SO4).  So we ignore those electrolytes!
              ! CMB (AER, Inc): floating-point inequality check again?
              IF (AqEquilibriaList(I,1) .LE. HowManyAqChems .AND. AqEquilibriaList(I,1) .ge. 1.0e-40) THEN ! GT 0 disallows water
                 !IF (AqEquilibriaList(I,1) .LE. HowManyAqChems .AND. AqEquilibriaList(I,1) .GT. 0.) THEN ! GT 0 disallows water

                 !! This is the equivalent molality of the associated electrolyte
                 !write(*,*) AqPhaseChemicalNames(INT(AqEquilibriaList(I,1))),I
                 !Write(*,*) "Temp Electrolytes: ", TempElectrolytes(AqEquilibriaList(I,1))
                 !Write(*,*) "InParticle%AqChems(1) ", InParticle%AqChems(1)
                 !Write(*,*) "moles ", moles
                 !Write(*,*) "grams ", grams
                 !Write(*,*) "AqMolecularMass(1) ", AqMolecularMass(1)			! cmb: changed int to nint?  not sure if necessary...
                 TempMolality = TempElectrolytes(nINT(AqEquilibriaList(I,1))) / InParticle%AqChems(1) / &
                      (moles / grams * AqMolecularMass(1)) * 1000.

                 !Mass    = Mass    + AqMolecularMass(INT(AqEquilibriaList(I,1)))*TempElectrolytes(INT(AqEquilibriaList(I,1)))
                 !ionMass = ionMass + AqMolecularMass(INT(AqEquilibriaList(I,1)))*TempElectrolytes(INT(AqEquilibriaList(I,1)))
                 Mass    = Mass    + AqMolecularMass(nINT(AqEquilibriaList(I,1)))*TempElectrolytes(nINT(AqEquilibriaList(I,1)))
                 ionMass = ionMass + AqMolecularMass(nINT(AqEquilibriaList(I,1)))*TempElectrolytes(nINT(AqEquilibriaList(I,1)))

                 IF(IsNaN (TempMolality)) THEN
                    WRITE(*,*) "Problem in Temp Molality"
                    WRITE(*,*) InParticle%NumberOfParticles, I, AqEquilibriaList(I,13)
                    WRITE(*,*) TempElectrolytes(INT(AqEquilibriaList(I,1))), TempMolality, InParticle%AqChems(1)
                 END IF

                 !WRITE(*,*) TempMolality, I
                 !! Binary density is an exponential function of the two parameters input into the system
                 IF(TempMolality .LT. 100.0) THEN

                    !WRITE(*,*) "Before Exp, Rxn: ", I
                    !WRITE(*,*) "TempMolality: ", TempMolality
                    BinaryDensity = AqEquilibriaList(I,12) + EXP(-1.*AqEquilibriaList(I,13)*TempMolality)*		&
                         (1.*grams/cm/cm/cm-AqEquilibriaList(I,12))

                    !This prevents the binary density of NH4OH from going below 0.9 g/cm3
                    IF(AqEquilibriaList(I,13) .LT. 0. .AND. BinaryDensity .LT. 0.9) BinaryDensity = 0.9

                 ELSE IF(TempMolality .GE. 100.0 .AND. AqEquilibriaList(I,13) .GT. 0.) THEN
                    !If the molality gets very large, and exponent goes to zero
                    !(This keeps EXP from giving an overflow error)
                    BinaryDensity = AqEquilibriaList(I,12)
                 ELSE !For NH4OH at high molality, keep it above 0.9 g/cm3
                    BinaryDensity = 0.9	
                 END IF

                 !! Employ Mixing Rule of Patwardhan and Kumar (AIChE J. 39, 711-714) to Mix these ions appropriately
                 PSIij =  1000. * PatKumWeights(I,1) + AqMolecularMass(INT(AqEquilibriaList(I,1))) * TempMolality
                 !write(*,*) AqPhaseChemicalNames(INT(AqEquilibriaList(I,1)))
                 !write(*,*) 'PatKumWeights(I,1) = ', PatKumWeights(I,1)
                 !write(*,*) 'AqMolecularMass(INT(AqEquilibriaList(I,1))) = ',AqMolecularMass(INT(AqEquilibriaList(I,1)))	
                 !write(*,*) 'TempMolality = ',TempMolality

                 IF(IsNaN(PatKumWeights(I,1))) THEN
                    CALL WARN("Problem in PatKumWeights! Get NaN!")
                 END IF
                 !! Composite Density is the Numerator Sum divided by the Denomenator Sum, although we will track it as Volume
                 NumerSum = NumerSum + PSIij
                 DenomSum = DenomSum + PSIij / BinaryDensity
              END IF
           END DO

           !WRITE(*,*) "After Ion loop"
           !PAUSE
           !! Now use the composite density to add volume to the system
           Volume = Volume + (ionMass + AqMolecularMass(1)*InParticle%AqChems(1)) * DenomSum / NumerSum
           Mass   = Mass   + AqMolecularMass(1)*InParticle%AqChems(1)
           !write(*,*) 'Volume = ', Volume
           !write(*,*) 'ionMass = ',ionMass	
           !write(*,*) 'AqMolecularMass(1) = ',AqMolecularMass(1)
           !write(*,*) 'InParticle%AqChems(1) = ',InParticle%AqChems(1)
           !write(*,*) 'DenomSum = ',DenomSum
           !write(*,*) 'NumerSum = ',NumerSum

           !WRITE(*,*) "Ion solution Volume: ", Volume, DenomSum, NumerSum

           !If for some reason the above doesn't work, set inorganic solution density to 1.0
           IF(ISNaN(Volume)) THEN
              CALL WARN("Setting ion solution density to 1.0 g/cm3")
              Volume = Mass 
              !Volume = 1.0 
           END IF

           !This only includes water and inorganic ions (needed for optical property calculation)
           InorgSolnDensity = Mass/Volume
           !WRITE(*,*) "Inorganic Solution Density: ", InorgSolnDensity
        ENDIF


	!Now add in contributions from Hydrophilic Organics (assumed constant density)
	!! This assumes a fixed density for these aqueous phase chemicals
	!WRITE(*,*) "Before Aq Org loop Volume = ", Volume
	DO I = 1, HowManyAqOrgChems
		!write(*,*)
		Mass   = Mass   + AqOrgMolecularMass(I)*InParticle%AqOrgChems(I)
		Volume = Volume + InParticle%AqOrgChems(I)*AqOrgMolecularMass(I)/AqOrgDensity(I) 
		!write(*,*) 'AqOrgPhaseChemicalNames(I) = ',AqOrgPhaseChemicalNames(I)
		!write(*,*) 'Volume =', Volume
		!write(*,*) 'Mass = ',Mass	
		!write(*,*) 'AqOrgMolecularMass(I) = ',AqOrgMolecularMass(I)
		!write(*,*) 'InParticle%AqOrgChems(I) = ',InParticle%AqOrgChems(I) 
		!write(*,*) 'AqOrgDensity(I)  = ',AqOrgDensity(I) 
		!write(*,*) 'DenomSum = ',DenomSum
		!write(*,*) 'NumerSum = ',NumerSum

	END DO
	
	!WRITE(*,*) "After Aq Org loop", Volume
! 28 September 2016 CMB: expanding this statement to let the user know why there is a problem
! 7 October 2016 CMB: just set them to 1 again...ugh...
!	IF(ISNaN(Volume)) STOP
!	if (IsNaN(Volume)) then
!		call error("In ParticleAttributes.h; volume is NaN after Aq Org loop, stop statement called")
!	endif
	if (isNaN(volume)) then
		call warn("Setting aq org solution density to 1.0 g/cm3")
		call flush(6)
		volume = mass
	endif

	!! This is purely the solution density (water, ions, non-solid electrolytes (e.g. NH3), aq organics)
	SolutionDensity =  Mass/Volume

	!!Calculate radius of the solution embryo
	InParticle%EmbryoRadius = (0.75 * Volume / pi)**0.33333333333
	!WRITE(*,*) "Embryo Radius: ", InParticle%EmbryoRadius
	
	!! Consider the contributions of Solid salts (except ions and water)
	InsolVolume = 0.
	InsolMass = 0.
	DO I = 2, HowManyAqChems
		SkipFlag = 0
		!Skip non-solid electrolytes (NH3, CO2, etc).
		DO C = 1, HowManyAqEqReactions
			! CMB
			if (dabs(aqequilibrialist(c,1) - I) .lt. 0.000001 .and. &
					dabs(AqEquilibriaList(c,6)) .lt. 1.0e-12) then
			!IF(AqEquilibriaList(C,1) .EQ. I .AND. AqEquilibriaList(C,6) .EQ. 0.0) THEN
				SkipFlag = 1
			END IF
		END DO
			
		IF(SkipFlag .EQ. 0) THEN
			InsolMass   = InsolMass   + AqMolecularMass(I)*InParticle%AqChems(I)
			InsolVolume = InsolVolume + InParticle%AqChems(I)*AqMolecularMass(I)/SolidSaltDensity(I)

			Mass   = Mass   + AqMolecularMass(I)*InParticle%AqChems(I)
			Volume = Volume + InParticle%AqChems(I)*AqMolecularMass(I)/SolidSaltDensity(I)
		END IF
	END DO
	
	!WRITE(*,*) "After Solid Salts loop", Volume
        ! CMB: Revise to add more info to exit
		! update again: set to 1 if NaN
!	IF(ISNaN(Volume)) STOP
        if (isnan(volume)) then
            call warn("volume is NaN after solid salts loop")
		call flush(6)
        	volume = mass	! cmb add
		endif

	 ! cmb add
 !778 continue
 	 ! end cmb add
	 
	!Now add in controbutions from hydrophobic organics	
	DO I = 1, HowManyOrgChems
		!write(*,*) 'In loop'
		!write(*,*) 'InsolMass = ',InsolMass
		!write(*,*) 'OrgMolecularMass(I) = ',OrgMolecularMass(I)	

		IF(ISNaN(InParticle%OrgChems(I))) THEN
			InParticle%OrgChems(I) = 0.
		END IF

		!write(*,*) 'InParticle%OrgChems(I) = ',InParticle%OrgChems(I),OrgPhaseChemicalNames(I)
		!write(*,*) 'InsolVolume = ',InsolVolume
		!write(*,*) 'OrgDensity(I) = ',OrgDensity(I)
	
		InsolMass   = InsolMass   + OrgMolecularMass(I)*InParticle%OrgChems(I)
		InsolVolume = InsolVolume + InParticle%OrgChems(I)*OrgMolecularMass(I)/OrgDensity(I) 
	
		!WRITE(*,*) I, InsolVolume
		Mass   = Mass   + OrgMolecularMass(I)*InParticle%OrgChems(I)

		Volume = Volume + InParticle%OrgChems(I)*OrgMolecularMass(I)/OrgDensity(I) 
	END DO
		
!	WRITE(*,*) "After Hydrophobic Org loop", Volume
        ! CMB: Revise to add more info to exit
		! update: make NaN into 1 and don't exit'
!	IF(ISNaN(Volume)) STOP
        if (isnan(volume)) then
            !call error("In ParticleAttributes.h; volume is NaN after hydrophobic org loop")
		    call warn("volume is NaN after hydrophobic org loop")
        	volume = mass
		endif
	
	!!This is the average density of the non-aqeous part (solid salts, BC, and hydrophobic organics)
	InsolDensity = InsolMass/InsolVolume
	InParticle%InsolubleRadius = (0.75 * InsolVolume / pi)**0.33333333333

	!! This is the average density of the particle
	ParticleDensity = Mass / Volume
	!IF(Mass .LT. 0.0 ) WRITE(*,*) "Mass (<0), Volume = ", Mass, Volume

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! RESET EFFECTIVE RADIUS, too. !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	InParticle%EffectiveRadius = (0.75 * Volume / pi)**0.33333333333


	InParticle%InorgSolnDensity = InorgSolnDensity 
	InParticle%SolutionDensity = SolutionDensity
	InParticle%InsolubleDensity = InsolDensity
	InParticle%ParticleDensity = ParticleDensity

	!WRITE(*,*) "End Particle Density", Volume, Mass
	


	RETURN
END FUNCTION ParticleDensity


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate Particle Mass !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION ParticleMass (InParticle)

	USE ModelParameters,     ONLY : grams, cm, aqchemscale, moles, pi, avogadro
	USE Chemistry,		 ONLY : AqMolecularMass, HowManyAqChems,	&
					HowManyAqAnions, HowManyAqCations,	&
					AqCationMass, AqAnionMass, &
					HowManyOrgChems, OrgMolecularMass, &
					HowManyAqOrgChems, AqOrgMolecularMass 

	IMPLICIT NONE

	!! External Variables
	TYPE(Particle),POINTER :: InParticle

	!! Internal Variables
	INTEGER :: I

	ParticleMass = 0. 

	!! Consider the contributions of Aqueous Chemicals (except ions)
	!! This assumes a fixed density for these aqueous phase chemicals
	DO I = 1, HowManyAqChems
		ParticleMass   = ParticleMass + AqMolecularMass(I)*InParticle%AqChems(I)
	END DO

	DO I = 1, HowManyAqCations
		ParticleMass   = ParticleMass + AqCationMass(I)*InParticle%AqChems(HowManyAqChems+I)
	END DO

	DO I = 1, HowManyAqAnions
		ParticleMass   = ParticleMass + AqAnionMass(I)*InParticle%AqChems(HowManyAqChems+HowManyAqCations+I)
	END DO

	DO I = 1, HowManyOrgChems
		ParticleMass   = ParticleMass + OrgMolecularMass(I)*InParticle%OrgChems(I)
	END DO

	DO I = 1, HowManyAqOrgChems
		ParticleMass   = ParticleMass + AqOrgMolecularMass(I)*InParticle%AqOrgChems(I)
	END DO

	RETURN
END FUNCTION ParticleMass


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Get a vector of weightings for the Patwardhan and Kumar	   !!
!! y_ij factors corrected for partial dissociation, incomplete !!
!! reaction sets, etc.										   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Returns 1). y_ij  2). Cation Index 3). Anion Index		   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Patwardhan and Kumar (AIChE J. 39, 711-714)				   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION GetPatwardhandAndKumarWeightings (InParticle)

	USE Chemistry, ONLY : HowManyAqEqReactions,	&
			      HowManyAqChems,		&
			      HowManyAqCations,		&
			      HowManyAqAnions,		&
			      AqEquilibriaList,		&
			      AqCationCharge,		&
			      AqAnionCharge,            &
			      AqPhaseChemicalNames

	USE InfrastructuralCode, ONLY : WARN, IsNaN

	IMPLICIT NONE

	!! External Variables
	REAL*8 :: GetPatwardhandAndKumarWeightings(HowManyAqEqReactions,3)
	TYPE(Particle),POINTER :: InParticle

	!! Internal Variables
	INTEGER :: K, CatIndex, AnIndex
	REAL*8  :: II, JJ,														  &
			   CationicStrFr(HowManyAqCations),AnionicStrFr(HowManyAqAnions), &
			   CationChgFrac(HowManyAqCations),AnionChgFrac(HowManyAqAnions)

	!! Update the ionic strength
	!WRITE(*,*) "CRL: Before CalculateIonic Strength... "
	!WRITE(*,*) "CRL: Ionic Strength: ", InParticle%IonicStr
	!WRITE(*,*) "CRL: Ionic Strength: ", InParticle%AqChems
	!CALL CalculateIonicStrength(InParticle)
	!WRITE(*,*) "CRL: Ionic Strength: ", InParticle%IonicStr
	!WRITE(*,*) "CRL: Ionic Strength: ", InParticle%AqChems
	IF(IsNAN(InParticle%IonicStr)) CALL WARN("Ionic Strength Error in PKweightings!")

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Make Charge Fraction Vectors !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!! Calculate Cationic First:
	II = 0.
	DO K = 1, HowManyAqCations
		CationChgFrac(K) = CationMolality(K, InParticle) * AqCationCharge(K)
		II               = II + CationChgFrac(K)

	END DO
	DO K = 1, HowManyAqCations
		CationChgFrac(K) = CationChgFrac(K) /II
		IF(IsNAN(CationChgFrac(K))) CALL WARN("CationChgFrac error in PKweightings!")
		IF(IsNAN(CationChgFrac(K))) WRITE(*,*)  CationMolality(K, InParticle),AqCationCharge(K),II
	END DO

	!! And Anionic Next:
	II = 0.
	DO K = 1, HowManyAqAnions
		AnionChgFrac(K) = AnionMolality(K, InParticle) * AqAnionCharge(K)
		II              = II + AnionChgFrac(K)
	END DO
	DO K = 1, HowManyAqAnions
		AnionChgFrac(K) = AnionChgFrac(K) / II
		IF(IsNAN(AnionChgFrac(K))) CALL WARN("AnionChgFrac error in PKweightings!")
	END DO

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Make Ionic Strength Fraction Vectors !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	IF (Inparticle%IonicStr .EQ. 0.) THEN
		DO K = 1, HowManyAqCations
			CationicStrFr(K) = 0.
		END DO
		DO K = 1, HowManyAqAnions
			AnionicStrFr(K)  = 0.
		END DO
	ELSE
		DO K = 1, HowManyAqCations
			CationicStrFr(K) = CationMolality(K, InParticle) * AqCationCharge(K) * AqCationCharge(K) / 2. / InParticle%IonicStr
			IF(IsNAN(CationicStrFr(K))) CALL WARN("CationicStrFr error in PKweightings!")

		END DO
		DO K = 1, HowManyAqAnions
			AnionicStrFr(K)  = AnionMolality(K, InParticle)  * AqAnionCharge(K) * AqAnionCharge(K) /2. / InParticle%IonicStr
			IF(IsNAN(AnionicStrFr(K))) CALL WARN("AnionicStrFr error in PKweightings!")
		END DO
	END IF

	
	!! This will be the sum of all the y_ij Patwardhan and Kumar mixing rule terms
	JJ = 0.

	!! Loop over all reactions
	DO K = 1, HowManyAqEqReactions
		

		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! For partial dissociation, consider the reactions as if the partially !!
		!! dissociated form exchanges directly with the fully dissociated form, !!
		!! which is not how it would have been written.							!!
		IF(AqEquilibriaList(K,1) .LE. HowManyAqChems) THEN ! electrolyte uncharged
			CatIndex = AqEquilibriaList(K,2)
			AnIndex  = AqEquilibriaList(K,3)
		ELSE											   ! electrolyte charged

			!! Is the partially dissociated "electrolyte" a cation or an anion?
			IF (AqEquilibriaList(K,1) .LE. HowManyAqChems+HowManyAqCations) THEN  ! It's a Cation
				CatIndex = AqEquilibriaList(K,1)-HowManyAqChems
				AnIndex  = AqEquilibriaList(K,3)
			ELSE																  ! It's an Anion
				CatIndex = AqEquilibriaList(K,2)
				AnIndex  = AqEquilibriaList(K,1)-HowManyAqChems-HowManyAqCations
			END IF
		END IF
			
		!! Store the indices used to calculate this	
		GetPatwardhandAndKumarWeightings(K,2) = REAL(CatIndex)
		GetPatwardhandAndKumarWeightings(K,3) = REAL(AnIndex)

		
		!Exclude Levitocite reaction from this calculation
		! cmb
		if (dabs(AqEquilibriaList(K, 22)) .gt. 1e-40) then
		!IF(AqEquilibriaList(K,22) .NE. 0) THEN
			IF (AqPhaseChemicalNames(INT(AqEquilibriaList(K,1))) .EQ. "(NH4)3H(SO4)2") THEN
				GetPatwardhandAndKumarWeightings(K,1) = 0.
			END IF
		
		!Normal Calculation
		ELSE
			!! Calculate the Appropriate Weightings for the Patwardhan and Kumar (AIChE J. 39, 711-714)
			!! Mixed Activity Coefficient:
			!!
			!!   ln(Aw_mix) = Sum[(i,j), (Isf[i]*Cf[j]+Isf[j]*Cf[i])*ln(Aw[ij])]
			GetPatwardhandAndKumarWeightings(K,1) = AnionicStrFr(AnIndex)*CationChgFrac(CatIndex) + &
												CationicStrFr(CatIndex)*AnionChgFrac(AnIndex)
		
		END IF
		
		JJ = JJ + GetPatwardhandAndKumarWeightings(K,1)
	END DO

	!! Normalize the y_ij estimates to equal one
	DO K = 1, HowManyAqEqReactions
		GetPatwardhandAndKumarWeightings(K,1) = GetPatwardhandAndKumarWeightings(K,1) / JJ 
	END DO

	RETURN
	!write(*,*) 'GetPatwardhandAndKumarWeightings =',GetPatwardhandAndKumarWeightings 
END FUNCTION GetPatwardhandAndKumarWeightings


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Estimate Dissociated Electrolyte Concentrations !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! For some routines, it is necessary to find the hypothetical
!! electrolyte concentrations that a set of dissociated ions 
!! would make up.  This function provides a vector of such electrolyte
!! abundances, in the form of a vector matching the AqChems vector
!! in the particle structure.
!!
!! The algorithm is taken from Jacobson's Fundamentals of Atmospheric
!! Modeling, 1999 page 494
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION HypotheticalElectrolyteConcentrations (InParticle)

	USE Chemistry, ONLY : HowManyAqChems, HowManyAqCations, &
                              HowManyAqAnions, HowManyAqEqReactions, &
                              AqEquilibriaList, & 
                              AqPhaseChemicalNames, AqAnionNames, &
                              AqCationNames

	USE ModelParameters, ONLY : AqThermoNumbAllowedIterationsForHypLoop, &
                                    moles,WaterEquilibriumIndex

	! CMB: add warn
	USE InfrastructuralCode, ONLY : ERROR, REAL2STR, TRANSCRIPT, warn

	IMPLICIT NONE

	!! External Variables
	REAL*8 :: HypotheticalElectrolyteConcentrations(HowManyAqChems)
	TYPE(Particle),POINTER :: InParticle

	!! Internal Variables
	REAL*8  :: TempChem(HowManyAqChems+HowManyAqCations+HowManyAqAnions), DeltaMoles, Sum, ContentSum, ExcludeRxnInit
	INTEGER :: I, J, K, C, ExcludeRxn 

	!! counters and settings
	K = 0
	ExcludeRxn = 0

	!! First set all electrolytes to zero, unless they are non-solid
	!write(*,*) '(HypotheticalElectrolyteConcentrations) Particle ID: ', inparticle%particleid
10	ContentSum=0
	DO I = 1, HowManyAqChems
		TempChem(I) = 0.
		!write(*,*) 'inparticle%aqchems(', i, ') = ', inparticle%aqchems(i)
		if (inparticle%aqchems(i) .lt. 0) then
			call warn("InParticle%AqChems("//trim(real2str(real(i)))//") < 0 ("//&
					trim(real2str(inparticle%aqchems(i)))//"in HypoElectConc; ID="//trim(real2str(&
					real(inparticle%particleid))))
		endif
		!Check for non-solid electrolyte (NH3, CO2, etc), and add in.
		DO C = 1, HowManyAqEqReactions
			if ((abs(AqEquilibriaList(C, 1) - I) .lt. 0.000001) .and.(&
					abs(AqEquilibriaList(c, 6)) .lt. 0.000001)) then
			!IF(AqEquilibriaList(C,1) .EQ. I .AND. AqEquilibriaList(C,6) .EQ. 0.0) THEN
				!if (inparticle%aqchems(i) .lt. 0) then
				!	call warn("InParticle%AqChems("//trim(real2str(real(i)))//") < 0 ("//&
				!			trim(real2str(inparticle%aqchems(i)))//"in HypoElectConc")
				!endif
				!write(*,*) 'assigning to temp chem for C index using value = ', c, i, inparticle%aqchems(i)
				TempChem(I) = InParticle%AqChems(I)
			END IF
		END DO			
		 
		ContentSum = ContentSum + InParticle%AqChems(I)
	END DO

	!WRITE(*,*) "Check 1 ContentSum: ", ContentSum
	!call flush(6)

	!! Then take the Ion Concentrations as is.
	DO I = HowManyAqChems+1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
		TempChem(I) = InParticle%AqChems(I)
		ContentSum  = ContentSum + InParticle%AqChems(I)
	END DO

	!WRITE(*,*) "Check 2 ContentSum: ", ContentSum

	!WRITE(*,*) "Check 3"

	!! Loop over the reactions and adjust the equations
	DO J = 1, 2
	DO I = 1, HowManyAqEqReactions
		  !WRITE(*,*) I, ExcludeRxn, AqEquilibriaList(I,1), HowManyAqChems
		IF (I .EQ. ExcludeRxn) THEN
			CYCLE

		!! Some reactions may have ions serving as electrolytes
		!! We must identify these, and then push them the other direction
		ELSE IF (AqEquilibriaList(I,1) .GT. HowManyAqChems) THEN

			!WRITE(*,*) TempChem (AqEquilibriaList(I,1))
			!WRITE(*,*) AqCationNames(AqEquilibriaList(I,2)), TempChem (HowManyAqChems+AqEquilibriaList(I,2)), TempChem(AqEquilibriaList(I,1))*AqEquilibriaList(I,4)
			!WRITE(*,*) AqAnionNames(AqEquilibriaList(I,3)), TempChem (HowManyAqChems+HowManyAqCations + AqEquilibriaList(I,3)), TempChem(AqEquilibriaList(I,1))*AqEquilibriaList(I,5)
		
			TempChem (HowManyAqChems+INT(AqEquilibriaList(I,2))) =&
                            TempChem (HowManyAqChems+INT(AqEquilibriaList(I,2))) + &
                            TempChem(INT(AqEquilibriaList(I,1)))*AqEquilibriaList(I,4)
			TempChem (HowManyAqChems+HowManyAqCations + INT(AqEquilibriaList(I,3))) = &
			    TempChem (HowManyAqChems+HowManyAqCations + INT(AqEquilibriaList(I,3))) + &
			    TempChem(INT(AqEquilibriaList(I,1)))*AqEquilibriaList(I,5)
			TempChem (INT(AqEquilibriaList(I,1))) = TempChem(INT(AqEquilibriaList(I,1))) - &
                            TempChem(INT(AqEquilibriaList(I,1)))*AqEquilibriaList(I,5)

			!WRITE(*,*) AqCationNames(AqEquilibriaList(I,2)), TempChem (HowManyAqChems+AqEquilibriaList(I,2))
			!WRITE(*,*) AqAnionNames(AqEquilibriaList(I,3)), TempChem (HowManyAqChems+HowManyAqCations + AqEquilibriaList(I,3))
			!WRITE(*,*) TempChem (AqEquilibriaList(I,1))


		!! Other Reactions are only water dissociations
		ELSE IF (AqEquilibriaList(I,1) .LE. 0) THEN

			!! The adjustment parameter is limited by the concentration of each ion divided
			!! by the appropriate stoicheometric coefficient.
			DeltaMoles  = MIN(TempChem(HowManyAqChems+INT(AqEquilibriaList(I,2)))/AqEquilibriaList(I,4),&
                                TempChem(HowManyAqChems+HowManyAqCations + INT(AqEquilibriaList(I,3)))/&
                                    AqEquilibriaList(I,5))

			!! Now adjust all of the concentrations appropriately
			IF (DeltaMoles .GT. 0) THEN
				!! Don't adjust back towards water at all, just get rid of it
				TempChem (HowManyAqChems+INT(AqEquilibriaList(I,2))) = &
					TempChem (HowManyAqChems+INT(AqEquilibriaList(I,2))) - &
                                        DeltaMoles*AqEquilibriaList(I,4)
				TempChem (HowManyAqChems+HowManyAqCations + INT(AqEquilibriaList(I,3))) = &
					TempChem (HowManyAqChems+HowManyAqCations + &
                                                INT(AqEquilibriaList(I,3))) - DeltaMoles * AqEquilibriaList(I,5)
			END IF

		!! And still others are canonical
		ELSE
			! CMB: more floating-point equality fixes
			!if (abs(AqEquilibriaList(I,22)) .le. 1.0e-40) then
			IF (AqEquilibriaList(I,22) .EQ. 0) THEN
			
				!! The adjustment parameter is limited by the concentration of each ion divided
				!! by the appropriate stoicheometric coefficient.
				DeltaMoles  = MIN(TempChem (HowManyAqChems+INT(AqEquilibriaList(I,2))) /&
                                        AqEquilibriaList(I,4),	& 
                                        TempChem (HowManyAqChems+HowManyAqCations + INT(AqEquilibriaList(I,3))) /&
                                        AqEquilibriaList(I,5))
                                !WRITE(*,*) "Delta Moles: ", DeltaMoles
				!! Now adjust all of the concentrations appropriately
				IF (DeltaMoles .GT. 0) THEN
					TempChem (INT(AqEquilibriaList(I,1))) = TempChem(INT(AqEquilibriaList(I,1))) +&
                                            DeltaMoles
                                        TempChem (HowManyAqChems+INT(AqEquilibriaList(I,2))) = &
                                            TempChem (HowManyAqChems+INT(AqEquilibriaList(I,2))) - &
                                            DeltaMoles * AqEquilibriaList(I,4)
					TempChem (HowManyAqChems+HowManyAqCations + INT(AqEquilibriaList(I,3))) =&
                                            TempChem (HowManyAqChems+HowManyAqCations + &
                                                    INT(AqEquilibriaList(I,3))) - DeltaMoles * AqEquilibriaList(I,5)
				END IF
		

			ELSE
				!Levitocite case: Skip this
				!write(*,*) 'line 830 in particleattributes.h'
				CYCLE						
			END IF
		
		END IF
	END DO
	END DO

	!WRITE(*,*) "Check 4"

	!! Check to make sure have accounted for all of the ions
	SUM = 0.
	DO I = HowManyAqChems+1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
		!IF(IsNan(TempChem(I))) THEN
		!	WRITE(*,*) "Error in HypElec!"
		!	STOP
		!END IF
		SUM = SUM + TempChem(I)
	END DO

	!! This is an error
	!! I changed zero limit to 1.0e-3 from 1.0e-12 to eliminate some problems (Matt Alvarado, 5/10/2006)
	!! Note that this does not error so long as the amount of ions left is less than 1/1000 of the total aqueous moles!
	!! For dry particles, just runs!
	IF (SUM .GT. ContentSum*1.e-3) THEN
	!! It should not be caught in an infinite loop!

		!WRITE(*,*) "Sum, ContentSum: ", SUM, ContentSum*1.0e-3

		!! Show us what was wrong if the calculation doesn't work right.  (Probably the ContentSum*1.e-12 is not right below.)
		CALL TRANSCRIPT("")
		CALL TRANSCRIPT("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		CALL TRANSCRIPT("+++ There are problems in Hypothetical Electrolyte Concentrations +++")
		CALL TRANSCRIPT("")

		DO I = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions

		IF (TempChem(I) .gt. 0) then
		IF (I .LE. HowManyAqChems) CALL TRANSCRIPT(TRIM(AqPhaseChemicalNames(I))//&
                        TRIM(REAL2STR(TempChem(I)/moles)))
		IF (I .GT. HowManyAqChems .AND. I .LE. HowManyAqChems+HowManyAqCations) &
		    CALL TRANSCRIPT(TRIM(AqCationNames(I-HowManyAqChems))//TRIM(REAL2STR(TempChem(I)/moles)))
		IF (I .GT. HowManyAqChems+HowManyAqCations .AND. I .LE. HowManyAqChems+HowManyAqCations+HowManyAqAnions) &
		    CALL TRANSCRIPT(TRIM(AqAnionNames(I-HowManyAqChems-HowManyAqCations))//&
                            TRIM(REAL2STR(TempChem(I)/moles)))
		END IF

		END DO

		CALL DumpParticleContentsAtError(InParticle, "FailedHypotheticalElectrolyteConcentrationsParticle.txt",	&
						 InRelativeHumidity=.TRUE., InRadius=.TRUE., InpH=.TRUE., &
                                                 InDoWaterActivity=.FALSE.,InDensity=.FALSE.,InDoSurfaceTension=.FALSE., &
                                                 InIonicStrength=.TRUE.)
		!Removed this so that routines warns of charge imbalance 
		!in transcript, but does not stop
		! CMB: replace with a warn
		! never mind, it gives a seg fault?
		!CALL warn ("HypotheticalElectrolyteConcentrations() is stuck in what appears to be an infinite loop.")
		!CALL ERROR ("HypotheticalElectrolyteConcentrations() is stuck in what appears to be an infinite loop.")

	END IF

	!! Transfer the good stuff into a shorter vector for output
	DO I = 1, HowManyAqChems
		!WRITE(*,*) I, TempChem(I)
		HypotheticalElectrolyteConcentrations(I) = TempChem(I)
		!WRITE(*,*) 'HypotheticalElectrolyteConcentrations =',TempChem(I), I,AqPhaseChemicalNames
	END DO

	RETURN
END FUNCTION HypotheticalElectrolyteConcentrations

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate the TERMINAL FALL VELOCITY for particles !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION TerminalVelocity (InParticle)

	USE GridPointFields, ONLY : GetAirDensity,GetDynamicViscosityOfAir

	IMPLICIT NONE

	!! Type input variables
	TYPE (Particle),POINTER :: InParticle

	!! Internal Vairables
	REAL*8 :: Dum

	TerminalVelocity = ReynoldsNumber (InParticle) * 0.5             &
					   * GetDynamicViscosityOfAir()  &
					   / GetAirDensity()		 &
					   / InParticle%EffectiveRadius

	RETURN
END FUNCTION TerminalVelocity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Recalculate the Radius of the Particle !!
SUBROUTINE RecalculateRadius (InParticle)

	USE Chemistry, ONLY : AqMolecularMass, HowManyAqChems, IFAQ, OrgPhaseChemicalNames
	!USE outputroutines, only : spillbeans
	IMPLICIT NONE 

	!! Type input variables
	TYPE (Particle),POINTER :: InParticle

	!! Type local variables
	REAL*8  :: II

	! CMB (AER, Inc): Floating-point inequality check
	if (dabs(inparticle%numberofparticles) .le. 1.0e-40) then
	!IF (InParticle%NumberOfParticles .EQ. 0. ) THEN
		IF (InParticle%Edges(1) .GT. 0.) THEN
			InParticle%EmbryoRadius    = InParticle%Edges(1)
			InParticle%InsolubleRadius = 0.
			InParticle%EffectiveRadius = InParticle%Edges(1)
		ELSE
			InParticle%EmbryoRadius    = 0.01*InParticle%Edges(2)
			InParticle%InsolubleRadius = 0.
			InParticle%EffectiveRadius = 0.01*InParticle%Edges(2)
		END IF

		InParticle%ParticleDensity = 1.
		InParticle%SolutionDensity = 1.
		InParticle%InsolubleDensity = 1.
		InParticle%InorgSolnDensity = 1.
		RETURN
	END IF

	!! Calling the Density Routine Resets the Radius, so just do that
	!Write(*,*) "In RecalculateRadius InParticle%AqChems(1) ", InParticle%AqChems(1)
	!Write(*,*) "In RecalculateRadius InParticle%OrgChems(14) ", InParticle%OrgChems(14),OrgPhaseChemicalNames(14)
	!Write(*,*) "Just before calling ParticleDensity spilling beans"
	!call spillbeans()
	II = ParticleDensity (InParticle)
	!Write(*,*) "End ParticleDensity() "
        !! Reset the Surface Tension as well
        II = SurfaceTension (InParticle)	
	!Write(*,*) "End SurfaceTension() "
	!write(*,*) 'Exiting RecalculateRadius()'
	RETURN
END SUBROUTINE RecalculateRadius

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! The correction to pressure from curvature, used for water activity and in condensation. !!
!! It is also known as the "Kelvin Effect"						   !!
!!											   !!
!! The optional argument PerturbWater is a multiple by which to change the water content   !!
!! so we can get a sense of the derivative of the water residual.			   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION CurvatureCorrection (InParticle) 

	USE Chemistry, ONLY : AqMolecularMass, HowManyAqChems, HowManyAqAnions, HowManyAqCations
	USE ModelParameters, ONLY : Rstar, micron, ThermoBulkMode, IgnoreCorrectionstoDiffusivity
	USE GridPointFields, ONLY : SurfaceTensionofWater

	IMPLICIT NONE 

	!! Type input variables
	TYPE (Particle),POINTER  :: InParticle
	REAL*8 :: TotalMoles
	INTEGER :: I

	!!Type Internal variables
	REAL*8 :: TotalDens
	!WRITE(*,*) "Curvature Started", SurfaceTension(InParticle), ParticleDensity(InParticle)
	
	!Ignore this for bulk mode, no particles, or if ordered to.
	IF (ThermoBulkMode .OR. InParticle%NumberofParticles .LE. 0.0D0 .OR. IgnoreCorrectionstoDiffusivity) THEN
	   CurvatureCorrection = 1.
	   !write(*,*) 'since InParticle%NumberofParticles should be le zero, CurvatureCorrection = 1',CurvatureCorrection, InParticle%NumberofParticles 	
	ELSE
	   !! Do this in stages because ParticleDensity forces a radius + solution density update...
	   !WRITE(*,*) "Before Dens"
	   TotalDens  = InParticle%SolutionDensity
	   !WRITE(*,*) "After Dens" 
		
	   !WRITE(*,*) "Curvature Correction: ", InParticle%NumberofParticles
	   !write(*,*) 'InParticle%NumberofParticles should be gt zero', InParticle%NumberofParticles
	   CurvatureCorrection = 2.*AqMolecularMass(1)/(InParticle%Temperature*Rstar*InParticle%SolutionDensity)
	   !WRITE(*,*) CurvatureCorrection, InParticle%EmbryoRadius, SurfaceTension(InParticle)
		
	   CurvatureCorrection = EXP(CurvatureCorrection*InParticle%SurfaceTension/InParticle%EmbryoRadius)
	   !WRITE(*,*) "Curvature Correction: ", CurvatureCorrection
	END IF

	!WRITE(*,*) "Curvature Okay", CurvatureCorrection 
	RETURN
	
END FUNCTION CurvatureCorrection

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Surface Tension depends on Temperature and Electrolytic Properties  !!
!! of the Aerosol's constituent species.			       !!
!!								       !!
!! Follow Li and Lu (2001): Chem. Eng. Sci., 56: 2879-2888, who use    !!
!! a Gibbs dividing surface, define surface excess, etc.  It is a 2-   !!
!! parameter-per-electrolyte system and the data is taken in when the  !!
!! electrolytes are defined.					       !!
!!								       !!
!! Many other systems seem to use the sig = sig_w + Beta I formulation !!
!! from Pruppacher and Klett (1997), but that is significantly less    !!
!! effective than this technique.				       !!
!!                                                                     !!
!! Units are dyn / cm                                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION SurfaceTension (InParticle)

	USE GridPointFields, ONLY : SurfaceTensionOfWater, &
                                    GetRelativeHumidity

	USE Chemistry,       ONLY : AqEquilibriaList,      &
                                    HowManyAqEqReactions
	
        USE ModelParameters, ONLY : UseDonnanSurfaceTension, &
                                    Rstar,                   &
                                    HowManyAqChems,          &
                                    HowManyAqAnions,         &
                                    HowManyAqCations   
	
        USE InfrastructuralCode, ONLY : Error, REAL2STR
	IMPLICIT NONE

	!! External Variables
	TYPE (Particle),POINTER :: InParticle

	!! Internal Variables
	INTEGER :: I
	REAL*8  :: II, Sum1, DenomSum, TotalIonMolality, Test

	! CMB: equality check
	if (inparticle%numberofparticles .le. 1.0e-40) return
	!IF (InParticle%NumberOfParticles .EQ. 0.0) RETURN

	!! CALL FindElectrolyteEquilibrium (InParticle)
	CALL CalculateIonicStrength(InParticle)
        CALL KusikMeissner(InParticle)

	IF(UseDonnanSurfaceTension) THEN !If true, Donnan; If False, Jacobson formulas
		!!Donnan's Parameterization
		!! 1. + SUM (K_i * Activity_i)
		DenomSum = 1.
		DO I = 1, HowManyAqEqReactions
			DenomSum = DenomSum + MeanActivity(I, InParticle) * AqEquilibriaList(I,15)
		END DO

!		IF (Denomsum .GT. 1.0e40) THEN
!			WRITE(*,*) DenomSum, InParticle%numberofParticles
!			CALL DumpParticleContentsAtError(InParticle, "BadActivitiesInSurfaceTensionParticle.txt", &
!                                                 InRelativeHumidity=.TRUE., InRadius=.TRUE., InpH=.TRUE., &
!                                                 InDoWaterActivity=.TRUE.,InDensity=.TRUE.,InDoSurfaceTension=.FALSE., &
!                                                 InIonicStrength=.TRUE.)
!			CALL ERROR ("In SurfaceTension(), the mean activity coefficients appear to have grown without bounds.  "// &
!					"This is likely because the water is leaching out of the system towards zero.  Perhaps the "// &
!					"relative humidity for the particle is too small?  I think it is "//trim(real2str(100*         &
!					GetRelativeHumidity()))//  &
!					"%.  Another possibility is that you are in the water equilibration routine and it is dumping "// &
!					"water from the particle to try to reduce the water activity, and it got too dry.")
!		END IF

		!! Loop over each reaction and consider the appropriate equation
		Sum1 = 0.
		DO I = 1, HowManyAqEqReactions

			!! Ignore dissociation of 
			!IF ((AqEquilibriaList(I,1) .GT. HowManyAqChems) .OR. (AqEquilibriaList(I,1) .LT. 0)) CYCLE

			!! AqIonGrid(I,14) is Surface Excess and AqIonGrid(I,15) is Equilibrium Adsorption Coefficient
			!WRITE(*,*) "Mean Act.: ", MeanActivity(I, InParticle)
			!WRITE(*,*) "Ka: ", AqEquilibriaList(I,15)
			!WRITE(*,*) "DenomSum: ", DenomSum
			Sum1 = Sum1 + AqEquilibriaList(I,14) * LOG (1. - AqEquilibriaList(I,15)	&  
				* MeanActivity(I, InParticle) / DenomSum) 

		END DO

		!! Surface Tesion of the Solution is the Surface Tension of Water + R * T * Sum1
		SurfaceTension = SurfaceTensionOfWater() &
	                 + Rstar*InParticle%Temperature*Sum1
	
		!Now, add in effect of aqueous organics using formula in
		! Jacobson, "Fund. of Atm. Modeling", 2nd ed., p. 534
                !Based on Fachini et al., Nature, 401, 257-259, 1999.
 		SurfaceTension = SurfaceTension - 0.0187*InParticle%Temperature*LOG(1+628.14*MolalityAqCarbon(InParticle))
	
	ELSE !Use formulas from Jacobson, "Fund. of Atm. Modeling", 2nd ed., p. 534		
		SurfaceTension = SurfaceTensionOfWater() 
		!Ions
		TotalIonMolality = 0.
		DO I = HowManyAqChems+1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
			TotalIonMolality = TotalIonMolality + Molality(I, InParticle)
		END DO
                !Based on Pruppacher and Klett, 1997.
		SurfaceTension = SurfaceTension + 1.7*TotalIonMolality
		!Organics
                !Based on Fachini et al., Nature, 401, 257-259, 1999.
		SurfaceTension = SurfaceTension - 0.0187*InParticle%Temperature*LOG(1+628.14*MolalityAqCarbon(InParticle))

	END IF

	!!Force a minimum Surface tension
	!!Use the surface tension of Benzene at 298 K 
	!!(Also used in OrgCurvatureCorrection as the surface tension)		
	!!See Seinfeld and Pandis (1998), Table 10.2, p. 558
	IF(SurfaceTension .LT. 28.21) SurfaceTension = 28.21

	!! Units are dyn / cm
	InParticle%SurfaceTension = SurfaceTension

	RETURN
END FUNCTION SurfaceTension

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate the PARTICLE REYNOLDS NUMBER !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Uses a three-regime (slip, continuum, and transition) !!
!! from Beard's parameterization. (see Jacobson, 1997)   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Tested against Fig. 16.3 in Jacobson, 1999 and it     !!
!! matches well.					 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION ReynoldsNumber (InParticle)

	USE GridPointFields, ONLY : GetAirDensity,GetDynamicViscosityOfAir
	USE ModelParameters, ONLY : g, grams, cm
	IMPLICIT NONE

	!! Type input variables
	TYPE (Particle),POINTER :: InParticle

	!! Type local variables
	REAL*8 :: Dv, Cs, AirDensity, X, Y

	Dv         = GetDynamicViscosityOfAir()
	Cs         = CunninghamSlipCorrectionFactor(InParticle)
	AirDensity = GetAirDensity()

	!! Estimate the Reynolds Number and then Calculate the Real Value, 
	!! this is the estimated terminal velocity
	ReynoldsNumber = 2*InParticle%EffectiveRadius*InParticle%EffectiveRadius*(InParticle%ParticleDensity-AirDensity)*g*Cs/(9*DV)

	!! The Estimated Reynolds Number
	ReynoldsNumber = 2.*InParticle%EffectiveRadius*AirDensity*ReynoldsNumber/DV 

	!! The temporary is the Reynolds Number unless Re > 0.01
	IF (ReynoldsNumber > 0.01 .and. ReynoldsNumber < 300) THEN

		X = dLog((InParticle%ParticleDensity-AirDensity)*AirDensity*InParticle%EffectiveRadius**3.*32.*g/(3.*DV*DV))

		!! The Values here taken from Prupp and Klett, p.417 and originate from Beard's paper (1976)
		Y = -.318657d1+(0.992696d0)*X-(0.153193d-2)*X**2.-(0.987059d-3)*X**3.	&
			-(0.578878d-3)*X**4.+(0.855176d-4)*X**5.-(0.327815d-5)*X**6.

		ReynoldsNumber = Cs*dExp(Y)

	ENDIF

	IF (ReynoldsNumber .GE. 300.) THEN

		!! Use X to signify Np**1/6
		X = (InParticle%SurfaceTension**3.*AirDensity*AirDensity/(DV**4. * (InParticle%ParticleDensity-AirDensity) * g))**(1./6.)

		Y = dLog(16./3.*X*InParticle%EffectiveRadius*InParticle%EffectiveRadius*	&
				 (InParticle%ParticleDensity-AirDensity)*g/InParticle%SurfaceTension)

		ReynoldsNumber = Cs*X*dExp(-5.00015+5.23778*Y-2.04914*Y**2.+0.475294*Y**3.-0.0542819*Y**4.+0.00238449*Y**5.)

	END IF

	RETURN
END FUNCTION ReynoldsNumber

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! The correction to pressure from curvature of the organic phase.			   !!
!! It is also known as the "Kelvin Effect"						   !!
!! This uses approximate values for molecular mass and surface tension                     !!
!! of the organic phase. Since the curvature correction is mainly for                      !!
!! water above RH = 98%, this should be of minor importance.				   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION OrgCurvatureCorrection (InParticle) 

	USE Chemistry, ONLY : AqMolecularMass, HowManyAqChems, HowManyAqAnions, HowManyAqCations
	USE ModelParameters, ONLY : Rstar, micron, ThermoBulkMode, IgnoreCorrectionstoDiffusivity
	USE GridPointFields, ONLY : SurfaceTensionofWater

	IMPLICIT NONE 

	!! Type input variables
	TYPE (Particle),POINTER  :: InParticle
	REAL*8 :: TotalMoles
	INTEGER :: I

	!!Type Internal variables
	REAL*8 :: TotalDens, SurfTens, Radius
	
	!Ignore this for bulk mode or no insoluble species
	IF (ThermoBulkMode .OR. IgnoreCorrectionstoDiffusivity .OR. InParticle%NumberofParticles .EQ. 0.0) THEN
		OrgCurvatureCorrection = 1.

	ELSE
		TotalDens = InParticle%InsolubleDensity	
                !write(*,*) "In OrgCurvatureCorrection ..."
                !write(*,*) "TotalDens = ", TotalDens
		!!Assume a molecular mass of 200 g/mol
		OrgCurvatureCorrection = 2.*200.0/(InParticle%Temperature*Rstar*InParticle%InsolubleDensity)

                !write(*,*) "OrgCurvatureCorrection =",OrgCurvatureCorrection
                !write(*,*) "InParticle%Temperature = ",InParticle%Temperature
		!write(*,*) "Rstar = ",Rstar
                !write(*,*) "InParticle%InsolubleDensity = ",InParticle%InsolubleDensity 

		!Keep radius above 1 nm for calculation
		IF(InParticle%InsolubleRadius .LE. 0.001*micron) THEN
			Radius = 0.001*micron
		ELSE
			Radius = InParticle%InsolubleRadius
		END IF
		
		!!Assume surface Tension of Benzene, Seinfeld and Pandis, p. 522
		SurfTens = 28.21
		OrgCurvatureCorrection = EXP(OrgCurvatureCorrection*SurfTens/Radius)
                !write(*,*) ' SurfTens = ',SurfTens
                !write(*,*) 'Radius = ',Radius

	END IF

	!WRITE(*,*) "OrgCurvature Okay"
	RETURN
	
END FUNCTION OrgCurvatureCorrection

!Molality of carbon (mol/kg H2O) in aqueous phase for surface tension calculation (Matt Alvarado, 6/2006)
REAL*8 FUNCTION MolalityAqCarbon(InParticle)

	USE Chemistry,				ONLY :	AqOrgNumCarbon
	USE ModelParameters
	
	IMPLICIT NONE

	REAL*8  :: Total
	INTEGER :: I
	TYPE(PARTICLE),POINTER :: InParticle

	
	Total = 0.0
	
	DO I = 1, HowManyAqOrgChems
		Total = Total + InParticle%AqOrgChems(I)*AqOrgNumCarbon(I)
	END DO
	
	MolalityAqCarbon = Total / InParticle%AqChems(1) / moles * grams / WaterMolecMass * 1000.

	RETURN

END FUNCTION MolalityAqCarbon

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate Mass of aqueous phase (g/particle) !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION AqueousSolutionMass (InParticle)

	USE ModelParameters,     ONLY : grams, cm, aqchemscale, moles, pi, avogadro
	USE Chemistry,			 ONLY : AqMolecularMass, HowManyAqChems,	&
									HowManyAqAnions, HowManyAqCations,	&
									AqCationMass, AqAnionMass, &
									HowManyAqOrgChems, AqOrgMolecularMass, &
									HowManyAqEqReactions, AqEquilibriaList

	IMPLICIT NONE

	!! External Variables
	TYPE(Particle),POINTER :: InParticle

	!! Internal Variables
	INTEGER :: I, C

	AqueousSolutionMass = 0. 

	!Water
	AqueousSolutionMass   = AqueousSolutionMass + AqMolecularMass(1)*InParticle%AqChems(1)

	!Non-solid Electrolytes
	DO I = 2, HowManyAqChems
		!Check for non-solid electrolyte (NH3, CO2, etc), and add in.
		DO C = 1, HowManyAqEqReactions
			IF(AqEquilibriaList(C,1) .EQ. I .AND. AqEquilibriaList(C,6) .EQ. 0.0) THEN
				AqueousSolutionMass = AqueousSolutionMass + AqMolecularMass(I)*InParticle%AqChems(I)
			END IF
		END DO					 
	END DO

	DO I = 1, HowManyAqCations
		AqueousSolutionMass   = AqueousSolutionMass + AqCationMass(I)*InParticle%AqChems(HowManyAqChems+I)
	END DO

	DO I = 1, HowManyAqAnions
		AqueousSolutionMass   = AqueousSolutionMass + AqAnionMass(I)*InParticle%AqChems(HowManyAqChems+HowManyAqCations+I)
	END DO

	DO I = 1, HowManyAqOrgChems
		AqueousSolutionMass   = AqueousSolutionMass + AqOrgMolecularMass(I)*InParticle%AqOrgChems(I)
	END DO
	!write(*,*) 'AqueousSolutionMass   = ',AqueousSolutionMass 
	RETURN
END FUNCTION AqueousSolutionMass

!! This subroutine calculates the volume-average refractive indices for
!! the reflective shell and then calls DMILAY for the efficienies and 
!! assymetry parameter
SUBROUTINE ShellRefIndAndRad(InParticle)
! 08/30/2010 Changed for the Photolyis Wavebands MJA
! 02/27/2012 Changed to wavelength-dependent refractive indices
! 08/17/2012 Added capability for external mixtures
! 11/08/2012 Added Maxwell Garnett mixing rule
	USE ModelParameters
	USE Chemistry, ONLY : AqMolecularMass, &
			        SolidSaltDensity, AqPhaseChemicalNames, &
				OrgMolecularMass, OrgDensity, &
				AqOrgMolecularMass, AqOrgDensity, &
                                CationIonicRefraction, &
                                AnionIonicRefraction,  &
                                RefIndFlag

	USE InfrastructuralCode, ONLY : INT2STR, REAL2STR, ERROR
	
	IMPLICIT NONE
	
	!! External Variables
	TYPE(Particle),POINTER :: InParticle

	!!External Subroutines 
        EXTERNAL  DMiLay

	!! Internal Variables
	INTEGER, PARAMETER  :: MaxAng = 100
 	INTEGER	:: I, JJ, NumAng, ib, w_ind_low, sizebin, ratiobin
	REAL*8	:: SingleParticleVolume, VolumeConc, Factor,xxx,radratio
	REAL*8	:: Er, Ei, SumEr, SumEi, Wavelength, Sum
	REAL*8	:: SolnMolWeight, MolalRefraction, MolVol, & 
                   TotalVolConc, Rshell, ShellVolumeConc
        REAL*8  :: ErBC, EiBC, ErVA, EiVA, bcvolrat, ErMG, EiMG

	REAL*8, ALLOCATABLE	:: AqVolConc(:), &
				OrgVolConc(:), AqOrgVolConc(:), MolFrac(:), &
                                wvln_array(:)
	INTEGER, ALLOCATABLE	:: wvln_rounddown(:)
        INTEGER :: numbin
        LOGICAL :: FASTTUV

        !These are refractive indices given as data statements
        !!Wavelength-dependent refractive indices from version 3.1a
        !!of Optical Properties of Aerosols and Clouds (OPAC) 
        !! M. Hess, P. Koepke, and I. Schult (1998): 
        !!Optical Properties of Aerosols and clouds: 
        !!The software package OPAC, Bull. Am. Met. Soc., 79, 831-844.
        REAL*8 :: WAVE(61), SOOT_REAL(61), SOOT_IMAG(61)
        REAL*8 :: WATER_REAL(61), WATER_IMAG(61)
        REAL*8 :: WASO_REAL(61), WASO_IMAG(61)
        REAL*8 :: INSO_REAL(61), INSO_IMAG(61)
        REAL*8 :: SEASALT_REAL(61), SEASALT_IMAG(61)
        REAL*8 :: SULF_REAL(61), SULF_IMAG(61)
        REAL*8 :: DUST_REAL(61), DUST_IMAG(61)
        REAL*8 :: a_w, b_w, n_real, n_imag, V_molal, R_molal

	COMPLEX :: SolarShellInd, SolarCoreInd, IRShellInd, IRCoreInd, VAInd
        COMPLEX :: Eshell, EBC, Factor2, EMG, MGInd
	REAL*8	:: Angles(MaxAng), CosAng(MaxAng), M1( MaxAng, 2), &
			M2( MaxAng, 2 ), S21( MaxAng, 2 ), &
		D21( MaxAng, 2 ), SolutionRefracIndex(2), SolutionVolConc

!First specify the wavelength-dependent refractive indices of
!soot, water, organics, sulfate,
!and sea salt.
!From OPAC v3.1a, downloaded from http://www.lrz.de/~uh234an/www/radaer/opac-des.html
!Reference: M. Hess, P. Koepke, and I. Schult (1998): 
!Optical Properties of Aerosols and clouds: 
!The software package OPAC, Bull. Am. Met. Soc., 79, 831-844.
!All values at 0% RH
!Following GEOS-Chem Aerosol Optics memo of Collette Heald, (Jan. 30, 2010)
!Model BC using OPAC SOOT, organics using WASO, sulfate and nitrate salts
!as sulfate, chloride salts as sea salt.

        !Wavelength in nm
        WAVE = 1.0e3*(/ &
        2.500e-01,  3.000e-01,  3.500e-01,  4.000e-01,  4.500e-01, & 
        5.000e-01,  5.500e-01,  6.000e-01,  6.500e-01,  7.000e-01, &
        7.500e-01,  8.000e-01,  9.000e-01,  1.000e+00,  1.250e+00, &
        1.500e+00,  1.750e+00,  2.000e+00,  2.500e+00,  3.000e+00, &
        3.200e+00,  3.390e+00,  3.500e+00,  3.750e+00,  4.000e+00, &
        4.500e+00,  5.000e+00,  5.500e+00,  6.000e+00,  6.200e+00, &
        6.500e+00,  7.200e+00,  7.900e+00,  8.200e+00,  8.500e+00, &
        8.700e+00,  9.000e+00,  9.200e+00,  9.500e+00,  9.800e+00, &
        1.000e+01,  1.060e+01,  1.100e+01,  1.150e+01,  1.250e+01, &
        1.300e+01,  1.400e+01,  1.480e+01,  1.500e+01,  1.640e+01, &
        1.720e+01,  1.800e+01,  1.850e+01,  2.000e+01,  2.130e+01, &
        2.250e+01,  2.500e+01,  2.790e+01,  3.000e+01,  3.500e+01, &
        4.000e+01/)

        SOOT_REAL=(/ & !From file soot00 for black carbon
        1.620e+00,  1.740e+00,  1.750e+00,  1.750e+00,  1.750e+00, & 
        1.750e+00,  1.750e+00,  1.750e+00,  1.750e+00,  1.750e+00, &
        1.750e+00,  1.750e+00,  1.750e+00,  1.760e+00,  1.760e+00, &
        1.770e+00,  1.790e+00,  1.800e+00,  1.820e+00,  1.840e+00, &
        1.860e+00,  1.870e+00,  1.880e+00,  1.900e+00,  1.920e+00, &
        1.940e+00,  1.970e+00,  1.990e+00,  2.020e+00,  2.030e+00, &
        2.040e+00,  2.060e+00,  2.120e+00,  2.130e+00,  2.150e+00, &
        2.160e+00,  2.170e+00,  2.180e+00,  2.190e+00,  2.200e+00, &
        2.210e+00,  2.220e+00,  2.230e+00,  2.240e+00,  2.270e+00, &
        2.280e+00,  2.310e+00,  2.330e+00,  2.330e+00,  2.360e+00, &
        2.380e+00,  2.400e+00,  2.410e+00,  2.450e+00,  2.460e+00, &
        2.480e+00,  2.510e+00,  2.540e+00,  2.570e+00,  2.630e+00, &
        2.690e+00/) 
        SOOT_IMAG=(/ &
       -4.500e-01, -4.700e-01, -4.650e-01, -4.600e-01, -4.550e-01, &
       -4.500e-01, -4.400e-01, -4.350e-01, -4.350e-01, -4.300e-01, &
       -4.300e-01, -4.300e-01, -4.350e-01, -4.400e-01, -4.500e-01, &
       -4.600e-01, -4.800e-01, -4.900e-01, -5.100e-01, -5.400e-01, &
       -5.400e-01, -5.495e-01, -5.600e-01, -5.700e-01, -5.800e-01, &
       -5.900e-01, -6.000e-01, -6.100e-01, -6.200e-01, -6.250e-01, &
       -6.300e-01, -6.500e-01, -6.700e-01, -6.800e-01, -6.900e-01, &
       -6.900e-01, -7.000e-01, -7.000e-01, -7.100e-01, -7.150e-01, &
       -7.200e-01, -7.300e-01, -7.300e-01, -7.400e-01, -7.500e-01, &
       -7.600e-01, -7.750e-01, -7.900e-01, -7.900e-01, -8.100e-01, &
       -8.200e-01, -8.250e-01, -8.300e-01, -8.500e-01, -8.600e-01, &
       -8.700e-01, -8.900e-01, -9.100e-01, -9.300e-01, -9.700e-01, &
       -1.000e+00/) 

        WATER_REAL=(/ & !Refractive index from cuma00 file
        1.362e+00,  1.349e+00,  1.343e+00,  1.339e+00,  1.337e+00, &
        1.335e+00,  1.333e+00,  1.332e+00,  1.331e+00,  1.331e+00, &
        1.330e+00,  1.329e+00,  1.328e+00,  1.327e+00,  1.323e+00, &
        1.321e+00,  1.313e+00,  1.306e+00,  1.261e+00,  1.371e+00, &
        1.478e+00,  1.423e+00,  1.400e+00,  1.369e+00,  1.351e+00, &
        1.332e+00,  1.325e+00,  1.298e+00,  1.265e+00,  1.363e+00, &
        1.339e+00,  1.312e+00,  1.294e+00,  1.286e+00,  1.278e+00, &
        1.272e+00,  1.262e+00,  1.255e+00,  1.243e+00,  1.229e+00, &
        1.218e+00,  1.179e+00,  1.153e+00,  1.126e+00,  1.123e+00, &
        1.146e+00,  1.210e+00,  1.258e+00,  1.270e+00,  1.346e+00, &
        1.386e+00,  1.423e+00,  1.443e+00,  1.480e+00,  1.491e+00, &
        1.506e+00,  1.531e+00,  1.549e+00,  1.551e+00,  1.532e+00, &
        1.519e+00/) 
        WATER_IMAG=(/ &
       -3.350e-08, -1.600e-08, -6.500e-09, -1.860e-09, -1.020e-09, & 
       -1.000e-09, -1.960e-09, -1.090e-08, -1.640e-08, -3.350e-08, &
       -1.560e-07, -1.250e-07, -4.860e-07, -2.890e-06, -8.700e-06, &
       -2.000e-04, -1.000e-04, -1.100e-03, -1.740e-03, -2.720e-01, &
       -9.240e-02, -2.314e-02, -9.400e-03, -3.500e-03, -4.600e-03, &
       -1.340e-02, -1.240e-02, -1.160e-02, -1.070e-01, -8.800e-02, &
       -3.920e-02, -3.210e-02, -3.390e-02, -3.510e-02, -3.670e-02, &
       -3.790e-02, -3.990e-02, -4.150e-02, -4.440e-02, -4.790e-02, &
       -5.080e-02, -6.740e-02, -9.680e-02, -1.420e-01, -2.590e-01, &
       -3.050e-01, -3.700e-01, -3.960e-01, -4.020e-01, -4.270e-01, &
       -4.290e-01, -4.260e-01, -4.210e-01, -3.930e-01, -3.790e-01, &
       -3.700e-01, -3.560e-01, -3.390e-01, -3.280e-01, -3.360e-01, &
       -3.850e-01/) 

        !From file waso00, water-soluble sulfates, nitrates, and organic
        WASO_REAL=(/ & 
        1.530e+00,  1.530e+00,  1.530e+00,  1.530e+00,  1.530e+00, &
        1.530e+00,  1.530e+00,  1.530e+00,  1.530e+00,  1.530e+00, &
        1.530e+00,  1.520e+00,  1.520e+00,  1.520e+00,  1.510e+00, &
        1.510e+00,  1.470e+00,  1.420e+00,  1.420e+00,  1.420e+00, &
        1.430e+00,  1.430e+00,  1.450e+00,  1.452e+00,  1.455e+00, &
        1.460e+00,  1.450e+00,  1.440e+00,  1.410e+00,  1.430e+00, &
        1.460e+00,  1.400e+00,  1.200e+00,  1.010e+00,  1.300e+00, &
        2.400e+00,  2.560e+00,  2.200e+00,  1.950e+00,  1.870e+00, &
        1.820e+00,  1.760e+00,  1.720e+00,  1.670e+00,  1.620e+00, &
        1.620e+00,  1.560e+00,  1.440e+00,  1.420e+00,  1.750e+00, &
        2.080e+00,  1.980e+00,  1.850e+00,  2.120e+00,  2.060e+00, &
        2.000e+00,  1.880e+00,  1.840e+00,  1.820e+00,  1.920e+00, &
        1.860e+00/) 
        WASO_IMAG=(/ &
       -3.000e-02, -8.000e-03, -5.000e-03, -5.000e-03, -5.000e-03, &
       -5.000e-03, -6.000e-03, -6.000e-03, -7.000e-03, -7.000e-03, &
       -8.500e-03, -1.000e-02, -1.300e-02, -1.550e-02, -1.900e-02, &
       -2.250e-02, -1.750e-02, -8.000e-03, -1.200e-02, -2.200e-02, &
       -8.000e-03, -7.050e-03, -5.000e-03, -4.000e-03, -5.000e-03, &
       -1.300e-02, -1.200e-02, -1.800e-02, -2.300e-02, -2.700e-02, &
       -3.300e-02, -7.000e-02, -6.500e-02, -1.000e-01, -2.150e-01, &
       -2.900e-01, -3.700e-01, -4.200e-01, -1.600e-01, -9.500e-02, &
       -9.000e-02, -7.000e-02, -5.000e-02, -4.700e-02, -5.300e-02, &
       -5.500e-02, -7.300e-02, -1.000e-01, -2.000e-01, -1.600e-01, &
       -2.420e-01, -1.800e-01, -1.700e-01, -2.200e-01, -2.300e-01, &
       -2.400e-01, -2.800e-01, -2.900e-01, -3.000e-01, -4.000e-01, &
       -5.000e-01/) 

       !From file ssam00, for sea salt
       SEASALT_REAL = (/ &
       1.510e+00,  1.510e+00,  1.510e+00,  1.500e+00,  1.500e+00, &
       1.500e+00,  1.500e+00,  1.490e+00,  1.490e+00,  1.490e+00, &
       1.490e+00,  1.480e+00,  1.480e+00,  1.470e+00,  1.470e+00, &
       1.460e+00,  1.450e+00,  1.450e+00,  1.430e+00,  1.610e+00, &
       1.490e+00,  1.480e+00,  1.480e+00,  1.470e+00,  1.480e+00, &
       1.490e+00,  1.470e+00,  1.420e+00,  1.410e+00,  1.600e+00, &
       1.460e+00,  1.420e+00,  1.400e+00,  1.420e+00,  1.480e+00, &
       1.600e+00,  1.650e+00,  1.610e+00,  1.580e+00,  1.560e+00, &
       1.540e+00,  1.500e+00,  1.480e+00,  1.480e+00,  1.420e+00, &
       1.410e+00,  1.410e+00,  1.430e+00,  1.450e+00,  1.560e+00, &
       1.740e+00,  1.780e+00,  1.770e+00,  1.760e+00,  1.760e+00, &
       1.760e+00,  1.760e+00,  1.770e+00,  1.770e+00,  1.760e+00, &
       1.740e+00/) 
       SEASALT_IMAG = (/ &
      -5.000e-06, -2.000e-06, -3.240e-07, -3.000e-08, -2.430e-08, &
      -1.550e-08, -1.000e-08, -1.600e-08, -4.240e-08, -2.000e-07, &
      -1.080e-06, -1.950e-06, -4.240e-05, -1.410e-04, -3.580e-04, &
      -5.700e-04, -7.620e-04, -1.000e-03, -4.000e-03, -1.000e-02, &
      -3.000e-03, -2.050e-03, -1.600e-03, -1.400e-03, -1.400e-03, &
      -1.400e-03, -2.500e-03, -3.600e-03, -1.100e-02, -2.200e-02, &
      -5.000e-03, -7.000e-03, -1.300e-02, -2.000e-02, -2.600e-02, &
      -3.000e-02, -2.800e-02, -2.620e-02, -1.800e-02, -1.600e-02, &
      -1.500e-02, -1.400e-02, -1.400e-02, -1.400e-02, -1.600e-02, &
      -1.800e-02, -2.300e-02, -3.000e-02, -3.500e-02, -9.000e-02, &
      -1.200e-01, -1.300e-01, -1.350e-01, -1.520e-01, -1.650e-01, &
      -1.800e-01, -2.050e-01, -2.750e-01, -3.000e-01, -5.000e-01, &
      -1.000e+00/) 

       !From file suso00, for acidic (stratospheric) H2SO4 droplets
       SULF_REAL = (/ &
        1.484e+00,  1.469e+00,  1.452e+00,  1.440e+00,  1.432e+00, &
        1.431e+00,  1.430e+00,  1.429e+00,  1.429e+00,  1.428e+00, &
        1.427e+00,  1.426e+00,  1.425e+00,  1.422e+00,  1.413e+00, &
        1.403e+00,  1.394e+00,  1.384e+00,  1.344e+00,  1.293e+00, &
        1.311e+00,  1.350e+00,  1.376e+00,  1.396e+00,  1.385e+00, &
        1.385e+00,  1.360e+00,  1.337e+00,  1.425e+00,  1.424e+00, &
        1.370e+00,  1.210e+00,  1.140e+00,  1.200e+00,  1.370e+00, &
        1.530e+00,  1.650e+00,  1.600e+00,  1.670e+00,  1.910e+00, &
        1.890e+00,  1.720e+00,  1.670e+00,  1.890e+00,  1.740e+00, &
        1.690e+00,  1.640e+00,  1.610e+00,  1.590e+00,  1.520e+00, &
        1.724e+00,  1.950e+00,  1.927e+00,  1.823e+00,  1.780e+00, &
        1.870e+00,  1.930e+00,  1.920e+00,  1.920e+00,  1.900e+00, &
        1.890e+00/) 
       SULF_IMAG = (/ &
       -1.000e-08, -1.000e-08, -1.000e-08, -1.000e-08, -1.000e-08, &
       -1.000e-08, -1.000e-08, -1.470e-08, -1.670e-08, -2.050e-08, &
       -7.170e-08, -8.630e-08, -2.550e-07, -1.530e-06, -6.940e-06, &
       -1.200e-04, -4.160e-04, -1.260e-03, -3.760e-03, -9.550e-02, &
       -1.350e-01, -1.578e-01, -1.580e-01, -1.310e-01, -1.260e-01, &
       -1.200e-01, -1.210e-01, -1.830e-01, -1.950e-01, -1.650e-01, &
       -1.280e-01, -1.760e-01, -4.880e-01, -6.450e-01, -7.550e-01, &
       -7.720e-01, -6.330e-01, -5.860e-01, -7.500e-01, -6.800e-01, &
       -4.550e-01, -3.400e-01, -4.850e-01, -3.740e-01, -1.980e-01, &
       -1.950e-01, -1.950e-01, -2.050e-01, -2.110e-01, -4.140e-01, &
       -5.900e-01, -4.100e-02, -3.025e-02, -2.352e-02, -2.925e-01, &
       -3.150e-02, -2.000e-01, -1.800e-01, -1.800e-01, -1.900e-01, &
       -2.200e-01/) 

       !From file miam00, for mineral dust (Quartz and clay particles)
       DUST_REAL = (/ &
        1.530e+00,  1.530e+00,  1.530e+00,  1.530e+00,  1.530e+00, &
        1.530e+00,  1.530e+00,  1.530e+00,  1.530e+00,  1.530e+00, &
        1.530e+00,  1.530e+00,  1.530e+00,  1.530e+00,  1.530e+00, &
        1.530e+00,  1.530e+00,  1.530e+00,  1.520e+00,  1.520e+00, &
        1.510e+00,  1.510e+00,  1.510e+00,  1.500e+00,  1.500e+00, &
        1.500e+00,  1.480e+00,  1.460e+00,  1.440e+00,  1.430e+00, &
        1.420e+00,  1.460e+00,  1.220e+00,  1.120e+00,  1.060e+00, &
        1.190e+00,  1.850e+00,  2.220e+00,  2.940e+00,  2.910e+00, &
        2.570e+00,  1.910e+00,  1.830e+00,  1.810e+00,  1.740e+00, &
        2.000e+00,  1.630e+00,  1.540e+00,  1.510e+00,  1.470e+00, &
        1.490e+00,  1.770e+00,  2.050e+00,  2.200e+00,  2.390e+00, &
        2.690e+00,  2.990e+00,  2.570e+00,  2.420e+00,  2.420e+00, &
        2.340e+00/) 
       DUST_IMAG = (/ &
       -3.000e-02, -2.500e-02, -1.700e-02, -1.300e-02, -8.500e-03, &
       -7.800e-03, -5.500e-03, -4.500e-03, -4.500e-03, -4.000e-03, &
       -4.000e-03, -4.000e-03, -4.000e-03, -4.000e-03, -5.000e-03, &
       -5.700e-03, -6.400e-03, -7.600e-03, -1.400e-02, -3.900e-02, &
       -2.400e-02, -1.925e-02, -1.800e-02, -1.200e-02, -6.700e-03, &
       -8.700e-03, -1.800e-02, -3.600e-02, -5.500e-02, -6.300e-02, &
       -5.200e-02, -1.300e-01, -8.900e-02, -1.200e-01, -2.100e-01, &
       -2.900e-01, -4.400e-01, -5.400e-01, -6.500e-01, -6.500e-01, &
       -5.000e-01, -2.500e-01, -2.000e-01, -3.500e-01, -5.000e-01, &
       -3.500e-01, -2.200e-01, -2.400e-01, -2.600e-01, -3.200e-01, &
       -3.700e-01, -4.700e-01, -5.700e-01, -8.200e-01, -9.400e-01, &
       -1.000e+00, -8.000e-01, -7.800e-01, -6.700e-01, -6.200e-01, &
       -7.000e-01/) 


	!! BUILD N ARRAY for AqSpecies Volume Concentrations
	ALLOCATE (AqVolConc(HowManyAqChems+HowManyAqCations+HowManyAqAnions), STAT = I)
	IF (I > 0) CALL ERROR("In ShellRefIndAndRad(), couldn't allocate "// &
                               "AqVolConc at size "//	&
	     TRIM(INT2STR(HowManyAqChems+HowManyAqCations+HowManyAqAnions)))

	!! BUILD N ARRAY for Org Species Volume Concentrations
	ALLOCATE (OrgVolConc(HowManyOrgChems), STAT = I)
	IF (I > 0) CALL ERROR("In ShellRefIndAndRad(), couldn't allocate OrgVolConc at size "	&
						  //TRIM(INT2STR(HowManyOrgChems)))

	!! BUILD N ARRAY for Aq Org Species Volume Concentrations
	ALLOCATE (AqOrgVolConc(HowManyAqOrgChems), STAT = I)
	IF (I > 0) CALL ERROR("In ShellRefIndAndRad(), couldn't allocate AqOrgVolConc at size "	&
						  //TRIM(INT2STR(HowManyAqOrgChems)))

	!! BUILD N ARRAY for H2O and Ion Mol Fractions
	ALLOCATE (MolFrac(1+HowManyAqCations+HowManyAqAnions), STAT = I)
	IF (I > 0) CALL ERROR("In ShellRefIndAndRad(), couldn't allocate MolFrac at size "	&
						  //TRIM(INT2STR(HowManyAqOrgChems)))
	
	IF (InParticle%numberofparticles .GT. 1.0e-6) THEN
		SingleParticleVolume = (4.0*Pi/3.0)*(InParticle%effectiveradius)**3 !cm3
		
		VolumeConc = SingleParticleVolume*InParticle%numberofparticles !cm3/cm3
		
		!Convert aqueous concentrations to cm3 solution/cm3 air
		
		!Water
		AqVolConc(1)  = InParticle%AqChems(1) * &
                                InParticle%NumberOfParticles*AqMolecularMass(1)/ &
                                InParticle%SolutionDensity
		SolutionVolConc = AqVolConc(1)
		TotalVolConc = AqVolConc(1)
				
		!Solid Salts
		DO I = 2, HowManyAqChems
			AqVolConc(I)  = InParticle%AqChems(I) * &
                                        InParticle%NumberOfParticles*AqMolecularMass(I)/ &
                                        SolidSaltDensity(I)
			TotalVolConc = TotalVolConc + AqVolConc(I)
		END DO

		!Inorganic Ions
		DO I = HowManyAqChems+1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
			AqVolConc(I)  = InParticle%AqChems(I) * &
                                        InParticle%NumberOfParticles*AqMolecularMass(I)/ &
                                        InParticle%SolutionDensity
			TotalVolConc = TotalVolConc + AqVolConc(I)
			SolutionVolConc = SolutionVolConc + AqVolConc(I)	
		END DO

		!Convert org concentrations to cm3 solution/cm3 air
		DO I = 1, HowManyOrgChems
			OrgVolConc(I)  = InParticle%OrgChems(I) * &
                                         InParticle%NumberOfParticles*OrgMolecularMass(I)/ &
                                         OrgDensity(I)
			TotalVolConc = TotalVolConc + OrgVolConc(I)	
		END DO
	
		!Convert hydrophilic org concentrations to cm3 solution/cm3 air
		DO I = 1, HowManyAqOrgChems
			AqOrgVolConc(I)  = InParticle%AqOrgChems(I) * &
                                           InParticle%NumberOfParticles*AqOrgMolecularMass(I) &
                                           /AqOrgDensity(I)
			TotalVolConc = TotalVolConc + AqOrgVolConc(I)	
		END DO

		!Calculate shell volume conc (Total - BC)
		ShellVolumeConc = VolumeConc - OrgVolConc(1)
			
	ELSE    !No particles, use dummy numbers
		InParticle%ExtEff = 1.0
		InParticle%ScaEff = 1.0
		InParticle%BackScaEff = 1.0
		InParticle%AssymParam = 1.0
                RETURN          
	END IF

        
        !Calculate shell and core refractive index for each wavelength bin
        FASTTUV = .TRUE.
        IF(FASTTUV) THEN !Use center wavelengths for 7 tropospheric fastTUV bins
           numbin = 7
	   ALLOCATE (wvln_array(numbin), STAT = I)
	   IF (I > 0) CALL ERROR("In ShellRefIndAndRad(), couldn't allocate "// &
                               "wvln_array at size "//TRIM(INT2STR(numbin)))
	   ALLOCATE (wvln_rounddown(numbin), STAT = I)
	   IF (I > 0) CALL ERROR("In ShellRefIndAndRad(), couldn't allocate "// &
                               "wvln_rounddown at size "//TRIM(INT2STR(numbin)))
           wvln_array(1) = 297.7  !in nm 
           wvln_rounddown(1) = 1  
           wvln_array(2) = 309.5       
           wvln_rounddown(2) = 2 
           wvln_array(3) = 325.5       
           wvln_rounddown(3) = 2 
           wvln_array(4) = 378.8      
           wvln_rounddown(4) = 3 
           wvln_array(5) = 447.5    
           wvln_rounddown(5) = 4 
           wvln_array(6) = 602.0
           wvln_rounddown(6) = 8 
           wvln_array(7) = 736.3 
           wvln_rounddown(7) = 10 
        ELSE !For photolysis between 250 and 700 nm, every 1 nm
           numbin = 451
	   ALLOCATE (wvln_array(numbin), STAT = I)
	   IF (I > 0) CALL ERROR("In ShellRefIndAndRad(), couldn't allocate "// &
                               "wvln_array at size "//TRIM(INT2STR(numbin)))
           DO ib = 1, numbin 
              wvln_array(ib) = (250+(ib-1)*1) !in nm
           ENDDO
        ENDIF
 
        DO ib = 1, numbin

           Wavelength = wvln_array(ib) !in nm

           !Index of point on array WAVE closest to, but less than,
           !the wavelength for interpolation
           IF(FASTTUV) THEN
              w_ind_low = wvln_rounddown(ib)
           ELSE
              w_ind_low = INT(ib/50)+1
           ENDIF

           !Calculate weightings of wavelength-dependent refractive index points
           a_w = (WAVE(w_ind_low+1)-Wavelength) &
                  /(WAVE(w_ind_low+1)-WAVE(w_ind_low))
           b_w = 1-a_w

           !Calculate Core Refractive Index
           n_real = a_w*SOOT_REAL(w_ind_low)+b_w*SOOT_REAL(w_ind_low+1)
           n_imag = -1.0*(a_w*SOOT_IMAG(w_ind_low)+b_w*SOOT_IMAG(w_ind_low+1))
	   SolarCoreInd = CMPLX(n_real, -1.0*n_imag)
           !IF (ib .eq. 419 .or. ib .eq. 420) WRITE(*,*) Wavelength, SolarShellInd, SolarCoreInd

           !We assume the core is BC, the shell everything else
           !BUT we need to handle the cases where the particle is all BC (no shell)
           ! or has no BC (all core). All core handled here.
           IF (ShellVolumeConc .LE. 0 .OR. InParticle%Dry) THEN
                !Shell and core refractive indices equal.
                SolarShellInd = SolarCoreInd
                !Calculate Core radius
	        InParticle%AbsCoreRad = OrgVolConc(1)/InParticle%NumberofParticles
	        InParticle%AbsCoreRad = (0.75*InParticle%AbsCoreRad/Pi)**(1.0/3.0) !In cm
                !Set Shell radius to above, and core a little less than shell
                InParticle%EffectiveRadius = InParticle%AbsCoreRad
                InParticle%AbsCoreRad = 0.99*InParticle%AbsCoreRad
           ELSE !Some or all shell, all shell handled at bottom
	     !Calculate shell refractive index
	     !Use the volume average dielectric constant mixing rule
	     !See Jacobson, "Fund. of Atm. Modeling", 2nd ed, p. 310
 	     SumEr = 0.
	     SumEi = 0.
 
             !Water refractive index and molal Refraction
             n_real = a_w*WATER_REAL(w_ind_low)+b_w*WATER_REAL(w_ind_low+1)
             n_imag = -1.0*(a_w*WATER_IMAG(w_ind_low)+b_w*WATER_IMAG(w_ind_low+1))
             V_molal = 18.0674 !Tang, 1997: density = 0.9971 g/cm3, MW = 18.015 g/mol, V = MW/dens
	     R_molal = V_molal*(n_real**2 - 1.0)/(n_real**2 + 2) 

             !Calculate average refractive index for ionic solution
	     !Following Tang et. al., JGR V.102, D19, pp. 23269-23275, 1997. 
             !and Tang, JGR V.102, D2, pp. 1883-1893, 1997.)
	     !Step 1: Calculate mole fractions in solution
	     Sum = InParticle%AqChems(1)
	     DO I = HowManyAqChems+1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
	   	Sum = Sum + InParticle%AqChems(I)
	     END DO
	     MolFrac(1) = InParticle%AqChems(1)/Sum
	     DO I = 2, 1+HowManyAqCations+HowManyAqAnions
	   	MolFrac(I) = InParticle%AqChems(HowManyAqChems+I-1)/Sum
	     END DO
		
	     !Step 2: Calculate average molecular weight	
	     SolnMolWeight = MolFrac(1)*AqMolecularMass(1)
	     DO I = 2, 1+HowManyAqCations+HowManyAqAnions
	   	SolnMolWeight = SolnMolWeight + MolFrac(I)*AqMolecularMass(HowManyAqChems+I-1)
	     END DO
		
	     !Step 3: Calculate average molal refraction
	     MolalRefraction = MolFrac(1)*R_molal !See Tang et. al, 1997
	     DO I = 2, 1+HowManyAqCations
                MolalRefraction = MolalRefraction + MolFrac(I)*CationIonicRefraction(I-1)
                !WRITE(*,*) MolFrac(I), CationIonicRefraction(I-1)
	     END DO
	     DO I = 2+HowManyAqCations, 1+HowManyAqCations+HowManyAqAnions
	   	MolalRefraction = MolalRefraction + MolFrac(I)*AnionIonicRefraction(I-HowManyAqCations-1)
                !WRITE(*,*) MolFrac(I), AnionIonicRefraction(I-HowManyAqCations-1)
	     END DO

	     !Step 4: Calculate solution molal volume
	     MolVol = SolnMolWeight/InParticle%InorgSolnDensity

	     !Step 5: Calculate Refractive Index
	     xxx = MolalRefraction/MolVol
	     SolutionRefracIndex(1) = SQRT((2.*xxx+1.)/(1.-xxx))
	     SolutionRefracIndex(2) = n_imag !Imaginary index for water
	  
	     !Step 6: Calculate Er and Ei
	     Er = SolutionRefracIndex(1)**2 - SolutionRefracIndex(2)**2
	     Ei = 2*SolutionRefracIndex(1)*SolutionRefracIndex(2)
             SumEr = SumEr + (SolutionVolConc/ShellVolumeConc)*Er
	     SumEi = SumEi + (SolutionVolConc/ShellVolumeConc)*Ei

             !Solid Salts
             DO I = 2, HowManyAqChems
                IF(RefIndFlag(I) .EQ. 0) THEN
                   n_real = a_w*WATER_REAL(w_ind_low)+b_w*SULF_REAL(w_ind_low+1)
                   n_imag = -1.0*(a_w*WATER_IMAG(w_ind_low)+b_w*SULF_IMAG(w_ind_low+1))
                ELSE IF(RefIndFlag(I) .EQ. 1) THEN
                   n_real = a_w*SULF_REAL(w_ind_low)+b_w*SULF_REAL(w_ind_low+1)
                   n_imag = -1.0*(a_w*SULF_IMAG(w_ind_low)+b_w*SULF_IMAG(w_ind_low+1))
                ELSE IF(RefIndFlag(I) .EQ. 2) THEN
                   n_real = a_w*SEASALT_REAL(w_ind_low)+b_w*SEASALT_REAL(w_ind_low+1)
                   n_imag = -1.0*(a_w*SEASALT_IMAG(w_ind_low)+b_w*SEASALT_IMAG(w_ind_low+1))
                ELSE IF(RefIndFlag(I) .EQ. 3) THEN
                   n_real = a_w*DUST_REAL(w_ind_low)+b_w*DUST_REAL(w_ind_low+1)
                   n_imag = -1.0*(a_w*DUST_IMAG(w_ind_low)+b_w*DUST_IMAG(w_ind_low+1))
                ENDIF
	        Er = n_real**2 - n_imag**2
	        Ei = 2*n_real*n_imag
                SumEr = SumEr + (AqVolConc(I)/ShellVolumeConc)*Er
		SumEi = SumEi + (AqVolConc(I)/ShellVolumeConc)*Ei
	     END DO
         
             !Insoluble organics, skipping BC core
	     n_real = a_w*WASO_REAL(w_ind_low)+b_w*WASO_REAL(w_ind_low+1)
             n_imag = -1.0*(a_w*WASO_IMAG(w_ind_low)+b_w*WASO_IMAG(w_ind_low+1))
	     Er = n_real**2 - n_imag**2
	     Ei = 2*n_real*n_imag
             DO I = 2, HowManyOrgChems
		SumEr = SumEr + (OrgVolConc(I)/ShellVolumeConc)*Er
		SumEi = SumEi + (OrgVolConc(I)/ShellVolumeConc)*Ei
	     END DO

             !Hydrophilic Organics
	     n_real = a_w*WASO_REAL(w_ind_low)+b_w*WASO_REAL(w_ind_low+1)
             n_imag = -1.0*(a_w*WASO_IMAG(w_ind_low)+b_w*WASO_IMAG(w_ind_low+1))
	     Er = n_real**2 - n_imag**2
	     Ei = 2*n_real*n_imag
             DO I = 1, HowManyAqOrgChems
		SumEr = SumEr + (AqOrgVolConc(I)/ShellVolumeConc)*Er
		SumEi = SumEi + (AqOrgVolConc(I)/ShellVolumeConc)*Ei
	     END DO

            !!For core-in-shell mixing rule
	    Factor = SQRT(SumEr**2 + SumEi**2)
	    InParticle%ShellRealRefracInd(ib) = SQRT((Factor + SumEr)/2.0)
	    InParticle%ShellImagRefracInd(ib) = SQRT((Factor - SumEr)/2.0)
	    SolarShellInd = CMPLX(InParticle%ShellRealRefracInd(ib),-1.0*InParticle%ShellImagRefracInd(ib))

            !We assume the core is BC, the shell everything else
            !BUT we need to handle the cases where the particle is all BC (all core, no shell)
            ! or has no BC (all shell, no core). No core, all shell handeled here.

            IF (OrgVolConc(1) .LE. 0.0) THEN !No core
               !Set core and shell refractive indices equal
               SolarCoreInd = SolarShellInd
               !Set small core radius
               InParticle%AbsCoreRad = 0.01*InParticle%EffectiveRadius
            ELSE !Both Shell and Core
	       !Calculate Core radius
	       InParticle%AbsCoreRad = OrgVolConc(1)/InParticle%NumberofParticles
	       InParticle%AbsCoreRad = (0.75*InParticle%AbsCoreRad/Pi)**(1.0/3.0) !In cm
	       !This shouldn't happen, but be safe
               IF(InParticle%AbsCoreRad .GT. 0.9999*InParticle%EffectiveRadius) & 
                InParticle%AbsCoreRad = 0.9999*InParticle%EffectiveRadius
            ENDIF

            !!For volume average dielectric constant mixing rule
!	     n_real = a_w*SOOT_REAL(w_ind_low)+b_w*SOOT_REAL(w_ind_low+1)
!            n_imag = -1.0*(a_w*SOOT_IMAG(w_ind_low)+b_w*SOOT_IMAG(w_ind_low+1))
!            ErBC = n_real**2 - n_imag**2
!            EiBC = 2*n_real*n_imag
!            BCVolRat = (InParticle%AbsCoreRad/InParticle%EffectiveRadius)**3.0
!            ErVA = ErBC*BCVolRat+SumEr*(1.0-BCVolRat)
!            EiVA = EiBC*BCVolRat+SumEi*(1.0-BCVolRat)
!            Factor = SQRT(ErVA**2 + EiVA**2)
!	    InParticle%VARealRefracInd(ib) = SQRT((Factor + ErVA)/2.0)
!	    InParticle%VAImagRefracInd(ib) = SQRT((Factor - ErVA)/2.0)            
!            VAInd = CMPLX(InParticle%VARealRefracInd(ib),-1.0*InParticle%VAImagRefracInd(ib))

            !!For Maxwell-Garnett mixing rule
!            Eshell =  CMPLX(SumEr,-1.0*SumEi)
!            EBC =  CMPLX(ErBC,-1.0*EiBC)
!            Factor2 = (EBC-Eshell)/(EBC+2.0*Eshell)
!            EMG = Eshell*(1.0+(3.0*BCVolRat*Factor2/(1.0-BCVolRat*Factor2)))
!            ErMG = REAL(EMG)
!            EiMG = AIMAG(EMG)
!            Factor = SQRT(ErMG**2 + EiMG**2)
!	     InParticle%MGRealRefracInd(ib) = SQRT((Factor + ErMG)/2.0)
!	     InParticle%MGImagRefracInd(ib) = SQRT((Factor - ErMG)/2.0)            
!            MGInd = CMPLX(InParticle%MGRealRefracInd(ib),-1.0*InParticle%MGImagRefracInd(ib))
          ENDIF

	  NumAng = 0.

          !DMiLay needs 2*Pi divided by the wavelength
          Wavelength = 2*PI/((Wavelength)*nm)

!    DMiLay Documentation				
!    I N P U T   A R G U M E N T S

!    (Definition:  size parameter = sphere circumference / wavelength )

!      Rshell      radius of shell

!      Rcore       radius of core

!      WVNO        2*pi / wavelength

!      RindSh      COMPLEX refractive index of shell (negative
!                     imaginary part)

!      RindCo      COMPLEX refractive index of core (negative
!                     imaginary part)

!      MU          array of cosines of scattering angles (angles between
!                     directions of incident and scattered radiation).
!                     For angles between 90 and 180 degrees, use the
!                     supplement (180-angle) of the angle instead, so
!                     0.le.MU.le.1 (see comments below on M1,M2,21,D21)

!      NumAng      Number of scattering angles for which computations
!                     are required; should not exceed MaxAng
!                     (NOTE:  NumAng=0 will suppress the calculation
!                      of the scattering matrix quantitities  M1, M2,
!                      S21, D21 and save a lot of computer time)

!      MaxAng      First dimension of M1,M2,21,D21 in calling program



!    O U T P U T   A R G U M E N T S

!      (Definitions for these arguments can be found in H.C. van de
!       Hulst, Light Scattering By Small Particles, Dover Press, New
!       York, 1981 (reprint of 1957 edition); abbreviated VDH below)

!      QEXT     Efficiency factor for extinction (VDH Sec 9.32)
!               (same as corresponding quantity in MIEV)

!      Qsca     Efficiency factor for scattering (VDH Sec 9.32)
!               (same as corresponding quantity in MIEV)

!      GQS!     average(cosine theta) * Qsca (VDH Sec 9.32)
!                  (<cos theta> is usually denoted by g, hence
!                   the name of the variable)
!               (same as corresponding quantity in MIEV)

!      QBS      Backscatter cross section.
!               ( Re(SBACK)**2 + Im(SBACK)**2 ) / (Pi*XSHELL**2)
!               where the corresponding quantity from MIEV is

!               SBACK = 0.5*sum(n=1 to inf)((-1)**(n+1)(2n+1)(an-bn))

!               and an,bn are ACOE,BCOE below.

!      M1(j,k)  Element M1 of scattering matrix F' (VDH Sec 5.14);
!                  M1(j,1) refers to angle with cosine MU(j); 
!                  M1(j,2) refers to supplement of that angle.
!               (Be sure to type REAL in calling program.)

!      M2(j,k)  Element M2 of scattering matrix F' (VDH Sec 5.14);
!                  M2(j,1) refers to angle with cosine MU(j); 
!                  M2(j,2) refers to supplement of that angle.
!               (Be sure to type REAL in calling program.)

!     S21(j,k)  Element S21 of scattering matrix F' (VDH Sec 5.14);
!                  S21(j,1) refers to angle with cosine MU(j); 
!                  S21(j,2) refers to supplement of that angle.

!     D21(j,k)  Element D21 of scattering matrix F' (VDH Sec 5.14);
!                  D21(j,1) refers to angle with cosine MU(j); 
!                  D21(j,2) refers to supplement of that angle.

          !MM DEBUG - CosAng, M1, M2, S21, and D21 are all zero arrays
          !WRITE(*,*) "InParticle%AbsCoreRad: ", InParticle%AbsCoreRad
          !WRITE(*,*) "InParticle%EffectiveRadius: ", InParticle%EffectiveRadius
          !WRITE(*,*) "Wavelength: ", Wavelength
          !WRITE(*,*) "SolarShellInd: ", SolarShellInd
          !WRITE(*,*) "SolarCoreInd: ", SolarCoreInd
          !WRITE(*,*) "CosAng: ", CosAng
          !WRITE(*,*) "NumAng: ", NumAng 
          !WRITE(*,*) "InParticle%ExtEff(ib): ", InParticle%ExtEff(ib)
          !WRITE(*,*) "InParticle%ScaEff(ib): ", InParticle%ScaEff(ib)
          !WRITE(*,*) "InParticle%BackScaEff(ib): ", InParticle%BackScaEff(ib)
          !WRITE(*,*) "InParticle%AssymParam(ib): ", InParticle%AssymParam(ib)
          !WRITE(*,*) "M1: ", M1
          !WRITE(*,*) "M2: ", M2			
          !WRITE(*,*) "S21: ", S21
          !WRITE(*,*) "D21: ", D21  
          !WRITE(*,*) "MAXANG: ", MAXANG
!          WRITE(*,*) "Call 1"


!MM-DEBUG	  CALL DMiLay( InParticle%AbsCoreRad, &
!MM-DEBUG			InParticle%EffectiveRadius, &
!MM-DEBUG			Wavelength, &
!MM-DEBUG			SolarShellInd, & 
!MM-DEBUG			SolarCoreInd, &
!MM-DEBUG			CosAng, &
!MM-DEBUG			NumAng, &
!MM-DEBUG			InParticle%ExtEff(ib), &
!MM-DEBUG			InParticle%ScaEff(ib), &
!MM-DEBUG			InParticle%BackScaEff(ib), &
!MM-DEBUG			InParticle%AssymParam(ib), & 
!MM-DEBUG                        M1, M2, S21, D21, MAXANG )

          !MM-DEBUG: Replace DMiLay call with if statements to assign values
          !  to InParticle%ExtEff(ib),InParticle%ScaEff(ib),
          !  InParticle%BackScaEff(ib),InParticle%AssymParam(ib) 
          !  for each wavelength bin  (ib), effective radius size bin (sizebin),
          !  and core radius ratio (ratiobin)

          !identify size bin for effective radius (cm)
          IF ((InParticle%EffectiveRadius .GT. 0.) .AND. &
              (InParticle%EffectiveRadius .LE. 2.5e-06)) THEN
            sizebin = 1
          ELSE IF ((InParticle%EffectiveRadius .GT. 2.5e-06) .AND. &
                   (InParticle%EffectiveRadius .LE. 4.0e-06)) THEN
            sizebin = 2
          ELSE IF ((InParticle%EffectiveRadius .GT. 4.0e-06) .AND. &
                   (InParticle%EffectiveRadius .LE. 6.3e-06)) THEN
            sizebin = 3
          ELSE IF ((InParticle%EffectiveRadius .GT. 6.3e-06) .AND. &
                   (InParticle%EffectiveRadius .LE. 1.0e-05)) THEN
            sizebin = 4
          ELSE IF ((InParticle%EffectiveRadius .GT. 1.0e-05) .AND. &
                   (InParticle%EffectiveRadius .LE. 1.59e-05)) THEN
            sizebin = 5
          ELSE IF ((InParticle%EffectiveRadius .GT. 1.59e-05) .AND. &
                   (InParticle%EffectiveRadius .LE. 2.52e-05)) THEN
            sizebin = 6
          ELSE IF ((InParticle%EffectiveRadius .GT. 2.52e-05) .AND. &
                   (InParticle%EffectiveRadius .LE. 4.01e-05)) THEN
            sizebin = 7
          ELSE IF ((InParticle%EffectiveRadius .GT. 4.01e-05) .AND. &
                   (InParticle%EffectiveRadius .LE. 6.36e-05)) THEN
            sizebin = 8
          ELSE IF ((InParticle%EffectiveRadius .GT. 6.36e-05) .AND. &
                   (InParticle%EffectiveRadius .LE. 1.01e-04)) THEN
            sizebin = 9
          ELSE IF ((InParticle%EffectiveRadius .GT. 1.01e-04) .AND. &
                   (InParticle%EffectiveRadius .LE. 1.0e5)) THEN
            sizebin = 10
          ELSE
            WRITE(*,*) "InParticle%EffectiveRadius = ",InParticle%EffectiveRadius
            call error("Unknown InParticle%EffectiveRadius range. STOP")
          END IF !InParticle%EffectiveRadius

          ! Find core radius ratio 
          radratio = InParticle%AbsCoreRad/InParticle%EffectiveRadius
          IF ((radratio .GT. 1e-09) .AND. (radratio .LE. 0.2)) THEN
            ratiobin = 1
          ELSE IF ((radratio .GT. 0.2) .AND. (radratio .LE. 0.4)) THEN
            ratiobin = 2
          ELSE IF ((radratio .GT. 0.4) .AND. (radratio .LE. 0.6)) THEN
            ratiobin = 3
          ELSE IF ((radratio .GT. 0.6) .AND. (radratio .LE. 0.8)) THEN
            ratiobin = 4
          ELSE IF ((radratio .GT. 0.8) .AND. (radratio .LE. 0.99)) THEN
            ratiobin = 5
          ELSE IF ((radratio .GT. 0.99)) THEN
            ratiobin=5
            WRITE(*,*) "Exceeded 0.99, force set"
            WRITE(*,*) "InParticle%AbsCoreRad = ",InParticle%AbsCoreRad
            WRITE(*,*) "InParticle%EffectiveRadius = ",InParticle%EffectiveRadius
            WRITE(*,*) "radratio = ",radratio
          ELSE
            WRITE(*,*) "InParticle%AbsCoreRad = ",InParticle%AbsCoreRad
            WRITE(*,*) "InParticle%EffectiveRadius = ",InParticle%EffectiveRadius
            WRITE(*,*) "radratio = ",radratio
            call error("Unknown radratio range. STOP")
          END IF !radratio          

	  ! DMiLay output lookup table by 
          ! wavelength bin (ib)
	  ! effective radius size bin (sizebin)
          ! core radius ratio (ratiobin)
          IF (ib .eq. 1) THEN !-------------------------------!MM-DEBUG
            IF (sizebin .EQ. 1) THEN
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =   6.338627489052332E-003
                InParticle%ScaEff(ib) =   1.241488234603144E-003
                InParticle%BackScaEff(ib) =   1.433321141549675E-004
                InParticle%AssymParam(ib) =   1.728520368530504E-005
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =   1.248333996193211E-002
                InParticle%ScaEff(ib) =   1.267648706378087E-003
                InParticle%BackScaEff(ib) =   1.463999009108589E-004
                InParticle%AssymParam(ib) =   1.747987158759047E-005
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =   3.545044720169458E-002
                InParticle%ScaEff(ib) =   1.377484024610864E-003
                InParticle%BackScaEff(ib) =   1.592507689480405E-004
                InParticle%AssymParam(ib) =   1.840337212544841E-005
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =   8.452503919162944E-002
                InParticle%ScaEff(ib) =   1.680195179957653E-003
                InParticle%BackScaEff(ib) =   1.944862060809522E-004
                InParticle%AssymParam(ib) =   2.162433372144669E-005
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =   0.162185165001693     
                InParticle%ScaEff(ib) =   2.364855797138718E-003
                InParticle%BackScaEff(ib) =   2.733733945103151E-004
                InParticle%AssymParam(ib) =   3.193401515405426E-005
              END IF !ratiobin
 
            ELSE IF (sizebin .EQ. 2) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =   7.324842246942065E-002
                InParticle%ScaEff(ib) =   5.714947852366278E-002
                InParticle%BackScaEff(ib) =   5.396783139919124E-003
                InParticle%AssymParam(ib) =   5.307740447917354E-003
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =   9.549220396408416E-002
                InParticle%ScaEff(ib) =   5.848915239976296E-002
                InParticle%BackScaEff(ib) =   5.540367277664674E-003
                InParticle%AssymParam(ib) =   5.367213844836445E-003
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =   0.178405027599639     
                InParticle%ScaEff(ib) =   6.393871254897240E-002
                InParticle%BackScaEff(ib) =   6.113377856825348E-003
                InParticle%AssymParam(ib) =   5.651958365539543E-003
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =   0.353384365393497     
                InParticle%ScaEff(ib) =   7.844087674108738E-002
                InParticle%BackScaEff(ib) =   7.567518690627478E-003
                InParticle%AssymParam(ib) =   6.685904844738792E-003
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =   0.625060843079395     
                InParticle%ScaEff(ib) =   0.109054267629784     
                InParticle%BackScaEff(ib) =   1.035389876492209E-002
                InParticle%AssymParam(ib) =   9.995562341826034E-003
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 3) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =   0.351291783830091     
                InParticle%ScaEff(ib) =   0.319748227133345     
                InParticle%BackScaEff(ib) =   1.898839944583400E-002
                InParticle%AssymParam(ib) =   7.789879551901473E-002
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =   0.402493239923126     
                InParticle%ScaEff(ib) =   0.325522247831505     
                InParticle%BackScaEff(ib) =   1.962465127322233E-002
                InParticle%AssymParam(ib) =   7.811687181908000E-002
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =   0.584608494330860     
                InParticle%ScaEff(ib) =   0.347368523258932     
                InParticle%BackScaEff(ib) =   2.174915876805711E-002
                InParticle%AssymParam(ib) =   8.009690764578657E-002
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =   0.938110307991438     
                InParticle%ScaEff(ib) =   0.404399010676221     
                InParticle%BackScaEff(ib) =   2.554968452519550E-002
                InParticle%AssymParam(ib) =   9.233244460421390E-002
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    1.43868408414879     
                InParticle%ScaEff(ib) =   0.515062275233131     
                InParticle%BackScaEff(ib) =   2.879153764666507E-002
                InParticle%AssymParam(ib) =   0.133777979016359 
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 4) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =    1.37342361010497     
                InParticle%ScaEff(ib) =    1.30034385267227     
                InParticle%BackScaEff(ib) =   3.763082625039913E-003
                InParticle%AssymParam(ib) =   0.800965518269872  
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =    1.43928900487113     
                InParticle%ScaEff(ib) =    1.27608250244860     
                InParticle%BackScaEff(ib) =   6.069771997363317E-003
                InParticle%AssymParam(ib) =   0.778293600986693 
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =    1.65044041551035     
                InParticle%ScaEff(ib) =    1.17975775055226     
                InParticle%BackScaEff(ib) =   9.242745474808526E-003
                InParticle%AssymParam(ib) =   0.704765161285258 
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    1.99089397695059     
                InParticle%ScaEff(ib) =    1.06950509863887     
                InParticle%BackScaEff(ib) =   3.094153704386623E-003
                InParticle%AssymParam(ib) =   0.641970780289727 
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    2.40116615784776     
                InParticle%ScaEff(ib) =    1.11349947132658     
                InParticle%BackScaEff(ib) =   4.742024510115556E-004
                InParticle%AssymParam(ib) =   0.665018755374916   
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 5) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =    3.38575015020377     
                InParticle%ScaEff(ib) =    3.25077565269196     
                InParticle%BackScaEff(ib) =   1.508446412793669E-002
                InParticle%AssymParam(ib) =    2.37766688732998  
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =    3.27409021012275     
                InParticle%ScaEff(ib) =    2.98143230924048     
                InParticle%BackScaEff(ib) =   5.598333667586948E-003
                InParticle%AssymParam(ib) =    2.10816040360092 
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =    3.01731868016347     
                InParticle%ScaEff(ib) =    2.34486286815046     
                InParticle%BackScaEff(ib) =   3.183202528173745E-003
                InParticle%AssymParam(ib) =    1.60144535705915 
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    2.70047140243071     
                InParticle%ScaEff(ib) =    1.54161483767586     
                InParticle%BackScaEff(ib) =   1.651629848283973E-003
                InParticle%AssymParam(ib) =    1.14005705699940  
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    2.76456387593862     
                InParticle%ScaEff(ib) =    1.33122951526543     
                InParticle%BackScaEff(ib) =   3.571474194449051E-003
                InParticle%AssymParam(ib) =    1.02423843119691 
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 6) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =    4.18136701800153     
                InParticle%ScaEff(ib) =    3.97025755785278     
                InParticle%BackScaEff(ib) =   0.108606663672993     
                InParticle%AssymParam(ib) =    2.91867188270452 
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =    4.06394415032188     
                InParticle%ScaEff(ib) =    3.67756794759130     
                InParticle%BackScaEff(ib) =   9.350836275077198E-002
                InParticle%AssymParam(ib) =    2.66177385186489     
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =    3.61794215147068     
                InParticle%ScaEff(ib) =    2.85166425306824     
                InParticle%BackScaEff(ib) =   5.061092198825547E-002
                InParticle%AssymParam(ib) =    1.97047089335660     
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    2.84701689988633     
                InParticle%ScaEff(ib) =    1.57140918123535     
                InParticle%BackScaEff(ib) =   1.004723256882462E-002
                InParticle%AssymParam(ib) =    1.21514996089066    
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    2.54452505400219     
                InParticle%ScaEff(ib) =    1.12185830190530     
                InParticle%BackScaEff(ib) =   2.046874128868622E-003
                InParticle%AssymParam(ib) =   0.979975694893448  
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 7) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =    1.87098454155520     
                InParticle%ScaEff(ib) =    1.57262814971042     
                InParticle%BackScaEff(ib) =   0.179455739020470     
                InParticle%AssymParam(ib) =   0.826587691449864     
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =    2.11331735314581     
                InParticle%ScaEff(ib) =    1.64429534097236     
                InParticle%BackScaEff(ib) =   0.184961775133334     
                InParticle%AssymParam(ib) =   0.893850657724618   
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =    2.65893542470647     
                InParticle%ScaEff(ib) =    1.83796372572409     
                InParticle%BackScaEff(ib) =   0.169899454145463     
                InParticle%AssymParam(ib) =    1.21605889053948    
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    2.83157834003745     
                InParticle%ScaEff(ib) =    1.63313473394288     
                InParticle%BackScaEff(ib) =   1.550668915559024E-002
                InParticle%AssymParam(ib) =    1.31003612229362
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    2.42950597324128     
                InParticle%ScaEff(ib) =    1.07532366873380     
                InParticle%BackScaEff(ib) =   7.534829502916470E-004
                InParticle%AssymParam(ib) =    1.01785913956035  
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 8) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =    2.75273381059395     
                InParticle%ScaEff(ib) =    2.29287661020004     
                InParticle%BackScaEff(ib) =   0.332542993634415     
                InParticle%AssymParam(ib) =    1.87727764549639  
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =    2.78522801824329     
                InParticle%ScaEff(ib) =    2.17747302639268     
                InParticle%BackScaEff(ib) =   0.384586401022679     
                InParticle%AssymParam(ib) =    1.76367656922482    
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =    2.32899233368725     
                InParticle%ScaEff(ib) =    1.43454170410489     
                InParticle%BackScaEff(ib) =   0.300743383637634     
                InParticle%AssymParam(ib) =    1.11524168831126  
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    2.38808091024147     
                InParticle%ScaEff(ib) =    1.17829949927142     
                InParticle%BackScaEff(ib) =   0.110810966562384     
                InParticle%AssymParam(ib) =    1.05724712475111   
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    2.41621435374200     
                InParticle%ScaEff(ib) =    1.14849985328176     
                InParticle%BackScaEff(ib) =   4.253146086490113E-003
                InParticle%AssymParam(ib) =    1.07742063648435 
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 9) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =    2.44858113182375     
                InParticle%ScaEff(ib) =    1.89873736543775     
                InParticle%BackScaEff(ib) =   1.756280136843464E-002
                InParticle%AssymParam(ib) =    1.55333735008986 
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =    2.55761296962624     
                InParticle%ScaEff(ib) =    1.89816188211727     
                InParticle%BackScaEff(ib) =   6.173942284503940E-002
                InParticle%AssymParam(ib) =    1.55863693680222   
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =    2.11680519282049     
                InParticle%ScaEff(ib) =    1.24821332222136     
                InParticle%BackScaEff(ib) =   8.819402620649464E-002
                InParticle%AssymParam(ib) =   0.960133085076927 
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    2.41818605856077     
                InParticle%ScaEff(ib) =    1.32259918890826     
                InParticle%BackScaEff(ib) =   4.781133523560435E-003
                InParticle%AssymParam(ib) =    1.23957418310785 
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    2.36876691735927     
                InParticle%ScaEff(ib) =    1.28176946657756     
                InParticle%BackScaEff(ib) =   9.201960147680920E-003
                InParticle%AssymParam(ib) =    1.13906085474838 
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 10) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =    2.13297337686945     
                InParticle%ScaEff(ib) =    1.35787060904386     
                InParticle%BackScaEff(ib) =   1.548551042971653E-003
                InParticle%AssymParam(ib) =    1.21616908455416 
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =    2.20928800758987     
                InParticle%ScaEff(ib) =    1.39399897848554     
                InParticle%BackScaEff(ib) =   8.294012329455452E-003
                InParticle%AssymParam(ib) =    1.25407274982272 
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =    2.22156244494289     
                InParticle%ScaEff(ib) =    1.32399449175411     
                InParticle%BackScaEff(ib) =   6.473783882947021E-003
                InParticle%AssymParam(ib) =    1.19914631263325 
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    2.18138867980128     
                InParticle%ScaEff(ib) =    1.19396817796879     
                InParticle%BackScaEff(ib) =   4.716963136782338E-003
                InParticle%AssymParam(ib) =    1.11352663602915   
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    2.18569477923327     
                InParticle%ScaEff(ib) =    1.20660859994311     
                InParticle%BackScaEff(ib) =   3.016967617235626E-003
                InParticle%AssymParam(ib) =    1.10586310013231 
              END IF !ratiobin
            END IF !sizebin
 
          ELSE IF (ib .eq. 2) THEN !-------------------------------!MM-DEBUG
            IF (sizebin .EQ. 1) THEN
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =   5.139038309134193E-003
                InParticle%ScaEff(ib) =   1.062224852832295E-003
                InParticle%BackScaEff(ib) =   1.229431498678993E-004
                InParticle%AssymParam(ib) =   1.368792181601226E-005
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =   1.103167815585350E-002
                InParticle%ScaEff(ib) =   1.085095451234660E-003
                InParticle%BackScaEff(ib) =   1.256286612884928E-004
                InParticle%AssymParam(ib) =   1.384556118672044E-005
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =   3.304657208718486E-002
                InParticle%ScaEff(ib) =   1.181018035298606E-003
                InParticle%BackScaEff(ib) =   1.368679159545767E-004
                InParticle%AssymParam(ib) =   1.459447195855821E-005
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =   8.003271280431039E-002
                InParticle%ScaEff(ib) =   1.444726916079719E-003
                InParticle%BackScaEff(ib) =   1.676181042182596E-004
                InParticle%AssymParam(ib) =   1.720509234474033E-005
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =   0.154227758362364     
                InParticle%ScaEff(ib) =   2.039345294924047E-003
                InParticle%BackScaEff(ib) =   2.363008157370650E-004
                InParticle%AssymParam(ib) =   2.554250517508710E-005
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 2) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
		InParticle%ExtEff(ib) =   6.174351265724396E-002
                InParticle%ScaEff(ib) =   4.900113880311652E-002
                InParticle%BackScaEff(ib) =   4.714337877280814E-003
                InParticle%AssymParam(ib) =   4.209664958810828E-003
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =   8.264246166699793E-002
                InParticle%ScaEff(ib) =   5.017934629863428E-002
                InParticle%BackScaEff(ib) =   4.841424949812438E-003
                InParticle%AssymParam(ib) =   4.258892657080572E-003
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =   0.160648446320607     
                InParticle%ScaEff(ib) =   5.497998991483244E-002
                InParticle%BackScaEff(ib) =   5.350354387040322E-003
                InParticle%AssymParam(ib) =   4.493809289003795E-003
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =   0.325671135931675     
                InParticle%ScaEff(ib) =   6.775085956735458E-002
                InParticle%BackScaEff(ib) =   6.647893983457620E-003
                InParticle%AssymParam(ib) =   5.338322497605263E-003
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =   0.582615408336245     
                InParticle%ScaEff(ib) =   9.476760901990416E-002
                InParticle%BackScaEff(ib) =   9.165054441082652E-003
                InParticle%AssymParam(ib) =   8.028235353939209E-003
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 3) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =   0.305464876517266     
                InParticle%ScaEff(ib) =   0.280575984311160     
                InParticle%BackScaEff(ib) =   1.789571418393630E-002
                InParticle%AssymParam(ib) =   6.261560271941657E-002
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =   0.353870970813702     
                InParticle%ScaEff(ib) =   0.286290024924826     
                InParticle%BackScaEff(ib) =   1.850049093473012E-002
                InParticle%AssymParam(ib) =   6.291873757478936E-002
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =   0.527508642710369     
                InParticle%ScaEff(ib) =   0.308021878809883     
                InParticle%BackScaEff(ib) =   2.058334313826214E-002
                InParticle%AssymParam(ib) =   6.496036939386565E-002
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =   0.868170462402061     
                InParticle%ScaEff(ib) =   0.363582058736275     
                InParticle%BackScaEff(ib) =   2.456105876110214E-002
                InParticle%AssymParam(ib) =   7.564743502150444E-002
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    1.35450097228911     
                InParticle%ScaEff(ib) =   0.470022020429068     
                InParticle%BackScaEff(ib) =   2.864258940056085E-002
                InParticle%AssymParam(ib) =   0.111199968822955 
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 4) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =    1.20301001368349     
                InParticle%ScaEff(ib) =    1.14643428189651     
                InParticle%BackScaEff(ib) =   3.728795170429543E-003
                InParticle%AssymParam(ib) =   0.683668626620534     
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =    1.27324128141990     
                InParticle%ScaEff(ib) =    1.13052692814481     
                InParticle%BackScaEff(ib) =   5.953346666972397E-003
                InParticle%AssymParam(ib) =   0.666327205880667 
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =    1.50912825191453     
                InParticle%ScaEff(ib) =    1.07089188764741     
                InParticle%BackScaEff(ib) =   9.428179277320462E-003
                InParticle%AssymParam(ib) =   0.615974866884009  
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    1.91045026916680     
                InParticle%ScaEff(ib) =    1.01633546416832     
                InParticle%BackScaEff(ib) =   4.507098562386480E-003
                InParticle%AssymParam(ib) =   0.587384510776133 
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    2.36268918025467     
                InParticle%ScaEff(ib) =    1.08522888394559     
                InParticle%BackScaEff(ib) =   1.196208984032071E-003
                InParticle%AssymParam(ib) =   0.626970551154550 
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 5) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =    3.17683312725760     
                InParticle%ScaEff(ib) =    3.06691675151759     
                InParticle%BackScaEff(ib) =   1.638071099334255E-002
                InParticle%AssymParam(ib) =    2.21469171886986     
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =    3.10316297314846     
                InParticle%ScaEff(ib) =    2.83578663035769     
                InParticle%BackScaEff(ib) =   6.341459097105847E-003
                InParticle%AssymParam(ib) =    1.99640121691916 
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =    2.92370396649955     
                InParticle%ScaEff(ib) =    2.26842320042052     
                InParticle%BackScaEff(ib) =   3.443317796962763E-003
                InParticle%AssymParam(ib) =    1.55074999195171 
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    2.70807959371097     
                InParticle%ScaEff(ib) =    1.54648151179000     
                InParticle%BackScaEff(ib) =   2.312392736724962E-003
                InParticle%AssymParam(ib) =    1.13771473209496 
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    2.77208430333744     
                InParticle%ScaEff(ib) =    1.34549800213284     
                InParticle%BackScaEff(ib) =   6.126951669554335E-003
                InParticle%AssymParam(ib) =    1.01984374784182
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 6) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =    4.23984112710875     
                InParticle%ScaEff(ib) =    4.05787742959576     
                InParticle%BackScaEff(ib) =   0.146014487633035     
                InParticle%AssymParam(ib) =    2.95175309258135     
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =    4.08513798438302     
                InParticle%ScaEff(ib) =    3.72078037219820     
                InParticle%BackScaEff(ib) =   0.125042913708215     
                InParticle%AssymParam(ib) =    2.63671204196320  
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =    3.61430185207521     
                InParticle%ScaEff(ib) =    2.86519233639460     
                InParticle%BackScaEff(ib) =   7.132297716563619E-002
                InParticle%AssymParam(ib) =    1.92704590242080  
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    2.89251558294860     
                InParticle%ScaEff(ib) =    1.58849906905663     
                InParticle%BackScaEff(ib) =   1.408075881283735E-002
                InParticle%AssymParam(ib) =    1.22953821185668  
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    2.56567682815404     
                InParticle%ScaEff(ib) =    1.13796710201124     
                InParticle%BackScaEff(ib) =   4.484374058573750E-004
                InParticle%AssymParam(ib) =   0.984674160656133  
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 7) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =    2.13135989024255     
                InParticle%ScaEff(ib) =    1.85793408755914     
                InParticle%BackScaEff(ib) =   0.241193857011960     
                InParticle%AssymParam(ib) =   0.971547494130008     
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =    2.40289721652576     
                InParticle%ScaEff(ib) =    1.94844833740984     
                InParticle%BackScaEff(ib) =   0.269022640075835     
                InParticle%AssymParam(ib) =    1.07243353795922 
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =    2.98625319288557     
                InParticle%ScaEff(ib) =    2.15647180032080     
                InParticle%BackScaEff(ib) =   0.216819365604938     
                InParticle%AssymParam(ib) =    1.41424419435123 
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    2.98413938438816     
                InParticle%ScaEff(ib) =    1.72378363078381     
                InParticle%BackScaEff(ib) =   2.513497377169763E-003
                InParticle%AssymParam(ib) =    1.38572260018926  
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    2.43456544414595     
                InParticle%ScaEff(ib) =    1.07471266757135     
                InParticle%BackScaEff(ib) =   9.341586908463232E-004
                InParticle%AssymParam(ib) =    1.01384372016860   
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 8) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =    2.98515002132892     
                InParticle%ScaEff(ib) =    2.58582574596983     
                InParticle%BackScaEff(ib) =   0.215299572699643     
                InParticle%AssymParam(ib) =    2.03039850011861     
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =    2.88156975776157     
                InParticle%ScaEff(ib) =    2.31574174515928     
                InParticle%BackScaEff(ib) =   0.265991554779424     
                InParticle%AssymParam(ib) =    1.78625929718757     
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =    2.33874988836856     
                InParticle%ScaEff(ib) =    1.46334441697541     
                InParticle%BackScaEff(ib) =   0.215984644393330     
                InParticle%AssymParam(ib) =    1.03370071092364     
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    2.41403552508252     
                InParticle%ScaEff(ib) =    1.22434799971637     
                InParticle%BackScaEff(ib) =   2.522146742431766E-002
                InParticle%AssymParam(ib) =    1.11515601761215 
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    2.41343216825036     
                InParticle%ScaEff(ib) =    1.13517179373856     
                InParticle%BackScaEff(ib) =   3.398463151616537E-003
                InParticle%AssymParam(ib) =    1.07000209174110 
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 9) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =    2.57288363479180     
                InParticle%ScaEff(ib) =    2.10963425218824     
                InParticle%BackScaEff(ib) =   1.204908203279876E-002
                InParticle%AssymParam(ib) =    1.73101183532299     
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =    2.49352858893840     
                InParticle%ScaEff(ib) =    1.90095836014293     
                InParticle%BackScaEff(ib) =   5.480643006265104E-002
                InParticle%AssymParam(ib) =    1.53622154748801    
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =    2.09685865156925     
                InParticle%ScaEff(ib) =    1.27034787442852     
                InParticle%BackScaEff(ib) =   8.413073035068880E-002
                InParticle%AssymParam(ib) =   0.957746167416338  
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    2.38692448995987     
                InParticle%ScaEff(ib) =    1.29576188661496     
                InParticle%BackScaEff(ib) =   1.134166856340495E-003
                InParticle%AssymParam(ib) =    1.20948776303605 
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    2.38285108786522     
                InParticle%ScaEff(ib) =    1.28014803616328     
                InParticle%BackScaEff(ib) =   8.369357332960934E-003
                InParticle%AssymParam(ib) =    1.13828009609778  
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 10) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =    2.23764411470239     
                InParticle%ScaEff(ib) =    1.53567895609872     
                InParticle%BackScaEff(ib) =   3.392344224483384E-002
                InParticle%AssymParam(ib) =    1.36350832165984     
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =    2.09578050154409     
                InParticle%ScaEff(ib) =    1.33898693533697     
                InParticle%BackScaEff(ib) =   4.483088457007030E-002
                InParticle%AssymParam(ib) =    1.16858182418855  
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =    2.10355282540924     
                InParticle%ScaEff(ib) =    1.23516698011055     
                InParticle%BackScaEff(ib) =   4.588265698179459E-002
                InParticle%AssymParam(ib) =    1.08635277102341  
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    2.24895160836331     
                InParticle%ScaEff(ib) =    1.26160576658518     
                InParticle%BackScaEff(ib) =   4.644847274671248E-003
                InParticle%AssymParam(ib) =    1.17875993716182   
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    2.18383476886376     
                InParticle%ScaEff(ib) =    1.20213075138530     
                InParticle%BackScaEff(ib) =   5.250066425087228E-003
                InParticle%AssymParam(ib) =    1.08774762885791  
              END IF !ratiobin
            END IF !sizebin

          ELSE IF (ib .eq. 3) THEN !-------------------------------!MM-DEBUG
            IF (sizebin .EQ. 1) THEN
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =   4.263139559487195E-003
                InParticle%ScaEff(ib) =   8.678188365096641E-004
                InParticle%BackScaEff(ib) =   1.007407481705718E-004
                InParticle%AssymParam(ib) =   1.011486858243403E-005
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =   9.824738790186843E-003
                InParticle%ScaEff(ib) =   8.865980196457613E-004
                InParticle%BackScaEff(ib) =   1.029492819769158E-004
                InParticle%AssymParam(ib) =   1.023204761827078E-005
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =   3.059731242702119E-002
                InParticle%ScaEff(ib) =   9.653255345575724E-004
                InParticle%BackScaEff(ib) =   1.121898933991893E-004
                InParticle%AssymParam(ib) =   1.078930308842013E-005
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =   7.490825388297630E-002
                InParticle%ScaEff(ib) =   1.181492663995140E-003
                InParticle%BackScaEff(ib) =   1.374515268544293E-004
                InParticle%AssymParam(ib) =   1.273109799514702E-005
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =   0.144812679292337     
                InParticle%ScaEff(ib) =   1.668311505157806E-003
                InParticle%BackScaEff(ib) =   1.938561160319200E-004
                InParticle%AssymParam(ib) =   1.892483759965921E-005
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 2) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =   5.057006973339805E-002
                InParticle%ScaEff(ib) =   4.010205846162691E-002
                InParticle%BackScaEff(ib) =   3.943011051862751E-003
                InParticle%AssymParam(ib) =   3.115401483151503E-003
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =   6.977025758065473E-002
                InParticle%ScaEff(ib) =   4.107031545316867E-002
                InParticle%BackScaEff(ib) =   4.048318605115704E-003
                InParticle%AssymParam(ib) =   3.152597047211763E-003
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =   0.141512475631317     
                InParticle%ScaEff(ib) =   4.502610859996364E-002
                InParticle%BackScaEff(ib) =   4.472038199127779E-003
                InParticle%AssymParam(ib) =   3.329597172343329E-003
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =   0.293670507177789     
                InParticle%ScaEff(ib) =   5.557005162808929E-002
                InParticle%BackScaEff(ib) =   5.560477419200620E-003
                InParticle%AssymParam(ib) =   3.960465083758706E-003
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =   0.531433841467238     
                InParticle%ScaEff(ib) =   7.800226375329375E-002
                InParticle%BackScaEff(ib) =   7.708573944335241E-003
                InParticle%AssymParam(ib) =   5.963814556042968E-003
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 3) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =   0.255624221481006     
                InParticle%ScaEff(ib) =   0.235302311794726     
                InParticle%BackScaEff(ib) =   1.622437278838499E-002
                InParticle%AssymParam(ib) =   4.699747032971522E-002
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =   0.299901356490615     
                InParticle%ScaEff(ib) =   0.240498741906031     
                InParticle%BackScaEff(ib) =   1.676200456901536E-002
                InParticle%AssymParam(ib) =   4.731412025801772E-002
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =   0.460186402890855     
                InParticle%ScaEff(ib) =   0.260466841879183     
                InParticle%BackScaEff(ib) =   1.867559148941026E-002
                InParticle%AssymParam(ib) =   4.915210649956857E-002
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =   0.779050110889278     
                InParticle%ScaEff(ib) =   0.311184486362578     
                InParticle%BackScaEff(ib) =   2.258814091761381E-002
                InParticle%AssymParam(ib) =   5.766194544095191E-002
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    1.23996624514364     
                InParticle%ScaEff(ib) =   0.408463054261145     
                InParticle%BackScaEff(ib) =   2.734083581546594E-002
                InParticle%AssymParam(ib) =   8.561118615865378E-002
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 4) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =    1.01940542275800     
                InParticle%ScaEff(ib) =   0.974707004521188     
                InParticle%BackScaEff(ib) =   6.393422006093996E-003
                InParticle%AssymParam(ib) =   0.543270423557770     
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =    1.09150328194008     
                InParticle%ScaEff(ib) =   0.966509038146904     
                InParticle%BackScaEff(ib) =   8.311286263230834E-003
                InParticle%AssymParam(ib) =   0.531535831562413 
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =    1.33852810361431     
                InParticle%ScaEff(ib) =   0.938731322831333     
                InParticle%BackScaEff(ib) =   1.162230282168859E-002
                InParticle%AssymParam(ib) =   0.502398103621807 
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    1.78429194952155     
                InParticle%ScaEff(ib) =   0.936569432400333     
                InParticle%BackScaEff(ib) =   8.035549488415072E-003
                InParticle%AssymParam(ib) =   0.506131319667848 
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    2.29222835390614     
                InParticle%ScaEff(ib) =    1.03395440348948     
                InParticle%BackScaEff(ib) =   3.580206860432335E-003
                InParticle%AssymParam(ib) =   0.568648666960359    
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 5) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =    2.82373243492915     
                InParticle%ScaEff(ib) =    2.73627678098817     
                InParticle%BackScaEff(ib) =   2.762871640954166E-002
                InParticle%AssymParam(ib) =    1.90456440796543     
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =    2.79171757067844     
                InParticle%ScaEff(ib) =    2.55402760711118     
                InParticle%BackScaEff(ib) =   1.875195187120392E-002
                InParticle%AssymParam(ib) =    1.75134734125139 
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =    2.70058090013906     
                InParticle%ScaEff(ib) =    2.08043513665998     
                InParticle%BackScaEff(ib) =   1.294640567268690E-002
                InParticle%AssymParam(ib) =    1.40223289260488   
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    2.67205845603090     
                InParticle%ScaEff(ib) =    1.52431176163280     
                InParticle%BackScaEff(ib) =   5.801715152035907E-003
                InParticle%AssymParam(ib) =    1.10659639703817 
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    2.78628380135758     
                InParticle%ScaEff(ib) =    1.36074130619291     
                InParticle%BackScaEff(ib) =   9.037377844402734E-003
                InParticle%AssymParam(ib) =    1.01355098208938 
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 6) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =    4.02361663987596     
                InParticle%ScaEff(ib) =    3.88486261377668     
                InParticle%BackScaEff(ib) =   8.836415432823291E-002
                InParticle%AssymParam(ib) =    2.84931045496053     
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =    3.81636850766442     
                InParticle%ScaEff(ib) =    3.49034762729065     
                InParticle%BackScaEff(ib) =   7.591212543258335E-002
                InParticle%AssymParam(ib) =    2.45544094723994 
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =    3.30904943558053     
                InParticle%ScaEff(ib) =    2.60632644752455     
                InParticle%BackScaEff(ib) =   3.168884184523759E-002
                InParticle%AssymParam(ib) =    1.69868864562301  
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    2.81698823339211     
                InParticle%ScaEff(ib) =    1.54303889094819     
                InParticle%BackScaEff(ib) =   8.778628905193386E-003
                InParticle%AssymParam(ib) =    1.17769066859975 
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    2.60591462359312     
                InParticle%ScaEff(ib) =    1.16408535385012     
                InParticle%BackScaEff(ib) =   1.110688542706662E-003
                InParticle%AssymParam(ib) =   0.995936288016512  
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 7) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =    2.13140593071492     
                InParticle%ScaEff(ib) =    1.91090021632114     
                InParticle%BackScaEff(ib) =   0.260019140372145     
                InParticle%AssymParam(ib) =   0.952632522132422     
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =    2.40423282417051     
                InParticle%ScaEff(ib) =    2.00506467553308     
                InParticle%BackScaEff(ib) =   0.245418726819286     
                InParticle%AssymParam(ib) =    1.05361578839653 
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =    2.88402673625340     
                InParticle%ScaEff(ib) =    2.08256710393457     
                InParticle%BackScaEff(ib) =   0.134760983236147     
                InParticle%AssymParam(ib) =    1.29297896211798   
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    2.90341053850232     
                InParticle%ScaEff(ib) =    1.66754570953437     
                InParticle%BackScaEff(ib) =   3.716049871002671E-002
                InParticle%AssymParam(ib) =    1.28930458133611 
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    2.44319477905224     
                InParticle%ScaEff(ib) =    1.07440846412043     
                InParticle%BackScaEff(ib) =   1.180109235404892E-004
                InParticle%AssymParam(ib) =    1.00793025770467  
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 8) THEN !!!!!!!!!!!!!!!!!!!!!
               IF (ratiobin .EQ. 1) THEN
                 InParticle%ExtEff(ib) =    2.81760720579818     
                 InParticle%ScaEff(ib) =    2.51729276797066     
                 InParticle%BackScaEff(ib) =   0.197318430580763     
                 InParticle%AssymParam(ib) =    2.04044574055271     
              ELSE IF (ratiobin .EQ. 2) THEN
                 InParticle%ExtEff(ib) =    2.57488480960150     
                 InParticle%ScaEff(ib) =    2.10668718533517     
                 InParticle%BackScaEff(ib) =   0.201193345762590     
                 InParticle%AssymParam(ib) =    1.63060098489225 
              ELSE IF (ratiobin .EQ. 3) THEN
                 InParticle%ExtEff(ib) =    2.04767334628316     
                 InParticle%ScaEff(ib) =    1.24793885026032     
                 InParticle%BackScaEff(ib) =   0.120656751271596     
                 InParticle%AssymParam(ib) =   0.890123772127963 
              ELSE IF (ratiobin .EQ. 4) THEN
                 InParticle%ExtEff(ib) =    2.45837383604067     
                 InParticle%ScaEff(ib) =    1.28060031799697     
                 InParticle%BackScaEff(ib) =   5.244784741941551E-002
                 InParticle%AssymParam(ib) =    1.14765636485615 
              ELSE IF (ratiobin .EQ. 5) THEN
                 InParticle%ExtEff(ib) =    2.41085234555214     
                 InParticle%ScaEff(ib) =    1.12036950583058     
                 InParticle%BackScaEff(ib) =   3.010789121369003E-003
                 InParticle%AssymParam(ib) =    1.06140943853389 
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 9) THEN !!!!!!!!!!!!!!!!!!!!!
               IF (ratiobin .EQ. 1) THEN
                 InParticle%ExtEff(ib) =    2.60877262142503     
                 InParticle%ScaEff(ib) =    2.19230972755405     
                 InParticle%BackScaEff(ib) =   3.817992804076797E-002
                 InParticle%AssymParam(ib) =    1.77911080621722     
              ELSE IF (ratiobin .EQ. 2) THEN
                 InParticle%ExtEff(ib) =    2.33891354601800     
                 InParticle%ScaEff(ib) =    1.78239959331719     
                 InParticle%BackScaEff(ib) =   9.912359476817148E-002
                 InParticle%AssymParam(ib) =    1.37900071570137   
              ELSE IF (ratiobin .EQ. 3) THEN
                 InParticle%ExtEff(ib) =    2.21209387722090     
                 InParticle%ScaEff(ib) =    1.39426224157774     
                 InParticle%BackScaEff(ib) =   0.130532300116988     
                 InParticle%AssymParam(ib) =    1.05801687815812  
              ELSE IF (ratiobin .EQ. 4) THEN
                 InParticle%ExtEff(ib) =    2.32551374554615     
                 InParticle%ScaEff(ib) =    1.24141997828906     
                 InParticle%BackScaEff(ib) =   1.670943852470345E-003
                 InParticle%AssymParam(ib) =    1.14203276923661 
              ELSE IF (ratiobin .EQ. 5) THEN
                 InParticle%ExtEff(ib) =    2.39805184963928     
                 InParticle%ScaEff(ib) =    1.27305434602451     
                 InParticle%BackScaEff(ib) =   9.189984740622667E-003
                 InParticle%AssymParam(ib) =    1.13577122111321 
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 10) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =    2.14964998564504     
                InParticle%ScaEff(ib) =    1.51295303925601     
                InParticle%BackScaEff(ib) =   5.334558397666915E-002
                InParticle%AssymParam(ib) =    1.31845874388974     
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =    2.13670642026765     
                InParticle%ScaEff(ib) =    1.43234784607461     
                InParticle%BackScaEff(ib) =   3.581665921291712E-002
                InParticle%AssymParam(ib) =    1.24015288508466  
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =    2.07566598312208     
                InParticle%ScaEff(ib) =    1.23342500902051     
                InParticle%BackScaEff(ib) =   5.009895263700193E-002
                InParticle%AssymParam(ib) =    1.07244466563173 
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    2.22863863650280     
                InParticle%ScaEff(ib) =    1.23323597862107     
                InParticle%BackScaEff(ib) =   9.968545595441249E-003
                InParticle%AssymParam(ib) =    1.1632950291983
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    2.17256371810077     
                InParticle%ScaEff(ib) =    1.18063013582568     
                InParticle%BackScaEff(ib) =   7.661975129484521E-003
                InParticle%AssymParam(ib) =    1.06154342330829  
              END IF !ratiobin
            END IF !sizebin

          ELSE IF (ib .eq. 4) THEN !-------------------------------!MM-DEBUG
            IF (sizebin .EQ. 1) THEN
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =   2.753364167085392E-003
                InParticle%ScaEff(ib) =   4.725331281362864E-004
                InParticle%BackScaEff(ib) =   5.525631584254712E-005
                InParticle%AssymParam(ib) =   4.071316893356675E-006
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =   7.425505208959161E-003
                InParticle%ScaEff(ib) =   4.827810529546341E-004
                InParticle%BackScaEff(ib) =   5.646618826504992E-005
                InParticle%AssymParam(ib) =   4.118701937006665E-006
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =   2.486664441636313E-002
                InParticle%ScaEff(ib) =   5.256939654885124E-004
                InParticle%BackScaEff(ib) =   6.152496371390335E-005
                InParticle%AssymParam(ib) =   4.344494468088785E-006
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =   6.203756948746719E-002
                InParticle%ScaEff(ib) =   6.431720886111547E-004
                InParticle%BackScaEff(ib) =   7.532868641386177E-005
                InParticle%AssymParam(ib) =   5.130076429292264E-006
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =   0.120590912593262     
                InParticle%ScaEff(ib) =   9.070519208908492E-004
                InParticle%BackScaEff(ib) =   1.061388237443969E-004
                InParticle%AssymParam(ib) =   7.627445733771990E-006
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 2) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =   2.864233227498378E-002
                InParticle%ScaEff(ib) =   2.186962766642239E-002
                InParticle%BackScaEff(ib) =   2.265966881205732E-003
                InParticle%AssymParam(ib) =   1.257057690430274E-003
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =   4.366756379127423E-002
                InParticle%ScaEff(ib) =   2.239247385477439E-002
                InParticle%BackScaEff(ib) =   2.324102770485585E-003
                InParticle%AssymParam(ib) =   1.272431377629956E-003
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =   9.988853168620168E-002
                InParticle%ScaEff(ib) =   2.454187372940983E-002
                InParticle%BackScaEff(ib) =   2.560548197627754E-003
                InParticle%AssymParam(ib) =   1.345304400266801E-003
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =   0.219654871327424     
                InParticle%ScaEff(ib) =   3.030065847186051E-002
                InParticle%BackScaEff(ib) =   3.178346543821798E-003
                InParticle%AssymParam(ib) =   1.600874676180828E-003
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =   0.408067820596703     
                InParticle%ScaEff(ib) =   4.272803842922383E-002
                InParticle%BackScaEff(ib) =   4.445949619160706E-003
                InParticle%AssymParam(ib) =   2.407762067481284E-003
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 3) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =   0.147288160966071     
                InParticle%ScaEff(ib) =   0.134465686275423     
                InParticle%BackScaEff(ib) =   1.097615308225398E-002
                InParticle%AssymParam(ib) =   1.948109916119950E-002
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =   0.180354610625155     
                InParticle%ScaEff(ib) =   0.137771689092346     
                InParticle%BackScaEff(ib) =   1.131710925779624E-002
                InParticle%AssymParam(ib) =   1.968055173124623E-002
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =   0.302347408115515     
                InParticle%ScaEff(ib) =   0.150817512007478     
                InParticle%BackScaEff(ib) =   1.260948478125482E-002
                InParticle%AssymParam(ib) =   2.068048504157519E-002
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =   0.553448056892250     
                InParticle%ScaEff(ib) =   0.184252437596009     
                InParticle%BackScaEff(ib) =   1.559196483684770E-002
                InParticle%AssymParam(ib) =   2.456498493945252E-002
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =   0.929974307432744     
                InParticle%ScaEff(ib) =   0.250380937135547     
                InParticle%BackScaEff(ib) =   2.033968194994756E-002
                InParticle%AssymParam(ib) =   3.699096772651982E-002
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 4) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =   0.652413758309375     
                InParticle%ScaEff(ib) =   0.626191688703548     
                InParticle%BackScaEff(ib) =   1.778549309279416E-002
                InParticle%AssymParam(ib) =   0.253968122802817     
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =   0.721025076713686     
                InParticle%ScaEff(ib) =   0.630468248902443     
                InParticle%BackScaEff(ib) =   1.893515198426141E-002
                InParticle%AssymParam(ib) =   0.251574796453445 
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =   0.956324844257156     
                InParticle%ScaEff(ib) =   0.646146505558791     
                InParticle%BackScaEff(ib) =   2.162672347371767E-002
                InParticle%AssymParam(ib) =   0.248984659038393     
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    1.39820340055190     
                InParticle%ScaEff(ib) =   0.705373479886065     
                InParticle%BackScaEff(ib) =   2.214629728784083E-002
                InParticle%AssymParam(ib) =   0.277005623997083 
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    1.98374226343678     
                InParticle%ScaEff(ib) =   0.840267636255933     
                InParticle%BackScaEff(ib) =   1.760556340296673E-002
                InParticle%AssymParam(ib) =   0.366428447602450  
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 5) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =    2.21481592015604     
                InParticle%ScaEff(ib) =    2.16240977794318     
                InParticle%BackScaEff(ib) =   3.898035959863895E-002
                InParticle%AssymParam(ib) =    1.33502856826726     
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =    2.24271052456984     
                InParticle%ScaEff(ib) =    2.06800362941874     
                InParticle%BackScaEff(ib) =   3.629342249738617E-002
                InParticle%AssymParam(ib) =    1.27799769373332  
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =    2.24314254153314     
                InParticle%ScaEff(ib) =    1.70972523872547     
                InParticle%BackScaEff(ib) =   2.928371976972988E-002
                InParticle%AssymParam(ib) =    1.06585367139569    
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    2.36865146176109     
                InParticle%ScaEff(ib) =    1.35898568682125     
                InParticle%BackScaEff(ib) =   1.235527109455989E-002
                InParticle%AssymParam(ib) =   0.906115632551481  
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    2.73420713362149     
                InParticle%ScaEff(ib) =    1.33278582460604     
                InParticle%BackScaEff(ib) =   9.967977942811165E-003
                InParticle%AssymParam(ib) =   0.925948717186274  
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 6) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =    4.13294149320973     
                InParticle%ScaEff(ib) =    4.03032111503355     
                InParticle%BackScaEff(ib) =   7.360828078472970E-002
                InParticle%AssymParam(ib) =    2.95878817637063     
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =    3.93657977690474     
                InParticle%ScaEff(ib) =    3.66770358783291     
                InParticle%BackScaEff(ib) =   6.378862352038903E-002
                InParticle%AssymParam(ib) =    2.61411380470694 
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =    3.45846369760482     
                InParticle%ScaEff(ib) =    2.76610072707306     
                InParticle%BackScaEff(ib) =   2.661771858894812E-002
                InParticle%AssymParam(ib) =    1.88735403979881
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    2.86381723273350     
                InParticle%ScaEff(ib) =    1.60383375935866     
                InParticle%BackScaEff(ib) =   5.473898240651396E-003
                InParticle%AssymParam(ib) =    1.22698810833429  
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    2.69304955122330     
                InParticle%ScaEff(ib) =    1.25033059908016     
                InParticle%BackScaEff(ib) =   6.660793655015043E-003
                InParticle%AssymParam(ib) =    1.02464397830354   
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 7) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =    2.95525036302991     
                InParticle%ScaEff(ib) =    2.81071944848402     
                InParticle%BackScaEff(ib) =   0.144722870003543     
                InParticle%AssymParam(ib) =    1.79750526281799     
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =    3.07622075305655     
                InParticle%ScaEff(ib) =    2.74767979573311     
                InParticle%BackScaEff(ib) =   0.114086144198558     
                InParticle%AssymParam(ib) =    1.75771823060325 
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =    3.13443353727213     
                InParticle%ScaEff(ib) =    2.38962795526458     
                InParticle%BackScaEff(ib) =   6.656313299481739E-002
                InParticle%AssymParam(ib) =    1.56192006850129  
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    2.81932474143716     
                InParticle%ScaEff(ib) =    1.59110473436932     
                InParticle%BackScaEff(ib) =   1.192733327901523E-002
                InParticle%AssymParam(ib) =    1.21178915781235  
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    2.47216837452084     
                InParticle%ScaEff(ib) =    1.07757752023344     
                InParticle%BackScaEff(ib) =   1.080634786921352E-003
                InParticle%AssymParam(ib) =   0.987493204096343  
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 8) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =    2.27724011917898     
                InParticle%ScaEff(ib) =    2.03534636405276     
                InParticle%BackScaEff(ib) =   0.296937265639801     
                InParticle%AssymParam(ib) =    1.33451093586322     
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =    2.06821111017703     
                InParticle%ScaEff(ib) =    1.64151843375264     
                InParticle%BackScaEff(ib) =   0.270707651280783     
                InParticle%AssymParam(ib) =   0.953114456017294 
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =    2.07020163745579     
                InParticle%ScaEff(ib) =    1.29124608353427     
                InParticle%BackScaEff(ib) =   0.170095675320752     
                InParticle%AssymParam(ib) =   0.702221463350847 
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    2.72568302842573     
                InParticle%ScaEff(ib) =    1.51534078721248     
                InParticle%BackScaEff(ib) =   4.813809898166456E-002
                InParticle%AssymParam(ib) =    1.25743137980625 
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    2.41019470324832     
                InParticle%ScaEff(ib) =    1.09009441391071     
                InParticle%BackScaEff(ib) =   1.305040745713129E-003
                InParticle%AssymParam(ib) =    1.04080910045869  
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 9) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =    2.04176676954650     
                InParticle%ScaEff(ib) =    1.74112647693926     
                InParticle%BackScaEff(ib) =   4.227385751022267E-002
                InParticle%AssymParam(ib) =    1.29350216586826     
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =    2.12983639968120     
                InParticle%ScaEff(ib) =    1.66122289119653     
                InParticle%BackScaEff(ib) =   6.941325831696683E-002
                InParticle%AssymParam(ib) =    1.23517403923357  
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =    2.65516028103601     
                InParticle%ScaEff(ib) =    1.88314938012702     
                InParticle%BackScaEff(ib) =   5.469645993818988E-002
                InParticle%AssymParam(ib) =    1.54336518258899  
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    2.13538718383939     
                InParticle%ScaEff(ib) =    1.05843352518967     
                InParticle%BackScaEff(ib) =   1.092278492080876E-004
                InParticle%AssymParam(ib) =   0.929098956036963   
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    2.42230784450767     
                InParticle%ScaEff(ib) =    1.22969039581765     
                InParticle%BackScaEff(ib) =   7.302126663631858E-003
                InParticle%AssymParam(ib) =    1.11814100390391  
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 10) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =    2.23747218460542     
                InParticle%ScaEff(ib) =    1.72439047618241     
                InParticle%BackScaEff(ib) =   6.409481896871368E-002
                InParticle%AssymParam(ib) =    1.47381594101870     
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =    2.06562574336866     
                InParticle%ScaEff(ib) =    1.45518814322546     
                InParticle%BackScaEff(ib) =   1.205943583104343E-002
                InParticle%AssymParam(ib) =    1.20942166076782  
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =    2.10833500508912     
                InParticle%ScaEff(ib) =    1.31014119463384     
                InParticle%BackScaEff(ib) =   7.973096319157684E-002
                InParticle%AssymParam(ib) =    1.10213226816747 
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    2.11640145406860     
                InParticle%ScaEff(ib) =    1.11670002015824     
                InParticle%BackScaEff(ib) =   1.780475953408376E-002
                InParticle%AssymParam(ib) =    1.03170971466352 
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    2.15041186500276     
                InParticle%ScaEff(ib) =    1.11393480479371     
                InParticle%BackScaEff(ib) =   7.529827215511730E-003
                InParticle%AssymParam(ib) =    1.04393547659690    
              END IF !ratiobin
            END IF !sizebin

          ELSE IF (ib .eq. 5) THEN !-------------------------------!MM-DEBUG
            IF (sizebin .EQ. 1) THEN
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =   2.156470752780501E-003
                InParticle%ScaEff(ib) =   2.423510895145302E-004
                InParticle%BackScaEff(ib) =   2.850558992752902E-005
                InParticle%AssymParam(ib) =   1.497584622399152E-006
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =   6.011100304721312E-003
                InParticle%ScaEff(ib) =   2.475458340153061E-004
                InParticle%BackScaEff(ib) =   2.912078589763700E-005
                InParticle%AssymParam(ib) =   1.514855723801651E-006
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =   2.039744966943594E-002
                InParticle%ScaEff(ib) =   2.692705667022480E-004
                InParticle%BackScaEff(ib) =   3.169077807365558E-005
                InParticle%AssymParam(ib) =   1.597228459102876E-006
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =   5.105526786922314E-002
                InParticle%ScaEff(ib) =   3.285750068981799E-004
                InParticle%BackScaEff(ib) =   3.868990015467943E-005
                InParticle%AssymParam(ib) =   1.883101191165021E-006
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =   9.935247288747233E-002
                InParticle%ScaEff(ib) =   4.615071333409716E-004
                InParticle%BackScaEff(ib) =   5.430775616233305E-005
                InParticle%AssymParam(ib) =   2.788408085006948E-006
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 2) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =   1.669458937281732E-002
                InParticle%ScaEff(ib) =   1.120209690850707E-002
                InParticle%BackScaEff(ib) =   1.209098321487180E-003
                InParticle%AssymParam(ib) =   4.628219868201358E-004
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =   2.837001252781295E-002
                InParticle%ScaEff(ib) =   1.146191869364921E-002
                InParticle%BackScaEff(ib) =   1.238543761667780E-003
                InParticle%AssymParam(ib) =   4.684342739501339E-004
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =   7.205891157290933E-002
                InParticle%ScaEff(ib) =   1.253454013758639E-002
                InParticle%BackScaEff(ib) =   1.359193147791398E-003
                InParticle%AssymParam(ib) =   4.950143410363731E-004
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =   0.165305049113183     
                InParticle%ScaEff(ib) =   1.541858216571215E-002
                InParticle%BackScaEff(ib) =   1.678107315956416E-003
                InParticle%AssymParam(ib) =   5.873438238030849E-004
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =   0.312501798369903     
                InParticle%ScaEff(ib) =   2.170880988992300E-002
                InParticle%BackScaEff(ib) =   2.350547480252358E-003
                InParticle%AssymParam(ib) =   8.776614745490395E-004
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 3) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =   8.048939535738715E-002
                InParticle%ScaEff(ib) =   7.043825104650552E-002
                InParticle%BackScaEff(ib) =   6.468355714037413E-003
                InParticle%AssymParam(ib) =   7.268310940692170E-003
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =   0.104420575985344     
                InParticle%ScaEff(ib) =   7.218307088083935E-002
                InParticle%BackScaEff(ib) =   6.653368078192959E-003
                InParticle%AssymParam(ib) =   7.353345246323822E-003
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =   0.193572106344506     
                InParticle%ScaEff(ib) =   7.919845083348445E-002
                InParticle%BackScaEff(ib) =   7.380076781413513E-003
                InParticle%AssymParam(ib) =   7.762654644333494E-003
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =   0.381053745280321     
                InParticle%ScaEff(ib) =   9.748035902239889E-002
                InParticle%BackScaEff(ib) =   9.168685947886111E-003
                InParticle%AssymParam(ib) =   9.244252821987917E-003
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =   0.670302709765434     
                InParticle%ScaEff(ib) =   0.135022061326670     
                InParticle%BackScaEff(ib) =   1.243710257185139E-002
                InParticle%AssymParam(ib) =   1.391330405960398E-002
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 4) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =   0.399820651387827     
                InParticle%ScaEff(ib) =   0.379883274074818     
                InParticle%BackScaEff(ib) =   2.020361226829106E-002
                InParticle%AssymParam(ib) =   0.103873792829933     
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =   0.455403541898221     
                InParticle%ScaEff(ib) =   0.386772613189943     
                InParticle%BackScaEff(ib) =   2.098751310955096E-002
                InParticle%AssymParam(ib) =   0.104049427405662   
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =   0.651763858720752     
                InParticle%ScaEff(ib) =   0.412126744622607     
                InParticle%BackScaEff(ib) =   2.344456978061815E-002
                InParticle%AssymParam(ib) =   0.106475504137779  
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    1.02870049467777     
                InParticle%ScaEff(ib) =   0.476681422317556     
                InParticle%BackScaEff(ib) =   2.713394500818145E-002
                InParticle%AssymParam(ib) =   0.123147638691622 
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    1.55707178531904     
                InParticle%ScaEff(ib) =   0.598791865162430     
                InParticle%BackScaEff(ib) =   2.878177527669097E-002
                InParticle%AssymParam(ib) =   0.177857561447019  
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 5) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =    1.62099816410536     
                InParticle%ScaEff(ib) =    1.57416424958781     
                InParticle%BackScaEff(ib) =   9.008312064238813E-003
                InParticle%AssymParam(ib) =   0.984514280597314     
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =    1.68427078358998     
                InParticle%ScaEff(ib) =    1.53976687599924     
                InParticle%BackScaEff(ib) =   1.108278493238152E-002
                InParticle%AssymParam(ib) =   0.955625625132237 
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =    1.86220978318272     
                InParticle%ScaEff(ib) =    1.38391969898805     
                InParticle%BackScaEff(ib) =   1.271187494500821E-002
                InParticle%AssymParam(ib) =   0.846451907574883    
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    2.12280361098686     
                InParticle%ScaEff(ib) =    1.18823762138111     
                InParticle%BackScaEff(ib) =   3.356462294672066E-003
                InParticle%AssymParam(ib) =   0.736746348599184 
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    2.50059318282653     
                InParticle%ScaEff(ib) =    1.20340374509707     
                InParticle%BackScaEff(ib) =   1.311716856288350E-003
                InParticle%AssymParam(ib) =   0.746498201158639    
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 6) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =    3.54078596830612     
                InParticle%ScaEff(ib) =    3.45847685303443     
                InParticle%BackScaEff(ib) =   3.198360178172965E-002
                InParticle%AssymParam(ib) =    2.51481322206534     
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =    3.38232496270729     
                InParticle%ScaEff(ib) =    3.13603488550063     
                InParticle%BackScaEff(ib) =   2.177095831839554E-002
                InParticle%AssymParam(ib) =    2.17920072508563 
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =    3.05185203148685     
                InParticle%ScaEff(ib) =    2.41454769078163     
                InParticle%BackScaEff(ib) =   1.359885498022188E-002
                InParticle%AssymParam(ib) =    1.61285574825469 
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    2.71020889915230     
                InParticle%ScaEff(ib) =    1.56121405249262     
                InParticle%BackScaEff(ib) =   4.120296335153626E-003
                InParticle%AssymParam(ib) =    1.15128600049189    
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    2.79775087376873     
                InParticle%ScaEff(ib) =    1.34534074511176     
                InParticle%BackScaEff(ib) =   1.266605735950606E-003
                InParticle%AssymParam(ib) =    1.05237794231973 
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 7) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =    3.79839753863319     
                InParticle%ScaEff(ib) =    3.67542762340395     
                InParticle%BackScaEff(ib) =   6.231530104948533E-002
                InParticle%AssymParam(ib) =    2.64329165893099     
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =    3.70555289052490     
                InParticle%ScaEff(ib) =    3.39630630241535     
                InParticle%BackScaEff(ib) =   5.870560457402267E-002
                InParticle%AssymParam(ib) =    2.39134165736987  
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =    3.32135560077935     
                InParticle%ScaEff(ib) =    2.61225061502338     
                InParticle%BackScaEff(ib) =   1.745859699296794E-002
                InParticle%AssymParam(ib) =    1.72559001421307  
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    2.76597986500317     
                InParticle%ScaEff(ib) =    1.52787613178405     
                InParticle%BackScaEff(ib) =   2.020718458469559E-003
                InParticle%AssymParam(ib) =    1.15470793231692 
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    2.52951896303713     
                InParticle%ScaEff(ib) =    1.10551449919514     
                InParticle%BackScaEff(ib) =   4.493611108713668E-003
                InParticle%AssymParam(ib) =   0.977491723508472  
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 8) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =    1.98269153957613     
                InParticle%ScaEff(ib) =    1.75070724654365     
                InParticle%BackScaEff(ib) =   0.596312156974750     
                InParticle%AssymParam(ib) =   0.814896613634202     
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =    2.13994561469170     
                InParticle%ScaEff(ib) =    1.71505990064196     
                InParticle%BackScaEff(ib) =   0.561875869574628     
                InParticle%AssymParam(ib) =   0.811718423690809  
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =    2.71924620632794     
                InParticle%ScaEff(ib) =    1.91427232338367     
                InParticle%BackScaEff(ib) =   0.450655248663453     
                InParticle%AssymParam(ib) =    1.11990749586247 
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    2.99844439004540     
                InParticle%ScaEff(ib) =    1.72913526143464     
                InParticle%BackScaEff(ib) =   7.159558927873695E-002
                InParticle%AssymParam(ib) =    1.39133971455852  
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    2.42621425011645     
                InParticle%ScaEff(ib) =    1.07919795900881     
                InParticle%BackScaEff(ib) =   1.832738097184638E-004
                InParticle%AssymParam(ib) =    1.02579094906008 
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 9) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =    2.51437206079227     
                InParticle%ScaEff(ib) =    2.23110602644525     
                InParticle%BackScaEff(ib) =   0.210648672235153     
                InParticle%AssymParam(ib) =    1.75170270019410     
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =    2.71636878080886     
                InParticle%ScaEff(ib) =    2.25775825646040     
                InParticle%BackScaEff(ib) =   0.194005985128595     
                InParticle%AssymParam(ib) =    1.79854133717791 
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =    2.45775870854036     
                InParticle%ScaEff(ib) =    1.66858194943098     
                InParticle%BackScaEff(ib) =   7.243187861157524E-002
                InParticle%AssymParam(ib) =    1.31182512900495  
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    2.25680078236096     
                InParticle%ScaEff(ib) =    1.10370970849717     
                InParticle%BackScaEff(ib) =   1.730107992007260E-002
                InParticle%AssymParam(ib) =   0.978909629046713  
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    2.41970148126958     
                InParticle%ScaEff(ib) =    1.16623188123296     
                InParticle%BackScaEff(ib) =   4.823501579071294E-003
                InParticle%AssymParam(ib) =    1.08759811723919 
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 10) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =    2.37604420245805     
                InParticle%ScaEff(ib) =    1.90518173079963     
                InParticle%BackScaEff(ib) =   6.002949337120361E-002
                InParticle%AssymParam(ib) =    1.59334484865810     
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =    2.13488567918209     
                InParticle%ScaEff(ib) =    1.55166581678293     
                InParticle%BackScaEff(ib) =   1.588575944612918E-002
                InParticle%AssymParam(ib) =    1.24950632183058 
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =    2.40555145061762     
                InParticle%ScaEff(ib) =    1.60765353248846     
                InParticle%BackScaEff(ib) =   9.845384462499929E-002
                InParticle%AssymParam(ib) =    1.35685469246507  
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    2.29395827727994     
                InParticle%ScaEff(ib) =    1.28530950905199     
                InParticle%BackScaEff(ib) =   5.914862772105029E-003
                InParticle%AssymParam(ib) =    1.18513281110847 
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    2.19533391582167     
                InParticle%ScaEff(ib) =    1.14522369107149     
                InParticle%BackScaEff(ib) =   1.248136068711735E-003
                InParticle%AssymParam(ib) =    1.11011499130817  
              END IF !ratiobin
            END IF !sizebin

          ELSE IF (ib .eq. 6) THEN !-------------------------------!MM-DEBUG
            IF (sizebin .EQ. 1) THEN
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =   1.745990730877233E-003
                InParticle%ScaEff(ib) =   7.391174987435961E-005
                InParticle%BackScaEff(ib) =   8.751096344873538E-006
                InParticle%AssymParam(ib) =   2.526572623971385E-007
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =   4.440939533503097E-003
                InParticle%ScaEff(ib) =   7.546527233444613E-005
                InParticle%BackScaEff(ib) =   8.935722060557604E-006
                InParticle%AssymParam(ib) =   2.555310567509572E-007
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =   1.449796849247914E-002
                InParticle%ScaEff(ib) =   8.192307596275294E-005
                InParticle%BackScaEff(ib) =   9.702702563658139E-006
                InParticle%AssymParam(ib) =   2.692444256782151E-007
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =   3.593629234104791E-002
                InParticle%ScaEff(ib) =   9.935145909251058E-005
                InParticle%BackScaEff(ib) =   1.176989932640274E-005
                InParticle%AssymParam(ib) =   3.163832137897422E-007
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =   6.974674587368791E-002
                InParticle%ScaEff(ib) =   1.380002289731614E-004
                InParticle%BackScaEff(ib) =   1.634264801108306E-005
                InParticle%AssymParam(ib) =   4.634947479846343E-007
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 2) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =   8.003059671339458E-003
                InParticle%ScaEff(ib) =   3.403547228659068E-003
                InParticle%BackScaEff(ib) =   3.844268516783187E-004
                InParticle%AssymParam(ib) =   7.809869126870609E-005
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =   1.562305355258488E-002
                InParticle%ScaEff(ib) =   3.478762198184714E-003
                InParticle%BackScaEff(ib) =   3.931456551665856E-004
                InParticle%AssymParam(ib) =   7.902165898845314E-005
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =   4.411449864602950E-002
                InParticle%ScaEff(ib) =   3.789242520650595E-003
                InParticle%BackScaEff(ib) =   4.289828086256152E-004
                InParticle%AssymParam(ib) =   8.339861704319302E-005
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =   0.104963438764179     
                InParticle%ScaEff(ib) =   4.620191164186989E-003
                InParticle%BackScaEff(ib) =   5.240207479496873E-004
                InParticle%AssymParam(ib) =   9.839738585336499E-005
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =   0.201198483519751     
                InParticle%ScaEff(ib) =   6.436651325675930E-003
                InParticle%BackScaEff(ib) =   7.281489100011894E-004
                InParticle%AssymParam(ib) =   1.449333341115026E-004
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 3) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =   2.958896551755912E-002
                InParticle%ScaEff(ib) =   2.161440615214868E-002
                InParticle%BackScaEff(ib) =   2.241389981346530E-003
                InParticle%AssymParam(ib) =   1.235340318836642E-003
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =   4.362928843777036E-002
                InParticle%ScaEff(ib) =   2.212396060954086E-002
                InParticle%BackScaEff(ib) =   2.298035921877446E-003
                InParticle%AssymParam(ib) =   1.250376698653324E-003
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =   9.617834969558516E-002
                InParticle%ScaEff(ib) =   2.420001626704581E-002
                InParticle%BackScaEff(ib) =   2.526221633178067E-003
                InParticle%AssymParam(ib) =   1.321561677527488E-003
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =   0.208229242962370     
                InParticle%ScaEff(ib) =   2.966670186588255E-002
                InParticle%BackScaEff(ib) =   3.111995027755152E-003
                InParticle%AssymParam(ib) =   1.567100051923124E-003
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =   0.384882292364829     
                InParticle%ScaEff(ib) =   4.126765149541759E-002
                InParticle%BackScaEff(ib) =   4.294187520066071E-003
                InParticle%AssymParam(ib) =   2.325522236346111E-003
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 4) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =   0.147289932438826     
                InParticle%ScaEff(ib) =   0.132262363149071     
                InParticle%BackScaEff(ib) =   1.083488882099055E-002
                InParticle%AssymParam(ib) =   1.900024659026915E-002
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =   0.178214429058716     
                InParticle%ScaEff(ib) =   0.135513949560788     
                InParticle%BackScaEff(ib) =   1.116901832652692E-002
                InParticle%AssymParam(ib) =   1.920097283927704E-002
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =   0.292463937989071     
                InParticle%ScaEff(ib) =   0.148240061353616     
                InParticle%BackScaEff(ib) =   1.242417396618586E-002
                InParticle%AssymParam(ib) =   2.019830645151947E-002
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =   0.528374841368533     
                InParticle%ScaEff(ib) =   0.180210496949550     
                InParticle%BackScaEff(ib) =   1.526841985223957E-002
                InParticle%AssymParam(ib) =   2.394712939964825E-002
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =   0.883780761355177     
                InParticle%ScaEff(ib) =   0.242377608224781     
                InParticle%BackScaEff(ib) =   1.973428625847923E-002
                InParticle%AssymParam(ib) =   3.563398716507265E-002
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 5) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =   0.655866354883852     
                InParticle%ScaEff(ib) =   0.625022675951519     
                InParticle%BackScaEff(ib) =   1.773540365132158E-002
                InParticle%AssymParam(ib) =   0.253551650833792     
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =   0.721449871409231     
                InParticle%ScaEff(ib) =   0.630362557417295     
                InParticle%BackScaEff(ib) =   1.888538677092272E-002
                InParticle%AssymParam(ib) =   0.251659400008581 
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =   0.947940060498512     
                InParticle%ScaEff(ib) =   0.649338032078007     
                InParticle%BackScaEff(ib) =   2.157892979403054E-002
                InParticle%AssymParam(ib) =   0.250707767889194 
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    1.37699128271348     
                InParticle%ScaEff(ib) =   0.710893696073668     
                InParticle%BackScaEff(ib) =   2.209180435872176E-002
                InParticle%AssymParam(ib) =   0.280300293979218 
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    1.95333915266329     
                InParticle%ScaEff(ib) =   0.843706535935239     
                InParticle%BackScaEff(ib) =   1.734450453694417E-002
                InParticle%AssymParam(ib) =   0.370309187029721 
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 6) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =    2.21155597637799     
                InParticle%ScaEff(ib) =    2.15005296254870     
                InParticle%BackScaEff(ib) =   3.842774817988816E-002
                InParticle%AssymParam(ib) =    1.32837581497928     
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =    2.24185436275420     
                InParticle%ScaEff(ib) =    2.06553593504804     
                InParticle%BackScaEff(ib) =   3.562590509315052E-002
                InParticle%AssymParam(ib) =    1.27824338108958  
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =    2.25389419640571     
                InParticle%ScaEff(ib) =    1.73367031909003     
                InParticle%BackScaEff(ib) =   2.875169477837651E-002
                InParticle%AssymParam(ib) =    1.08150467565139 
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    2.38277647808416     
                InParticle%ScaEff(ib) =    1.39719657679278     
                InParticle%BackScaEff(ib) =   1.276567136944368E-002
                InParticle%AssymParam(ib) =   0.927166406690145  
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    2.75039792073941     
                InParticle%ScaEff(ib) =    1.36715119187155     
                InParticle%BackScaEff(ib) =   9.849083369577715E-003
                InParticle%AssymParam(ib) =   0.947433996589785   
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 7) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =    4.12200574278341     
                InParticle%ScaEff(ib) =    4.00174704435463     
                InParticle%BackScaEff(ib) =   7.175774120730871E-002
                InParticle%AssymParam(ib) =    2.94243861164802     
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =    3.93526119047506     
                InParticle%ScaEff(ib) =    3.65772322316503     
                InParticle%BackScaEff(ib) =   6.286334426448967E-002
                InParticle%AssymParam(ib) =    2.61403874572427 
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =    3.47453652819076     
                InParticle%ScaEff(ib) =    2.78560024289033     
                InParticle%BackScaEff(ib) =   2.579655482451309E-002
                InParticle%AssymParam(ib) =    1.90668298076763  
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    2.90209025062548     
                InParticle%ScaEff(ib) =    1.65039032725206     
                InParticle%BackScaEff(ib) =   4.996210538022457E-003
                InParticle%AssymParam(ib) =    1.25721943136780   
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    2.71994191970323     
                InParticle%ScaEff(ib) =    1.27893743131628     
                InParticle%BackScaEff(ib) =   6.444624170997618E-003
                InParticle%AssymParam(ib) =    1.04453110948951  
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 8) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =    2.95992740299252     
                InParticle%ScaEff(ib) =    2.79135476125437     
                InParticle%BackScaEff(ib) =   0.140381758860169     
                InParticle%AssymParam(ib) =    1.79503415913209     
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =    3.07041978015654     
                InParticle%ScaEff(ib) =    2.72580740171228     
                InParticle%BackScaEff(ib) =   0.111857495874945     
                InParticle%AssymParam(ib) =    1.75174584094881  
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =    3.11849224338763     
                InParticle%ScaEff(ib) =    2.36906649499724     
                InParticle%BackScaEff(ib) =   6.687453488181348E-002
                InParticle%AssymParam(ib) =    1.54733378731707  
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    2.80619569449233     
                InParticle%ScaEff(ib) =    1.57642238421905     
                InParticle%BackScaEff(ib) =   1.038199592688869E-002
                InParticle%AssymParam(ib) =    1.19575528096164   
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    2.47860056315199     
                InParticle%ScaEff(ib) =    1.08036129578137     
                InParticle%BackScaEff(ib) =   1.181422890651614E-003
                InParticle%AssymParam(ib) =   0.987927800861476    
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 9) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =    2.26311639619112     
                InParticle%ScaEff(ib) =    1.98787413805281     
                InParticle%BackScaEff(ib) =   0.262094922520487     
                InParticle%AssymParam(ib) =    1.32193503285332     
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =    2.06766010618677     
                InParticle%ScaEff(ib) =    1.61556765233881     
                InParticle%BackScaEff(ib) =   0.241008695412522     
                InParticle%AssymParam(ib) =   0.962079949615984   
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =    2.07358117616456     
                InParticle%ScaEff(ib) =    1.27957179697689     
                InParticle%BackScaEff(ib) =   0.150505123505038     
                InParticle%AssymParam(ib) =   0.719201438179157  
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    2.71080656248684     
                InParticle%ScaEff(ib) =    1.49868666323326     
                InParticle%BackScaEff(ib) =   4.339057737712454E-002
                InParticle%AssymParam(ib) =    1.24784272622554 
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    2.41482825444001     
                InParticle%ScaEff(ib) =    1.09375391218603     
                InParticle%BackScaEff(ib) =   1.237477201533048E-003
                InParticle%AssymParam(ib) =    1.04390195376583     
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 10) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =    2.27037858457026     
                InParticle%ScaEff(ib) =    1.81605110214765     
                InParticle%BackScaEff(ib) =   5.576530762260055E-002
                InParticle%AssymParam(ib) =    1.48495089386339     
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =    2.05247974159767     
                InParticle%ScaEff(ib) =    1.47654312653192     
                InParticle%BackScaEff(ib) =   0.125019234074214     
                InParticle%AssymParam(ib) =    1.15154802069341     
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =    2.46409375074626     
                InParticle%ScaEff(ib) =    1.65435627485793     
                InParticle%BackScaEff(ib) =   0.170163724610554     
                InParticle%AssymParam(ib) =    1.38854140585308   
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    2.29582319104450     
                InParticle%ScaEff(ib) =    1.24745592736670     
                InParticle%BackScaEff(ib) =   1.843349596407969E-002
                InParticle%AssymParam(ib) =    1.13258828871739  
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    2.29684611193584     
                InParticle%ScaEff(ib) =    1.26036871862716     
                InParticle%BackScaEff(ib) =   5.638748954638020E-003
                InParticle%AssymParam(ib) =    1.13886041054344     
              END IF !ratiobin
            END IF !sizebin

          ELSE IF (ib .eq. 7) THEN !-------------------------------!MM-DEBUG  
            IF (sizebin .EQ. 1) THEN
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =   1.822525031820004E-003
                InParticle%ScaEff(ib) =   3.301443988157131E-005
                InParticle%BackScaEff(ib) =   3.919444865237387E-006
                InParticle%AssymParam(ib) =   7.547501792487873E-008
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =   3.976783558301614E-003
                InParticle%ScaEff(ib) =   3.370752387684337E-005
                InParticle%BackScaEff(ib) =   4.001932571094586E-006
                InParticle%AssymParam(ib) =   7.633503523585032E-008
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =   1.201545347492798E-002
                InParticle%ScaEff(ib) =   3.658147190747974E-005
                InParticle%BackScaEff(ib) =   4.343829425194149E-006
                InParticle%AssymParam(ib) =   8.043689139507312E-008
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =   2.915180842891664E-002
                InParticle%ScaEff(ib) =   4.430108466274274E-005
                InParticle%BackScaEff(ib) =   5.261370472705972E-006
                InParticle%AssymParam(ib) =   9.447242553291143E-008
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =   5.618368287014076E-002
                InParticle%ScaEff(ib) =   6.134127290403799E-005
                InParticle%BackScaEff(ib) =   7.283375779648253E-006
                InParticle%AssymParam(ib) =   1.379899582095978E-007
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 2) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =   6.350211255988941E-003
                InParticle%ScaEff(ib) =   1.516929500974525E-003
                InParticle%BackScaEff(ib) =   1.745174041177409E-004
                InParticle%AssymParam(ib) =   2.332714440848649E-005
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =   1.227184437690569E-002
                InParticle%ScaEff(ib) =   1.549921682674426E-003
                InParticle%BackScaEff(ib) =   1.783785489322480E-004
                InParticle%AssymParam(ib) =   2.360067697563908E-005
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =   3.439992174863217E-002
                InParticle%ScaEff(ib) =   1.686098531844482E-003
                InParticle%BackScaEff(ib) =   1.942696076786347E-004
                InParticle%AssymParam(ib) =   2.489889273816199E-005
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =   8.164294863008409E-002
                InParticle%ScaEff(ib) =   2.049861918610712E-003
                InParticle%BackScaEff(ib) =   2.364616196182227E-004
                InParticle%AssymParam(ib) =   2.932714034242508E-005
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =   0.156339333036700     
                InParticle%ScaEff(ib) =   2.845401925807674E-003
                InParticle%BackScaEff(ib) =   3.276771175075854E-004
                InParticle%AssymParam(ib) =   4.299545171349382E-005
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 3) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =   1.778645879661125E-002
                InParticle%ScaEff(ib) =   9.630624272648958E-003
                InParticle%BackScaEff(ib) =   1.047199790254237E-003
                InParticle%AssymParam(ib) =   3.693591132836390E-004
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =   2.817176459979360E-002
                InParticle%ScaEff(ib) =   9.850994744658644E-003
                InParticle%BackScaEff(ib) =   1.072254685596641E-003
                InParticle%AssymParam(ib) =   3.738316430915966E-004
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =   6.703824339568021E-002
                InParticle%ScaEff(ib) =   1.075298035076666E-002
                InParticle%BackScaEff(ib) =   1.174049194434285E-003
                InParticle%AssymParam(ib) =   3.949759463885883E-004
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =   0.150070514170024     
                InParticle%ScaEff(ib) =   1.313740414208335E-002
                InParticle%BackScaEff(ib) =   1.438868893717062E-003
                InParticle%AssymParam(ib) =   4.671294722475529E-004
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =   0.281435205955893     
                InParticle%ScaEff(ib) =   1.825633648332662E-002
                InParticle%BackScaEff(ib) =   1.990015549345189E-003
                InParticle%AssymParam(ib) =   6.887746698199976E-004
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 4) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =   7.511537565154930E-002
                InParticle%ScaEff(ib) =   6.037039703762549E-002
                InParticle%BackScaEff(ib) =   5.661249225289418E-003
                InParticle%AssymParam(ib) =   5.763441149550745E-003
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =   9.608087097199536E-002
                InParticle%ScaEff(ib) =   6.185787194192215E-002
                InParticle%BackScaEff(ib) =   5.819799652834955E-003
                InParticle%AssymParam(ib) =   5.832797077916601E-003
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =   0.174315781725701     
                InParticle%ScaEff(ib) =   6.780390801172184E-002
                InParticle%BackScaEff(ib) =   6.439851071284188E-003
                InParticle%AssymParam(ib) =   6.163905103196722E-003
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =   0.339577525959736     
                InParticle%ScaEff(ib) =   8.305973787255339E-002
                InParticle%BackScaEff(ib) =   7.951594246883955E-003
                InParticle%AssymParam(ib) =   7.325304521773114E-003
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =   0.596353723504103     
                InParticle%ScaEff(ib) =   0.114051583657035     
                InParticle%BackScaEff(ib) =   1.072263296429436E-002
                InParticle%AssymParam(ib) =   1.088757311624777E-002
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 5) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =   0.367213460887450     
                InParticle%ScaEff(ib) =   0.338129863781664     
                InParticle%BackScaEff(ib) =   1.942458711762533E-002
                InParticle%AssymParam(ib) =   8.548555144248844E-002
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =   0.416908850227615     
                InParticle%ScaEff(ib) =   0.345048113260993     
                InParticle%BackScaEff(ib) =   2.015557361575887E-002
                InParticle%AssymParam(ib) =   8.586748845454591E-002
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =   0.594332181245844     
                InParticle%ScaEff(ib) =   0.370529389685645     
                InParticle%BackScaEff(ib) =   2.251890121891220E-002
                InParticle%AssymParam(ib) =   8.865139920209882E-002
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =   0.939935921136212     
                InParticle%ScaEff(ib) =   0.432391959185205     
                InParticle%BackScaEff(ib) =   2.636799177643430E-002
                InParticle%AssymParam(ib) =   0.103202829981282  
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    1.43145499993234     
                InParticle%ScaEff(ib) =   0.545537007405742     
                InParticle%BackScaEff(ib) =   2.893519960973154E-002
                InParticle%AssymParam(ib) =   0.149244833889450    
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 6) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =    1.45816059241186     
                InParticle%ScaEff(ib) =    1.38995743041681     
                InParticle%BackScaEff(ib) =   4.843928180701920E-003
                InParticle%AssymParam(ib) =   0.863975288788868     
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =    1.52390102005966     
                InParticle%ScaEff(ib) =    1.36972244072900     
                InParticle%BackScaEff(ib) =   6.940145129539475E-003
                InParticle%AssymParam(ib) =   0.843783357047611  
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =    1.73431959302154     
                InParticle%ScaEff(ib) =    1.27869184206174     
                InParticle%BackScaEff(ib) =   9.603752900582917E-003
                InParticle%AssymParam(ib) =   0.773201646713378   
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    2.06314074085207     
                InParticle%ScaEff(ib) =    1.16139320271965     
                InParticle%BackScaEff(ib) =   2.764357017880969E-003
                InParticle%AssymParam(ib) =   0.706554457231506    
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    2.45645089565132     
                InParticle%ScaEff(ib) =    1.19339683407815     
                InParticle%BackScaEff(ib) =   6.017506001740586E-004
                InParticle%AssymParam(ib) =   0.722065171180426     
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 7) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =    3.45951821027362     
                InParticle%ScaEff(ib) =    3.33678681914849     
                InParticle%BackScaEff(ib) =   1.893751813256706E-002
                InParticle%AssymParam(ib) =    2.44167514622609     
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =    3.34310639991427     
                InParticle%ScaEff(ib) =    3.06831050588359     
                InParticle%BackScaEff(ib) =   1.005370487819418E-002
                InParticle%AssymParam(ib) =    2.16582495070475  
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =    3.07849309828382     
                InParticle%ScaEff(ib) =    2.43162727982638     
                InParticle%BackScaEff(ib) =   6.015725161101985E-003
                InParticle%AssymParam(ib) =    1.65243326082859  
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    2.75540223704259     
                InParticle%ScaEff(ib) =    1.61943232037147     
                InParticle%BackScaEff(ib) =   1.834465461441918E-003
                InParticle%AssymParam(ib) =    1.18930306573657  
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    2.82286734244914     
                InParticle%ScaEff(ib) =    1.39140112703455     
                InParticle%BackScaEff(ib) =   2.767346080934591E-003
                InParticle%AssymParam(ib) =    1.07134563671328   
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 8) THEN !!!!!!!!!!!!!!!!!!!!!              
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =    4.05511210350117     
                InParticle%ScaEff(ib) =    3.86942445930530     
                InParticle%BackScaEff(ib) =   8.359854479853722E-002
                InParticle%AssymParam(ib) =    2.83651622731413     
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =    3.94244140337470     
                InParticle%ScaEff(ib) =    3.58637491300409     
                InParticle%BackScaEff(ib) =   7.429315165712792E-002
                InParticle%AssymParam(ib) =    2.58680033029294 
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =    3.51060204855541     
                InParticle%ScaEff(ib) =    2.77177769656355     
                InParticle%BackScaEff(ib) =   3.406103442355219E-002
                InParticle%AssymParam(ib) =    1.89292821121038     
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    2.80140587662272     
                InParticle%ScaEff(ib) =    1.54935066760967     
                InParticle%BackScaEff(ib) =   5.461204278651640E-003
                InParticle%AssymParam(ib) =    1.17934121311270     
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    2.54983639259369     
                InParticle%ScaEff(ib) =    1.12259967128459     
                InParticle%BackScaEff(ib) =   3.408981499861641E-003
                InParticle%AssymParam(ib) =   0.979790120132008   
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 9) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =    1.78841120186596     
                InParticle%ScaEff(ib) =    1.50285152352527     
                InParticle%BackScaEff(ib) =   0.256009325703574     
                InParticle%AssymParam(ib) =   0.724721746210302     
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =    2.00646138691121     
                InParticle%ScaEff(ib) =    1.55006876500885     
                InParticle%BackScaEff(ib) =   0.240575656297192     
                InParticle%AssymParam(ib) =   0.770042000363660 
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =    2.54595095593060     
                InParticle%ScaEff(ib) =    1.73068609919718     
                InParticle%BackScaEff(ib) =   0.198090514450071     
                InParticle%AssymParam(ib) =    1.07049084651018 
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    2.83394281860531     
                InParticle%ScaEff(ib) =    1.61873736873039     
                InParticle%BackScaEff(ib) =   3.734013425565764E-002
                InParticle%AssymParam(ib) =    1.28031725207590 
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    2.43717515107114     
                InParticle%ScaEff(ib) =    1.08346654076758     
                InParticle%BackScaEff(ib) =   5.314656322385932E-004
                InParticle%AssymParam(ib) =    1.02602074598043   
              END IF !ratiobin
            ELSE IF (sizebin .EQ. 10) THEN !!!!!!!!!!!!!!!!!!!!!
              IF (ratiobin .EQ. 1) THEN
                InParticle%ExtEff(ib) =    2.45011338214344     
                InParticle%ScaEff(ib) =    1.95615262558493     
                InParticle%BackScaEff(ib) =   3.891223075565580E-002
                InParticle%AssymParam(ib) =    1.59939995051104     
              ELSE IF (ratiobin .EQ. 2) THEN
                InParticle%ExtEff(ib) =    2.49067116114546     
                InParticle%ScaEff(ib) =    1.87809892686616     
                InParticle%BackScaEff(ib) =   0.105901825728398     
                InParticle%AssymParam(ib) =    1.52462700764817  
              ELSE IF (ratiobin .EQ. 3) THEN
                InParticle%ExtEff(ib) =    2.04981903293732     
                InParticle%ScaEff(ib) =    1.20916697342055     
                InParticle%BackScaEff(ib) =   0.108478335713496     
                InParticle%AssymParam(ib) =   0.912513517238510  
              ELSE IF (ratiobin .EQ. 4) THEN
                InParticle%ExtEff(ib) =    2.40722168261355     
                InParticle%ScaEff(ib) =    1.30811672639663     
                InParticle%BackScaEff(ib) =   5.825940715847113E-003
                InParticle%AssymParam(ib) =    1.22516460991257    
              ELSE IF (ratiobin .EQ. 5) THEN
                InParticle%ExtEff(ib) =    2.37166743841577     
                InParticle%ScaEff(ib) =    1.27322805612039     
                InParticle%BackScaEff(ib) =   9.146157931267608E-003
                InParticle%AssymParam(ib) =    1.13626397409937  
              END IF !ratiobin
            END IF !sizebin

          ELSE                                          !MM-DEBUG
            WRITE(*,*) "ib = ",ib
            call error("Unknown value for Wavelength ib. STOP")
          ENDIF !ib

	   !IF (ib .eq. 419 .or. ib .eq. 420) THEN
                !WRITE(*,*) "ExtEff: ", InParticle%ExtEff(ib)
                !WRITE(*,*) "ScaEff: ", InParticle%ScaEff(ib)
                !WRITE(*,*) "SSA: ", InParticle%ScaEff(ib)/InParticle%ExtEff(ib)
                !WRITE(*,*) "AssymParam: ", InParticle%AssymParam(ib)
           !ENDIF

	   !DMiLay gives the assymetry parameter times the scattering efficiency
	   !So divide by scattering efficiency here (MJA, 10/15/07)
	   InParticle%AssymParam(ib) = InParticle%AssymParam(ib)/InParticle%ScaEff(ib)

          !Now call for a hypothetical particle of
          !same effective size but all BC
!          WRITE(*,*) "Call 2"
!	  CALL DMiLay( InParticle%AbsCoreRad, &
!			InParticle%EffectiveRadius, &
!			Wavelength, &
!			SolarCoreInd, & 
!			SolarCoreInd, &
!			CosAng, &
!			NumAng, &
!			InParticle%ExtEffBC(ib), &
!			InParticle%ScaEffBC(ib), &
!			InParticle%BackScaEffBC(ib), &
!			InParticle%AssymParamBC(ib), & 
!                        M1, M2, S21, D21, MAXANG )

	   !DMiLay gives the assymetry parameter times the scattering eff.
	   !So divide by scattering efficiency here (MJA, 10/15/07)
!	   InParticle%AssymParamBC(ib) = &
!                InParticle%AssymParamBC(ib)/InParticle%ScaEffBC(ib)

          !Now call for a hypothetical particle of
          !same effective size but all shell
!          WRITE(*,*) "Call 3"
!	  CALL DMiLay( InParticle%AbsCoreRad, &
!			InParticle%EffectiveRadius, &
!			Wavelength, &
!			SolarShellInd, & 
!			SolarShellInd, &
!			CosAng, &
!			NumAng, &
!			InParticle%ExtEffShell(ib), &
!			InParticle%ScaEffShell(ib), &
!			InParticle%BackScaEffShell(ib), &
!			InParticle%AssymParamShell(ib), & 
!                        M1, M2, S21, D21, MAXANG )

	   !DMiLay gives the assymetry parameter times the scattering eff.
	   !So divide by scattering efficiency here (MJA, 10/15/07)
!	   InParticle%AssymParamShell(ib) = &
!                InParticle%AssymParamShell(ib)/InParticle%ScaEffShell(ib)

!          WRITE(*,*) "Call 4"
          !Now call for volume average dielectric constant mixing rule
!	  CALL DMiLay( InParticle%AbsCoreRad, &
!			InParticle%EffectiveRadius, &
!			Wavelength, &
!			VAInd, & 
!			VAInd, &
!			CosAng, &
!			NumAng, &
!			InParticle%ExtEffVA(ib), &
!			InParticle%ScaEffVA(ib), &
!			InParticle%BackScaEffVA(ib), &
!			InParticle%AssymParamVA(ib), & 
!                        M1, M2, S21, D21, MAXANG )

	   !DMiLay gives the assymetry parameter times the scattering eff.
	   !So divide by scattering efficiency here (MJA, 10/15/07)
!	   InParticle%AssymParamVA(ib) = &
!                InParticle%AssymParamVA(ib)/InParticle%ScaEffVA(ib)

          !Now call for Maxwell-Garnett mixing rule
!	  CALL DMiLay( InParticle%AbsCoreRad, &
!			InParticle%EffectiveRadius, &
!			Wavelength, &
!			MGInd, & 
!			MGInd, &
!			CosAng, &
!			NumAng, &
!			InParticle%ExtEffMG(ib), &
!			InParticle%ScaEffMG(ib), &
!			InParticle%BackScaEffMG(ib), &
!			InParticle%AssymParamMG(ib), & 
!                        M1, M2, S21, D21, MAXANG )

	   !DMiLay gives the assymetry parameter times the scattering eff.
	   !So divide by scattering efficiency here (MJA, 10/15/07)
!	   InParticle%AssymParamMG(ib) = &
!                InParticle%AssymParamMG(ib)/InParticle%ScaEffMG(ib)
          
	END DO

        DEALLOCATE(AqVolConc,OrgVolConc,AqOrgVolConc, MolFrac, wvln_array, STAT = I)

	RETURN
END SUBROUTINE ShellRefIndAndRad

!Calculates aerosol optical depth for a grid point at 550 nm.
!Can only be used in inversion mode!

! Modified by CMB (AER): Make Gfortran cooperate with type safety
!                        as described above                         
REAL*8 FUNCTION AerosolOpticalDepth ()
! 08/30/2010 Changed for photolysis radiation bins MJA
! 02/27/2012 Changed for photolysis radiation bins MJA	
	USE ModelParameters, ONLY : PI, InversionHeight

	IMPLICIT NONE
	
	Type(Particle),POINTER :: p_Particle

	AerosolOpticalDepth = 0.
	
	p_Particle => Particles%First

	DO WHILE (ASSOCIATED(p_Particle))

           AerosolOpticalDepth = AerosolOpticalDepth + &
                p_Particle%NumberofParticles*&
                    (PI*p_Particle%EffectiveRadius**2)*p_Particle%ExtEff(301)

           p_Particle => p_Particle%Next

	END DO
	
	AerosolOpticalDepth = AerosolOpticalDepth*InversionHeight
	
	RETURN
END FUNCTION AerosolOpticalDepth

!Calculates total aerosol extinction coefficient (m**-1) 
!average single scattering albedo, and 
!average assymetry factor 
!for a grid point at 18 radiative bands.
SUBROUTINE AerosolOptProp (ExtCoeff, SingScat, Assym, BackScatCoeff, &
                          SubExtCoeff, SubSingScat, Flag)!&
                          !ExtCoeffSep, SingScatSep, AssymSep, BackScatCoeffSep, &
                          !SubExtCoeffSep, SubSingScatSep)
! 08/30/2010 Changed for the photolysis wavebands MJA
! 02/27/2012 Changed for photolysis radiation bins MJA	
! 05/03/2012 Changed to calculate backscattering coefficient, units
!             for extcoeff and backscatcoeff now Mm^-1
! 08/16/2012 Fixed backscattering coefficient, added submicron extinction & SSA
! 08/17/2012 Added capability for external mixtures
! 11/08/2012 Added flag to decide between mixing rules
!            Flag = 0, core-in-shell (default)
!            Flag = 1, external mixture of BC and rest
!            Flag = 2, volume average dielectric constant mixing rule
!            Flag = 3, Maxwell-Garnett mixing rule
!
	USE ModelParameters, ONLY : PI, InversionHeight

	IMPLICIT NONE

        !Input Variables
        INTEGER :: Flag
	!Output Variables
	REAL*8, DIMENSION(451) :: ExtCoeff, SingScat, Assym, BackScatCoeff
	REAL*8, DIMENSION(451) :: SubExtCoeff, SubSingScat
        !Below are for external mixture of BC and "shell"
!	REAL*8, DIMENSION(451) :: ExtCoeffSep, SingScatSep, AssymSep, BackScatCoeffSep
!	REAL*8, DIMENSION(451) :: SubExtCoeffSep, SubSingScatSep

	!Internal variables
	INTEGER :: ib, numbin
        LOGICAL :: FASTTUV
	REAL*8, DIMENSION(451) :: ScatCoeff, SubScatCoeff
!	REAL*8, DIMENSION(451) :: ScatCoeffSep, SubScatCoeffSep
	REAL*8 :: AssymSum,xxx, hbs_frac, assym_sq, bc_vol_frac

        ! CMB: Gfortran doesn't like the additional definitions of the 
        !      public variables, nor does it like "Particle" being used
        !      both as a type and a variable name
	Type(Particle),POINTER :: p_Particle


        DO ib = 1, 451
		ExtCoeff(ib) = 0.0
		ScatCoeff(ib) = 0.0
		BackScatCoeff(ib) = 0.0
		SubExtCoeff(ib) = 0.0
		SubScatCoeff(ib) = 0.0
                AssymSum = 0.0
!		ExtCoeffSep(ib) = 0.0
!		ScatCoeffSep(ib) = 0.0
!		BackScatCoeffSep(ib) = 0.0
!		SubExtCoeffSep(ib) = 0.0
!		SubScatCoeffSep(ib) = 0.0
!                AssymSumSep = 0.0
		p_Particle => Particles%First

		DO WHILE (ASSOCIATED(p_Particle))
                     IF (p_Particle%EffectiveRadius .gt. 0) THEN 
                       bc_vol_frac = (p_Particle%AbsCoreRad/&
                               p_Particle%EffectiveRadius)**3
                     ELSE
                       bc_vol_frac = 0.0
                     ENDIF
		     xxx = p_Particle%NumberofParticles*&
                        (PI*p_Particle%EffectiveRadius**2)

                     IF (FLAG .EQ. 0) THEN !Core in Shell  
			ExtCoeff(ib) = ExtCoeff(ib) + &
                            xxx*p_Particle%ExtEff(ib)
			ScatCoeff(ib) = ScatCoeff(ib) + &
                            xxx*p_Particle%ScaEff(ib)
			AssymSum = AssymSum + &
                            xxx*p_Particle%ScaEff(ib)*p_Particle%AssymParam(ib)

                        !Particles with diameters less than 1 micron
                        IF (2*p_Particle%EffectiveRadius .LT. 1.0e-4) THEN
 			  SubExtCoeff(ib) = SubExtCoeff(ib) + &
                            xxx*p_Particle%ExtEff(ib)
			  SubScatCoeff(ib) = SubScatCoeff(ib) + &
                            xxx*p_Particle%ScaEff(ib)
                        ENDIF

                     ELSE IF (Flag .EQ. 1) THEN !External mixture
                        ExtCoeff(ib) = ExtCoeff(ib) + &
                            bc_vol_frac*xxx*p_Particle%ExtEffBC(ib)&
                            +(1.0-bc_vol_frac)*xxx*p_Particle%ExtEffShell(ib)
			ScatCoeff(ib) = ScatCoeff(ib) + &
                            bc_vol_frac*xxx*p_Particle%ScaEffBC(ib) &
                            +(1.0-bc_vol_frac)*xxx*p_Particle%ScaEffShell(ib)
			AssymSum = AssymSum + &
                            bc_vol_frac*xxx*p_Particle%ScaEffBC(ib)*&
                                p_Particle%AssymParamBC(ib) &
                                       +(1.0-bc_vol_frac)*xxx* &
                                       p_Particle%ScaEffShell(ib)*&
                                       p_Particle%AssymParamShell(ib) 
                        !Particles with diameters less than 1 micron
                        IF (2*p_Particle%EffectiveRadius .LT. 1.0e-4) THEN	
                          SubExtCoeff(ib) = SubExtCoeff(ib) + &
                                        bc_vol_frac*xxx*p_Particle%ExtEffBC(ib) &
                                        +(1.0-bc_vol_frac)*xxx*p_Particle%ExtEffShell(ib)
			  SubScatCoeff(ib) = SubScatCoeff(ib) + &
                                        bc_vol_frac*xxx*p_Particle%ScaEffBC(ib) &
                                        +(1.0-bc_vol_frac)*xxx*p_Particle%ScaEffShell(ib)
                        ENDIF
	              ELSE IF (FLAG .EQ. 2) THEN !Volume-averaged internal mixture
			ExtCoeff(ib) = ExtCoeff(ib) + xxx*p_Particle%ExtEffVA(ib)
			ScatCoeff(ib) = ScatCoeff(ib) + xxx*p_Particle%ScaEffVA(ib)
			AssymSum = AssymSum + xxx*p_Particle%ScaEffVA(ib)*p_Particle%AssymParamVA(ib)

                        !Particles with diameters less than 1 micron
                        IF (2*p_Particle%EffectiveRadius .LT. 1.0e-4) THEN
 			  SubExtCoeff(ib) = SubExtCoeff(ib) + xxx*p_Particle%ExtEffVA(ib)
			  SubScatCoeff(ib) = SubScatCoeff(ib) + xxx*p_Particle%ScaEffVA(ib)
                        ENDIF
                      ELSE IF (FLAG .EQ. 3) THEN !Maxwell-Garnett internal mixture
			ExtCoeff(ib) = ExtCoeff(ib) + xxx*p_Particle%ExtEffMG(ib)
			ScatCoeff(ib) = ScatCoeff(ib) + xxx*p_Particle%ScaEffMG(ib)
			AssymSum = AssymSum + xxx*p_Particle%ScaEffMG(ib)*&
                            p_Particle%AssymParamMG(ib)

                        !Particles with diameters less than 1 micron
                        IF (2*p_Particle%EffectiveRadius .LT. 1.0e-4) THEN
 			  SubExtCoeff(ib) = SubExtCoeff(ib) + xxx*p_Particle%ExtEffMG(ib)
			  SubScatCoeff(ib) = SubScatCoeff(ib) + xxx*p_Particle%ScaEffMG(ib)
                        ENDIF
                      ENDIF

                      p_Particle => p_Particle%Next
		END DO

		!Formula from p.327, eqn 9.103 of Jacobson, 2nd ed. "Fundamentals of Atm. Modeling"
		SingScat(ib) = ScatCoeff(ib)/ExtCoeff(ib)
		SubSingScat(ib) = SubScatCoeff(ib)/SubExtCoeff(ib)                
!		SingScatSep(ib) = ScatCoeffSep(ib)/ExtCoeffSep(ib)
!		SubSingScatSep(ib) = SubScatCoeffSep(ib)/SubExtCoeffSep(ib)   
		!IF(ib .EQ. 419 .OR. ib .EQ. 420) WRITE(*,*) "SingScat: ", SingScat(ib), ScatCoeff(ib), ExtCoeff(ib) 
		
		!Formula from p.329, eqn 9.117 of Jacobson, 2nd ed. "Fundamentals of Atm. Modeling"
		Assym(ib) = AssymSum/ScatCoeff(ib)
!		AssymSep(ib) = AssymSumSep/ScatCoeffSep(ib)
		!IF(ib .EQ. 61) WRITE(*,*) "Assym: ", Assym(1)

		!Convert ExtCoeff and BackScatCoeff from 1/cm to 1/Mm
		ExtCoeff(ib) = 1.0e8*ExtCoeff(ib)
		SubExtCoeff(ib) = 1.0e8*SubExtCoeff(ib)
!		ExtCoeffSep(ib) = 1.0e8*ExtCoeffSep(ib)
!		SubExtCoeffSep(ib) = 1.0e8*SubExtCoeffSep(ib)


		!hbs_frac = hemispherical backscattering fraction
                !See Wiscombe and Grams (1976) J. Atmos. Sci. 33(12), pp. 2440-2451.
                !Integrate Eq. 23 using Henyey-Greenstein phase function from Eq. 16
                assym_sq = assym(ib)**2
                hbs_frac = 0.5*((assym_sq-1)/assym(ib))*(1.0/SQRT(assym_sq+2*assym(ib)+1)-1.0/SQRT(assym_sq+1))
		BackScatCoeff(ib) = hbs_frac*ScatCoeff(ib)
                BackScatCoeff(ib) = 1.0e8*BackScatCoeff(ib)

!                assym_sq = assymSep(ib)**2
!                hbs_frac = 0.5*((assym_sq-1)/assymSep(ib))*(1.0/SQRT(assym_sq+2*assymSep(ib)+1)-1.0/SQRT(assym_sq+1))
!		BackScatCoeffSep(ib) = hbs_frac*ScatCoeffSep(ib)
!                BackScatCoeffSep(ib) = 1.0e8*BackScatCoeffSep(ib)
		!IF(ib .EQ. 61) WRITE(*,*) "ExtCoeff: ", ExtCoeff(1)

	END DO
	
        
	RETURN
END SUBROUTINE AerosolOptProp

!Calculates total aerosol extinction coefficient (m**-1) 
!average single scattering albedo, and 
!average assymetry factor 
!for a grid point at 18 radiative bands.
SUBROUTINE AerosolOptPropFASTTUV (ExtCoeff, SingScat, Assym, BackScatCoeff, &
                          SubExtCoeff, SubSingScat, Flag)!&
                          !ExtCoeffSep, SingScatSep, AssymSep, BackScatCoeffSep, &
                          !SubExtCoeffSep, SubSingScatSep)
! 08/30/2010 Changed for the photolysis wavebands MJA
! 02/27/2012 Changed for photolysis radiation bins MJA	
! 05/03/2012 Changed to calculate backscattering coefficient, units
!             for extcoeff and backscatcoeff now Mm^-1
! 08/16/2012 Fixed backscattering coefficient, added submicron extinction & SSA
! 08/17/2012 Added capability for external mixtures
! 11/08/2012 Added flag to decide between mixing rules
!            Flag = 0, core-in-shell (default)
!            Flag = 1, external mixture of BC and rest
!            Flag = 2, volume average dielectric constant mixing rule
!            Flag = 3, Maxwell-Garnett mixing rule
!
	USE ModelParameters, ONLY : PI, InversionHeight

	IMPLICIT NONE

        !Input Variables
        INTEGER :: Flag
	!Output Variables
	REAL*8, DIMENSION(7) :: ExtCoeff, SingScat, Assym, BackScatCoeff
	REAL*8, DIMENSION(7) :: SubExtCoeff, SubSingScat
        !Below are for external mixture of BC and "shell"
!	REAL*8, DIMENSION(7) :: ExtCoeffSep, SingScatSep, AssymSep, BackScatCoeffSep
!	REAL*8, DIMENSION(7) :: SubExtCoeffSep, SubSingScatSep

	!Internal variables
	INTEGER :: ib, numbin
        LOGICAL :: FASTTUV
	REAL*8, DIMENSION(7) :: ScatCoeff, SubScatCoeff
!	REAL*8, DIMENSION(7) :: ScatCoeffSep, SubScatCoeffSep
	REAL*8 :: AssymSum,xxx, hbs_frac, assym_sq, bc_vol_frac

        ! CMB: Gfortran doesn't like the additional definitions of the 
        !      public variables, nor does it like "Particle" being used
        !      both as a type and a variable name
	Type(Particle),POINTER :: p_Particle


        DO ib = 1, 7
		ExtCoeff(ib) = 0.0
		ScatCoeff(ib) = 0.0
		BackScatCoeff(ib) = 0.0
		SubExtCoeff(ib) = 0.0
		SubScatCoeff(ib) = 0.0
                AssymSum = 0.0
!		ExtCoeffSep(ib) = 0.0
!		ScatCoeffSep(ib) = 0.0
!		BackScatCoeffSep(ib) = 0.0
!		SubExtCoeffSep(ib) = 0.0
!		SubScatCoeffSep(ib) = 0.0
!                AssymSumSep = 0.0
		p_Particle => Particles%First

		DO WHILE (ASSOCIATED(p_Particle))

!                     p_Particle%Temperature = TEMP
                     ! WRITE(*,*)    'cur%Temperature = ',TEMP   
                     !Recalculate optical parameters
                     !WRITE(*,*) "Call Optical, Particle #", I
                     CALL ShellRefIndAndRad(p_Particle)

                     IF (p_Particle%EffectiveRadius .gt. 0) THEN 
                       bc_vol_frac = (p_Particle%AbsCoreRad/&
                               p_Particle%EffectiveRadius)**3
                     ELSE
                       bc_vol_frac = 0.0
                     ENDIF
		     xxx = p_Particle%NumberofParticles*&
                        (PI*p_Particle%EffectiveRadius**2)

                     IF (FLAG .EQ. 0) THEN !Core in Shell  
			ExtCoeff(ib) = ExtCoeff(ib) + &
                            xxx*p_Particle%ExtEff(ib)
			ScatCoeff(ib) = ScatCoeff(ib) + &
                            xxx*p_Particle%ScaEff(ib)
			AssymSum = AssymSum + &
                            xxx*p_Particle%ScaEff(ib)*p_Particle%AssymParam(ib)

                        !Particles with diameters less than 1 micron
                        IF (2*p_Particle%EffectiveRadius .LT. 1.0e-4) THEN
 			  SubExtCoeff(ib) = SubExtCoeff(ib) + &
                            xxx*p_Particle%ExtEff(ib)
			  SubScatCoeff(ib) = SubScatCoeff(ib) + &
                            xxx*p_Particle%ScaEff(ib)
                        ENDIF

                     ELSE IF (Flag .EQ. 1) THEN !External mixture
                        ExtCoeff(ib) = ExtCoeff(ib) + &
                            bc_vol_frac*xxx*p_Particle%ExtEffBC(ib)&
                            +(1.0-bc_vol_frac)*xxx*p_Particle%ExtEffShell(ib)
			ScatCoeff(ib) = ScatCoeff(ib) + &
                            bc_vol_frac*xxx*p_Particle%ScaEffBC(ib) &
                            +(1.0-bc_vol_frac)*xxx*p_Particle%ScaEffShell(ib)
			AssymSum = AssymSum + &
                            bc_vol_frac*xxx*p_Particle%ScaEffBC(ib)*&
                                p_Particle%AssymParamBC(ib) &
                                       +(1.0-bc_vol_frac)*xxx* &
                                       p_Particle%ScaEffShell(ib)*&
                                       p_Particle%AssymParamShell(ib) 
                        !Particles with diameters less than 1 micron
                        IF (2*p_Particle%EffectiveRadius .LT. 1.0e-4) THEN	
                          SubExtCoeff(ib) = SubExtCoeff(ib) + &
                                        bc_vol_frac*xxx*p_Particle%ExtEffBC(ib) &
                                        +(1.0-bc_vol_frac)*xxx*p_Particle%ExtEffShell(ib)
			  SubScatCoeff(ib) = SubScatCoeff(ib) + &
                                        bc_vol_frac*xxx*p_Particle%ScaEffBC(ib) &
                                        +(1.0-bc_vol_frac)*xxx*p_Particle%ScaEffShell(ib)
                        ENDIF
	              ELSE IF (FLAG .EQ. 2) THEN !Volume-averaged internal mixture
			ExtCoeff(ib) = ExtCoeff(ib) + xxx*p_Particle%ExtEffVA(ib)
			ScatCoeff(ib) = ScatCoeff(ib) + xxx*p_Particle%ScaEffVA(ib)
			AssymSum = AssymSum + xxx*p_Particle%ScaEffVA(ib)*p_Particle%AssymParamVA(ib)

                        !Particles with diameters less than 1 micron
                        IF (2*p_Particle%EffectiveRadius .LT. 1.0e-4) THEN
 			  SubExtCoeff(ib) = SubExtCoeff(ib) + xxx*p_Particle%ExtEffVA(ib)
			  SubScatCoeff(ib) = SubScatCoeff(ib) + xxx*p_Particle%ScaEffVA(ib)
                        ENDIF
                      ELSE IF (FLAG .EQ. 3) THEN !Maxwell-Garnett internal mixture
			ExtCoeff(ib) = ExtCoeff(ib) + xxx*p_Particle%ExtEffMG(ib)
			ScatCoeff(ib) = ScatCoeff(ib) + xxx*p_Particle%ScaEffMG(ib)
			AssymSum = AssymSum + xxx*p_Particle%ScaEffMG(ib)*&
                            p_Particle%AssymParamMG(ib)

                        !Particles with diameters less than 1 micron
                        IF (2*p_Particle%EffectiveRadius .LT. 1.0e-4) THEN
 			  SubExtCoeff(ib) = SubExtCoeff(ib) + xxx*p_Particle%ExtEffMG(ib)
			  SubScatCoeff(ib) = SubScatCoeff(ib) + xxx*p_Particle%ScaEffMG(ib)
                        ENDIF
                      ENDIF

                      p_Particle => p_Particle%Next
		END DO

		!Formula from p.327, eqn 9.103 of Jacobson, 2nd ed. "Fundamentals of Atm. Modeling"
		SingScat(ib) = ScatCoeff(ib)/ExtCoeff(ib)
		SubSingScat(ib) = SubScatCoeff(ib)/SubExtCoeff(ib)                
!		SingScatSep(ib) = ScatCoeffSep(ib)/ExtCoeffSep(ib)
!		SubSingScatSep(ib) = SubScatCoeffSep(ib)/SubExtCoeffSep(ib)   
		!IF(ib .EQ. 419 .OR. ib .EQ. 420) WRITE(*,*) "SingScat: ", SingScat(ib), ScatCoeff(ib), ExtCoeff(ib) 
		
		!Formula from p.329, eqn 9.117 of Jacobson, 2nd ed. "Fundamentals of Atm. Modeling"
		Assym(ib) = AssymSum/ScatCoeff(ib)
!		AssymSep(ib) = AssymSumSep/ScatCoeffSep(ib)
		!IF(ib .EQ. 61) WRITE(*,*) "Assym: ", Assym(1)

		!Convert ExtCoeff and BackScatCoeff from 1/cm to 1/Mm
		ExtCoeff(ib) = 1.0e8*ExtCoeff(ib)
		SubExtCoeff(ib) = 1.0e8*SubExtCoeff(ib)
!		ExtCoeffSep(ib) = 1.0e8*ExtCoeffSep(ib)
!		SubExtCoeffSep(ib) = 1.0e8*SubExtCoeffSep(ib)


		!hbs_frac = hemispherical backscattering fraction
                !See Wiscombe and Grams (1976) J. Atmos. Sci. 33(12), pp. 2440-2451.
                !Integrate Eq. 23 using Henyey-Greenstein phase function from Eq. 16
                assym_sq = assym(ib)**2
                hbs_frac = 0.5*((assym_sq-1)/assym(ib))*(1.0/SQRT(assym_sq+2*assym(ib)+1)-1.0/SQRT(assym_sq+1))
		BackScatCoeff(ib) = hbs_frac*ScatCoeff(ib)
                BackScatCoeff(ib) = 1.0e8*BackScatCoeff(ib)

!                assym_sq = assymSep(ib)**2
!                hbs_frac = 0.5*((assym_sq-1)/assymSep(ib))*(1.0/SQRT(assym_sq+2*assymSep(ib)+1)-1.0/SQRT(assym_sq+1))
!		BackScatCoeffSep(ib) = hbs_frac*ScatCoeffSep(ib)
!                BackScatCoeffSep(ib) = 1.0e8*BackScatCoeffSep(ib)
		!IF(ib .EQ. 61) WRITE(*,*) "ExtCoeff: ", ExtCoeff(1)

	END DO
	
        
	RETURN
END SUBROUTINE AerosolOptPropFASTTUV
