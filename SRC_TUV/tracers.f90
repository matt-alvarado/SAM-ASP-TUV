
module tracers

! This module serves as a template for adding tracer transport in the model. The tracers can be 
! chemical tracers, or bin microphysics drop/ice categories, etc. 
! The number of tracers is set by the parameter ntracers which is set in domain.f90.
! Also, the logical flag dotracers should be set to .true. in namelist (default is .false.).
! The model will transport the tracers around automatically (advection and SGS diffusion).
! The user must supply the initialization in the subroutine tracers_init() in this module.
! By default, the surface flux of all tracers is zero. Nonzero values can be set in tracers_flux().
! The local sinks/sources of tracers should be supplied in tracers_physics().

 use grid
 use StepASP, ONLY: InitializeASP, ASPInterface, &
      SAM_wrapper, StepASPOnce, OutputConcentrations
 use Chemistry
 use Aerosols
 use Condensation
 use OutputRoutines
 use Coagulation
 use GridPointFields
 use InfrastructuralCode
 implicit none

 real tracer  (dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, 0:ntracers)
 real tracerfire  (dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, 0:ntracers)
 real traceravg (nx,ny,nzm,ntracers) ! saved for 3D diagnostics - CRL
 real fluxbtr (nx, ny, 0:ntracers) ! surface flux of tracers
 real fluxttr (nx, ny, 0:ntracers) ! top boundary flux of tracers
 real trwle(nz,0:ntracers)  ! resolved vertical flux 
 real trwsb(nz,0:ntracers)  ! SGS vertical flux 
 real tradv(nz,0:ntracers)  ! tendency due to vertical advection 
 real trdiff(nz,0:ntracers)  ! tendency due to vertical diffusion 
 real trphys(nz,0:ntracers)  ! tendency due to physics 
 REAL*8 ExtCoeffArr(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, 7)
 REAL*8 SingScatArr(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, 7)
 REAL*8 AssymArr(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, 7)
 REAL*8 PhotoArr(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, 51)
!REAL*8 :: samtstart, samAOD(nzm,7), samSSA(nzm,7), samASYM(nzm,7), samvalj(nzm,51)
!REAL*8 :: samtstart, samAOD(nzm,7), samSSA(nzm,7), samASYM(nzm,7), samvalj(51,nzm) !AD DEBUG
! use same array dimensions as TUV4SAMASP
REAL*8 :: samtstart, samAOD(200,8), samSSA(200,8), samASYM(200,8), samvalj(150,200) !MM DEBUG
 character *20 tracername(0:ntracers)
 character *10 tracerunits(0:ntracers)
 character *20 aero_name
 real tracerbound(nzm,0:ntracers) ! boundary conditions for tracers
                                  ! if dotracerperiodic is used
      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      REAL*8 :: BOXVOL, BOXMASS, TEMPTM, BOXAREA
      REAL*8 :: PRESTM, RHTM 

      TYPE(Particle),    POINTER :: Current, Env

      ! Scalars   
      INTEGER, PARAMETER :: ICOMP = 3     
      INTEGER, PARAMETER :: IBINS = 1 ! Changed this from 15 KS
      INTEGER, PARAMETER :: IDIAG = 2
      INTEGER, PARAMETER :: SRTSO4 = 1
      INTEGER, PARAMETER :: SRTNH4 = 2
      INTEGER, PARAMETER :: SRTH2O = 3
      INTEGER, PARAMETER :: ngas = 617 ! # of gas species before aerosol in 
                                     ! tracer array


      !INTEGER, PARAMETER :: HowManyEvolveGasChems = 617

      INTEGER, PARAMETER :: HowManyBins = 10

      ! Tracer classification
      integer, parameter :: tbulk = 1 ! bulk tracers
      integer, parameter :: tnum = 2 ! aerosol number tracers
      integer, parameter :: tpmass = 3 ! aerosol prognostic mass tracers
      integer, parameter :: tdmass = 4 ! aerosol diagnostic mass tracers
      integer, parameter :: ipmass = (icomp-idiag)*ibins ! number of prog mass
                                             ! tracers
      integer, parameter :: idmass = idiag*ibins ! num of diagnostic mass tr.
      integer tclass(ntracers)
      integer tbin(ntracers) ! which bin does tracer correspond to (0 for non-TOMAS)
                             ! initialized in tracers_init
      real fluxsave(-1:nxp3,-1:nyp3,nz,2,3,ibins) ! used for advection
      real fnumsave(-2:nxp3,-2:nyp3,nz,2,ibins) ! used for advection
      real tracervar (nx,ny,nzm,ngas) ! variance saved for 3D diagnostics

	! Emissions - CRL

      real*8 :: emisgridy ! y-grid to emit stuff from
      real*8 :: xextent   ! length of fire, in boxes needed for stuff KS
      real*8 :: yextent   ! y extent of fire in boxes
      real*8 :: zextent   ! z extent of fire in boxes
      real*8 :: zlower    ! determines emission range in x for fire KS
      real*8 :: zhigher   ! high bound of box KS
      real*8 :: Mgas ! Air density from ASP units?
      real*8 :: tv ! days since vernal equinox
      real*8 :: lambda ! longitude of Earth from vernal equinox

      integer ASP_time_run
      real*8, parameter :: airmw = 28.97 ! molecular weight of air [g/mol]
      real*8, parameter :: mol_O3= 0.21 ! molecular weight of air [kg/mol]
      real*8, parameter :: gasc = 8.3144598 ! gas constant [m3 pa K-1 mol-1]
      real*8, parameter :: Av = 6.022e23 ! Avagadro's number [molec mol-1]
      real*8, parameter :: ug_to_kg = 1e-9 ! convert ug to kg
      real*8, parameter :: kg_to_ug = 1e9 ! convert kg to ug
      real*8, parameter :: m3_to_cm3 = 1e-6 ! convert m-3 to cm-3
      real*8, parameter :: cm3_to_m3 = 1e6 ! convert cm-3 to m-3
      real*8, parameter :: g_to_kg = 1e-3
      real*8, parameter :: pi = 3.1415926535897931
      real*8, parameter :: vernal_equinox  = 79 ! day of year
      real*8, parameter :: earth_tilt = 23.5 ! earth's tilt (degrees)

      REAL*8, DIMENSION(7):: ExtCoeff, SingScat, Assym
CONTAINS

FUNCTION Find_Chem(Chem) RESULT (Index)

  USE domain, ONLY: ntracers
  IMPLICIT NONE
  INTEGER				  :: i, Index, Phase
  CHARACTER*(*)		  :: Chem

  DO i = 1, ntracers
     IF (Trim(Chem) .eq. Trim(tracername(i))) THEN
        Index = i
        EXIT
     END IF
  ENDDO

  RETURN

END FUNCTION

 subroutine tracers_init()
  use vars, only: P, tabs
  use params
  use ModelParameters,	ONLY : HowManyBins
  use InfrastructuralCode	
  use GridPointFields
  use Chemistry
  use Aerosols		
  use Condensation
  use OutputRoutines
  use Coagulation   
  integer i,j,k,l,ntr,ntrtrack,NumAqChems ,NumOrgChems,NumAqOrgChems,NumAq,PART,B
  integer binnum,temp_num,aerosol_track

  Integer :: 	KNO3i, KCli, KHSO4i, &
			Nai, NaNO3i, NaCli, NaHSO4i, Na2SO4i, &
			NH3i, NH4i, NH4NO3i, NH4Cli, NH4HSO4i, NH42SO4i, &
			Cli, NO3i, SO4i, HSO4i, K2SO4i, PotassiumIndex, &
			K2SO4Index, FH, NumBins	, COi, Traceri, Cai, &
                        KOHi, NaOHi, Mgi, MgNO32i, MgCl2i, MgHSO42i, &
                        MgSO4i, MgOH2i, CaNO32i, CaCl2i, CaHSO42i, CaSO4i, &
                        CaOH2i, LEVi, CO2i

  character *2 ntrchar
  character *2 tmp_name

  integer, external :: lenstr
  real*8 :: airdens ! density of air [kg/m3]
  real*8 :: airdenstm ! density of air [molec/cm3]
  real*8 :: part_conc ! particle concentration [#/kg_burned]
  real*8 :: CO_conc ! CO concentration [g/cm3]
  TYPE(Particle), POINTER :: env, cur

      !! ------------- from ASP ----------------- External Variables
		REAL*8 :: temp_val
		REAL*8 :: Timestep !(in s)
		REAL*8 :: Temp  !(in K)
		REAL*8 :: Press !(in mbar)
		REAL*8 :: Dens  !(in kg/m3 dry air)
		REAL*8 :: GasConc(616) ! (in ppb)
		REAL*8 :: AeroNumConc(10) !(in # particles/cm3)
		REAL*8 :: AeroMassConc(10, 86) !(in ug/m3)
		REAL*8 :: AeroNumConcCheck(10) !(in # particles/cm3)
		REAL*8 :: AeroMassConcCheck(10, 86) !(in ug/m3)
		REAL*8 :: Radius(10)
		REAL*8 :: TermVel(10)
		REAL*8 :: EnvNumConc(10) !(in # particles/cm3)
		REAL*8 :: EnvMass(10, 86) !(in ug/m3)
                ! CB: Ifort complains unless these dimensions match
		!REAL*8, DIMENSION(18) :: ExtCoeff, SingScat, Assym
                REAL*8, DIMENSION(7):: ExtCoeff, SingScat, Assym
                REAL*8, allocatable :: SaltConc(:)
                REAL*8, allocatable :: SaltConcEnv(:)
		REAL*8 :: filler_value
		REAL*8 :: NumConcWeight(10)
		REAL*8 :: MassConcWeight(10)
	 	integer :: doAerosols
                real*8,parameter :: CO_init = 10000.0 ! ppb
                integer, parameter :: CO_index = 2
      ! ---------------------- end from ASP

 tracer = 0.
 tracerfire = 0.
 traceravg = 0. ! - CRL
 fluxbtr = 0.
 fluxttr = 0.
 filler_value = 1.0E-34
 doAerosols = 1
! Add your initialization code here. Default is to set to 0 in setdata.f90.

! if(nrestart.eq.0) then

! here ....

! endif

  call InitializeASP(doAerosols)
  call ReadBoundaryDistributionsOrganic
  call PopulateParticlesSectionsRightAway(.TRUE.)
  call SortBoundaryPart
  write(*,*) 'After InitializeASP.... '
! Initialize values
  Temp = GetTemp()
  Press = GetPress()
  Dens = Press*airmw/(gasc*Temp) ! in kg/m3
  Mgas = GetM()
  ExtCoeffArr = 1e-10 !AD DEBUG AOD !mja change 20210722
  SingScatArr = 0.99 !AD DEBUG SSA !mja change 20210722
  AssymArr = 0.61 !AD DEBUG ASYM !mja change 20210722
  samAOD = 1e-10 !MM DEBUG
  samSSA = 0.99 !MM DEBUG
  samASYM = 0.61 !MM DEBUG

IF(doAerosols.gt.0) then

   NumAq = HowManyAqChems+HowManyAqCations+HowManyAqAnions
  ALLOCATE (SaltConc(HowManyAqChems))
   !current => particles%first
   Env => boundaryparticles%first
   cur => particles%first
   NumAq = HowManyAqChems+HowManyAqCations+HowManyAqAnions
   I = 1
   DO WHILE(associated(Env))

     EnvNumConc(I) = Env%numberofparticles
     AeroNumConcCheck(I) = cur%numberofparticles
     SaltConc = HypotheticalElectrolyteConcentrations(Env)
     EnvMass(I,1) = Env%AqChems(1)*(&
		Env%Numberofparticles*AqMolecularMass(1))/(1.0e-12)

     write(*,*) 'Edges of bin ',I,' = ',Env%Edges, '(cm)'
     !Force minimum water
     IF(EnvMass(I,1) .LT. 0.001) EnvMass(I,1) = 0.001
 
     IF(Env%numberofparticles .GT. filler_value) then !0.) THEN		
	DO J = 2, HowManyAqChems
	   if (J .EQ. findchem('NH3', 1)) then
	      EnvMass(I,J) = Env%AqChems(J)*Env%Numberofparticles*AqMolecularMass(J)*1.0e-12	
	      AeroMassConcCheck(I,J) = cur%AqChems(J)*cur%Numberofparticles*AqMolecularMass(J)*1.0e-12			
	   else 
	      EnvMass(I,J) = (SaltConc(J)+Env%AqChems(J))*Env%Numberofparticles*AqMolecularMass(J)*1.0e-12
	      AeroMassConcCheck(I,J) = (SaltConc(J)+cur%AqChems(J))*cur%Numberofparticles*AqMolecularMass(J)*1.0e-12
	   endif
	ENDDO
	DO J = HowManyAqChems+1, NumAq
		EnvMass(I,J) = 0.
	ENDDO
				
	DO J = 1, HowManyOrgChems
		EnvMass(I,J+NumAq) = Env%OrgChems(J)*Env%Numberofparticles*OrgMolecularMass(J)*1e-12
		AeroMassConcCheck(I,J+NumAq) = cur%OrgChems(J)*cur%Numberofparticles*OrgMolecularMass(J)*1e-12
                write(*,*) 'Declaring Name ', OrgPhaseChemicalNames(J), J+NumAq
	ENDDO
				
	DO J = 1, HowManyAqOrgChems
		EnvMass(I,J+NumAq+HowManyOrgChems) = Env%AqOrgChems(J)*Env%Numberofparticles*AqOrgMolecularMass(J)*1e-12
		AeroMassConcCheck(I,J+NumAq+HowManyOrgChems) = cur%AqOrgChems(J)*cur%Numberofparticles*AqOrgMolecularMass(J)*1e-12
	ENDDO

    ELSE !No particles
	EnvNumConc(I) = 0.
	DO J = 2, NumAq+HowManyOrgChems+HowManyAqOrgChems
 		EnvMass(I,J) = 0.
	ENDDO			
    ENDIF
    I = I + 1
		
    Env => Env%next
    cur => cur%next

    ENDDO !END STEP PARTICLE DILUTION



  CALL OutputConcentrations(Temp, Press, Dens, GasConc, &
					AeroNumConc, AeroMassConc, &
					ExtCoeff, SingScat, Assym, Radius, TermVel)


!AD DEBUG remove rank.eq.5 criterion
if (rank.eq.0) then
!if (rank.eq.5) then
   write(*,*) 'Org Aerosols out of ASP'
do j=1, 23
   write(*,*) OrgPhaseChemicalNames(j)
   write(*,*) 'EnvMass(:,j) = ',EnvMass(:,j+NumAq)
   write(*,*) 'AeroMassConcCheck(:,j) = ',AeroMassConcCheck(:,j+NumAq)
   write(*,*) 'AeroMass(:,j) = ',AeroMassConc(:,j+NumAq)
enddo 
endif
endif! do aerosols
! Specify the tracers' default names:

   ! Default names are TRACER01, TRACER02, etc:

   do ntr = 1, HowManyEvolveGasChems! - 616 evolving Gas species CRL
   	tracername(ntr) = GasPhaseChemicalNames(ntr)
   	tracerunits(ntr) = '[ppb]'
   enddo

NumAq = HowManyAqChems+HowManyAqCations+HowManyAqAnions
write(*,*) NumAq
write(*,*) 'HowManyEvolveGasChems = ',HowManyEvolveGasChems
write(*,*) 'HowManyAqChems = ',HowManyAqChems
write(*,*) 'HowManyAqCations = ',HowManyAqCations
write(*,*) 'HowManyAqAnions = ',HowManyAqAnions
write(*,*) 'HowManyOrgChems = ',HowManyOrgChems
write(*,*) 'HowManyAqOrgChems = ',HowManyAqOrgChems

	!Find chemical indices for later
	PotassiumIndex = FindChem("K+", 1)
 
	KNO3i = FindChem("KNO3", 1)
	KCli = FindChem("KCl", 1)
	KHSO4i = FindChem("KHSO4", 1)
	K2SO4Index = FindChem("K2SO4", 1)
	K2SO4i = K2SO4Index
        KOHi = FindChem("KOH", 1)
		
	Nai = FindChem("Na+", 1)
	NaNO3i = FindChem("NaNO3", 1)
	NaCli = FindChem("NaCl", 1)
	NaHSO4i = FindChem("NaHSO4", 1)
	Na2SO4i = FindChem("Na2SO4", 1)
        NaOHi = FindChem("NaOH", 1)

	Mgi = FindChem("Mg++", 1)
	MgNO32i = FindChem("Mg(NO3)2", 1)
	MgCl2i = FindChem("MgCl2", 1)
	MgHSO42i = FindChem("Mg(HSO4)2", 1)
	MgSO4i = FindChem("MgSO4", 1)
        MgOH2i = FindChem("Mg(OH)2", 1)

	Cai = FindChem("Ca++", 1)
	CaNO32i = FindChem("Ca(NO3)2", 1)
	CaCl2i = FindChem("CaCl2", 1)
	CaHSO42i = FindChem("Ca(HSO4)2", 1)
	CaSO4i = FindChem("CaSO4*2H2O", 1)
        CaOH2i = FindChem("Ca(OH)2", 1)
		
	NH3i = FindChem("NH3", 1)
	NH4i = FindChem("NH4+", 1)
	NH4NO3i = FindChem("NH4NO3", 1)
	NH4Cli = FindChem("NH4Cl", 1)
	NH4HSO4i = FindChem("NH4HSO4", 1)
	NH42SO4i = FindChem("(NH4)2SO4", 1)
	LEVi = FindChem("(NH4)3H(SO4)2", 1)

	Cli = FindChem("Cl-", 1)
	NO3i = FindChem("NO3-", 1)
	SO4i = FindChem("SO4--", 1)
	HSO4i = FindChem("HSO4-", 1)

	Traceri = FindChem("Trace", 0)
	COi = FindChem("CO", 0)
	CO2i = FindChem("CO2", 0)

	!Printing chemical indices for later
	write(*,*) 'Printing chemical indices for later'
	write(*,*) PotassiumIndex, AqPhaseChemicalNames(PotassiumIndex), AqMolecularMass(PotassiumIndex)
 
	write(*,*)	KNO3i , AqPhaseChemicalNames(KNO3i), AqMolecularMass(KNO3i)
	write(*,*)	KCli , AqPhaseChemicalNames(KCli), AqMolecularMass(KCli)
	write(*,*)	KHSO4i, AqPhaseChemicalNames(KHSO4i), AqMolecularMass(KHSO4i)
	write(*,*)	K2SO4Index , AqPhaseChemicalNames(K2SO4Index), AqMolecularMass(K2SO4Index)
	write(*,*)	K2SO4i , AqPhaseChemicalNames(K2SO4i), AqMolecularMass(K2SO4i)
 	write(*,*)      KOHi , AqPhaseChemicalNames(KOHi), AqMolecularMass(KOHi)
		
	write(*,*)	Nai , AqPhaseChemicalNames(Nai), AqMolecularMass(Nai)
	write(*,*)	NaNO3i , AqPhaseChemicalNames(NaNO3i), AqMolecularMass(NaNO3i)
	write(*,*)	NaCli , AqPhaseChemicalNames(NaCli), AqMolecularMass(NaCli)
	write(*,*)	NaHSO4i , AqPhaseChemicalNames(NaHSO4i), AqMolecularMass(NaHSO4i)
	write(*,*)	Na2SO4i , AqPhaseChemicalNames(Na2SO4i), AqMolecularMass(Na2SO4i)
	write(*,*)      NaOHi , AqPhaseChemicalNames(NaOHi), AqMolecularMass(NaOHi)

	write(*,*)	Mgi, AqPhaseChemicalNames(Mgi), AqMolecularMass(Mgi)
	write(*,*)	MgNO32i , AqPhaseChemicalNames(MgNO32i), AqMolecularMass(MgNO32i)
	write(*,*)	MgCl2i , AqPhaseChemicalNames(MgCl2i), AqMolecularMass(MgCl2i)
	write(*,*)	MgHSO42i, AqPhaseChemicalNames(MgHSO42i), AqMolecularMass(MgHSO42i)
	write(*,*)	MgSO4i , AqPhaseChemicalNames(MgSO4i), AqMolecularMass(MgSO4i)
	write(*,*)      MgOH2i , AqPhaseChemicalNames(MgOH2i), AqMolecularMass(MgOH2i)

	write(*,*)	Cai, AqPhaseChemicalNames(Cai), AqMolecularMass(Cai)
	write(*,*)	CaNO32i, AqPhaseChemicalNames(CaNO32i), AqMolecularMass(CaNO32i)
	write(*,*)	CaCl2i , AqPhaseChemicalNames(CaCl2i), AqMolecularMass(CaCl2i)
	write(*,*)	CaHSO42i , AqPhaseChemicalNames(CaHSO42i), AqMolecularMass(CaHSO42i)
	write(*,*)	CaSO4i , AqPhaseChemicalNames(CaSO4i), AqMolecularMass(CaSO4i)
	write(*,*)      CaOH2i , AqPhaseChemicalNames(CaOH2i), AqMolecularMass(CaOH2i)
		
	write(*,*)	NH3i, AqPhaseChemicalNames(NH3i), AqMolecularMass(NH3i)
	write(*,*)	NH4i, AqPhaseChemicalNames(NH4i), AqMolecularMass(NH4i)
	write(*,*)	NH4NO3i, AqPhaseChemicalNames(NH4NO3i), AqMolecularMass(NH4NO3i)
	write(*,*)	NH4Cli , AqPhaseChemicalNames(NH4Cli), AqMolecularMass(NH4Cli)
	write(*,*)	NH4HSO4i , AqPhaseChemicalNames(NH4HSO4i), AqMolecularMass(NH4HSO4i)
	write(*,*)	NH42SO4i , AqPhaseChemicalNames(NH42SO4i), AqMolecularMass(NH42SO4i)
	write(*,*)	LEVi , AqPhaseChemicalNames(LEVi), AqMolecularMass(LEVi)

	write(*,*)	Cli , AqPhaseChemicalNames(Cli), AqMolecularMass(Cli)
	write(*,*)	NO3i , AqPhaseChemicalNames(NO3i), AqMolecularMass(NO3i)
	write(*,*)	SO4i, AqPhaseChemicalNames(SO4i), AqMolecularMass(SO4i)
	write(*,*)	HSO4i , AqPhaseChemicalNames(HSO4i), AqMolecularMass(HSO4i)

if(doAerosols.gt.0) then
! Aqueous Phase Chems
ntrtrack = HowManyEvolveGasChems+1
do binnum=1,HowManyBins
   do ntr = 1, NumAq! - 41 total AqPhaseChemicals  CRL

		if (binnum.lt.10) then
			write(tmp_name,'(i1)') binnum
			aero_name = trim(AqPhaseChemicalNames(ntr))//'0'//trim(tmp_name)
		else
			write(tmp_name,'(i2)') binnum
			aero_name = trim(AqPhaseChemicalNames(ntr))//trim(tmp_name)
		endif

   		tracername(ntrtrack) = aero_name
   		tracerunits(ntrtrack) = '[kg/kg_air]'
		ntrtrack = ntrtrack+1
   enddo


write(*,*)'!!!!!!!!!!!! '
write(*,*)'Organic Phase Chems Names '
! Organic Phase Chems
   do ntr = 1,HowManyOrgChems	
		if (binnum.lt.10) then
			write(tmp_name,'(i1)') binnum
			aero_name = trim(OrgPhaseChemicalNames(ntr))//'0'//trim(tmp_name)
		else
			write(tmp_name,'(i2)') binnum
			aero_name = trim(OrgPhaseChemicalNames(ntr))//trim(tmp_name)
		endif

   		tracername(ntrtrack) = aero_name
   		tracerunits(ntrtrack) = '[kg/kg_air]'
		ntrtrack = ntrtrack+1

   enddo

! Aqueous Organic Phase Chems

   do ntr = 1,HowManyAqOrgChems	

		if (binnum.lt.10) then
			write(tmp_name,'(i1)') binnum
			aero_name = trim('Aq'//AqOrgPhaseChemicalNames(ntr))//'0'//trim(tmp_name)
		else
			write(tmp_name,'(i2)') binnum
			aero_name = trim('Aq'//AqOrgPhaseChemicalNames(ntr))//trim(tmp_name)
		endif

   		tracername(ntrtrack) = aero_name
   		tracerunits(ntrtrack) = '[kg/kg_air]'
		ntrtrack = ntrtrack+1

   enddo

! NumConc



		if (binnum.lt.10) then
			write(tmp_name,'(i1)') binnum
			aero_name = 'NoCo'//'0'//trim(tmp_name)
		else
			write(tmp_name,'(i2)') binnum
			aero_name = 'NoCo'//trim(tmp_name)
		endif

   		tracername(1476+binnum) = aero_name
   		tracerunits(1476+binnum) = '[#/kg_air]'

enddo ! binnum

endif ! do aerosols

do i=1,nx
   do j=1,ny
      do l=1,nzm
            prestm = dble(pres(l)*100.d0) +p(i,j,l)! in Pa  
 
            temptm = dble(tabs(I,J,L)) ! in K
            boxarea = dble(dx*dy) ! in m2

            if (l.eq.1)then
               boxvol = dble(dx*dy*dz*1.0d6) ! in cm3
            else
               boxvol = dble(dx*dy*(zi(l)-zi(l-1))*1.0d6) ! in cm3
            endif

            airdens = prestm*(airmw*1e-3)/(gasc*temptm) ! in kg_air/m3
            airdenstm = prestm*Av/(gasc*temptm)/1e6 ! in molec/cm3
            boxmass = dble(airdens*(boxvol/1.0d6)) ! in kg_air!

            !massemis kg burned m-2 s-1
            ntrtrack = 1

            !Gas Phase Chems
            do ntr=1,HowManyEvolveGasChems

                ! assuming dry air
		tracer(i,j,l,ntr)= (GasPhaseBackground(ntr)*GasMolecularMass(ntr))/(Mgas*airmw) !kg/kg_air
  
                ! GasPhaseEmisFactor ! g/kg_burned
                ! massemis ! kg m-2 s-1 
		tracerfire(i,j,l,ntr) = (GasPhaseEmisFactor(ntr)*massemis*boxarea*dt*g_to_kg/boxmass)!kg/kg_air
		ntrtrack = ntrtrack+1

            enddo

            ! Initializing Aerosols
            CO_conc = (CO_init*1e-9 * GasMolecularMass(CO_index) / Av )*airdenstm ! ppbCO to g/cm3

            cur => particles%first
            env => boundaryparticles%first
            write(*,*)'Initializing Aerosols from ASP to SAM'
            write(*,*)'CO_conc =',CO_conc

            B = 1

            DO WHILE(associated(cur))
               SaltConc = HypotheticalElectrolyteConcentrations (cur)
               !SaltConcEnv =  HypotheticalElectrolyteConcentrations (env)

               part_conc = GasPhaseEmisFactor(CO_index)*cur%numberofparticles/CO_conc ! #/kg_burned

               ! water
               tracer(i,j,l,ntrtrack) = (env%AqChems(1))*(env%numberofparticles*AqMolecularMass(1))/(Mgas/Av*airmw)	! convert mol/cm3 to kg/kg_air
               tracerfire(i,j,l,ntrtrack) = (cur%AqChems(1))*AqMolecularMass(1)*part_conc*massemis*boxarea*dt*g_to_kg/boxmass
               ntrtrack = ntrtrack+1

               ! AqChems
               DO PART = 2, HowManyAqChems
                  if (PART .EQ. findchem('NH3', 1)) then
                     tracer(i,j,l,ntrtrack) = (env%AqChems(PART))*(env%numberofparticles*AqMolecularMass(PART))/(Mgas/Av*airmw)	! convert mol/cm3 to kg/kg_air
                     tracerfire(i,j,l,ntrtrack) = (cur%AqChems(PART))*AqMolecularMass(PART)*part_conc*massemis*boxarea*dt*g_to_kg/boxmass	
                  else
                     tracer(i,j,l,ntrtrack) = (env%AqChems(PART))*(env%numberofparticles*AqMolecularMass(PART))/(Mgas/Av*airmw)	! convert mol/cm3 to kg/kg_air
                     tracerfire(i,j,l,ntrtrack) = (SaltConc(PART)+cur%AqChems(PART))*AqMolecularMass(PART)*part_conc*massemis*boxarea*dt*g_to_kg/boxmass

                  end if ! Be careful with NH3 to avoid double counting

                  ntrtrack = ntrtrack+1
               ENDDO

               ! Anions and Cations
               !write(*,*) 'Anions and Cations'
               DO PART = HowManyAqChems+1, NumAq
                  tracer(i,j,l,ntrtrack) = 0.0
                  tracerfire(i,j,l,ntrtrack) = 0.0
                  ntrtrack = ntrtrack+1
               ENDDO
               
               !OrgChems
               !write(*,*) 'OrgChems'
               DO PART = 1, HowManyOrgChems
                  tracer(i,j,l,ntrtrack) = (env%OrgChems(PART))*(env%numberofparticles*OrgMolecularMass(PART))/(Mgas/Av*airmw)! convert ug/m3 to kg/kg_air
                  tracerfire(i,j,l,ntrtrack) = cur%OrgChems(PART)*OrgMolecularMass(PART)*part_conc*massemis*boxarea*dt*g_to_kg/boxmass
                  ntrtrack = ntrtrack+1
               ENDDO
        
               ! AqOrgChems
               !write(*,*) 'AqOrgChems'
               DO PART =1, HowManyAqOrgChems
                  tracer(i,j,l,ntrtrack) = (env%AqOrgChems(PART))*(env%numberofparticles*AqOrgMolecularMass(PART))/(Mgas/Av*airmw)! convert ug/m3 to kg/kg_air
                  tracerfire(i,j,l,ntrtrack) = cur%AqOrgChems(PART)*AqOrgMolecularMass(PART)*part_conc*massemis*boxarea*dt*g_to_kg/boxmass
                  ntrtrack = ntrtrack+1          
               ENDDO

               ! NumConc
               !write(*,*) 'NumConc'
               tracer(i,j,l,1476+B) = env%numberofparticles/(Mgas/Av*airmw*1e-3) ! convert #/cm3 to #/kg_air
               tracerfire(i,j,l,1476+B) = part_conc*massemis*boxarea*dt/boxmass ! convert #/kg_burned to #/kg_air

               B = B + 1
		
               cur => cur%next
               env => env%next				
            ENDDO !END STEP PARTICLE INITIALIZATION

      enddo ! nz
   enddo ! ny
enddo ! nx


   if(dotracerperiodic) then
   else
      do l=1,nzm
         do ntr=1,ntracers
            tracerbound(l,ntr)=tracer(1,1,l,ntr)
         enddo
      enddo
   endif

 end subroutine tracers_init

 subroutine tracers_flux()

! Set surface and top fluxes of tracers. Default is 0 set in setdata.f90

 end subroutine tracers_flux

 !---------------------------------
 subroutine tracers_physics()

 ! add here a call to a subroutine that does something to tracers besides advection and diffusion.
 ! The transport is done automatically. 
  use vars
  use GridPointFields
  use OutputRoutines
  use grid, only: day
  use params
  use Chemistry
      implicit none

      integer i,j,k,l,ntr,ntrtrack,mm,iup,idown,iden,zden,aerosol_track,binnum,ASP_time,ntrtrack2,irxn,nwvn,kk
      integer :: doAerosols

      REAL*8 :: prestm,rhtm, temptm, airdens, boxmass,  Transmissivity,np,svp,airdenstm,tstep,solar_zenith_angle,cur_dist      
      REAL*8 :: dels1,H
      REAL*8 :: airmols,boxmols
      REAL*8 :: wvmixsat,wvmixamb,waterconc,rhtm1
      REAL*8 :: track_totoa
      REAL*8 :: dels ! solar declination angle (seasonal effect) (radians)
      REAL*8 :: hour ! hour of day (UTC time)
      REAL*8 :: iday ! floor of day of year (i.e. no fractional part of day)
      REAL*8 :: sinlea ! sin of local elevation angle

      !! ------------- from ASP ----------------- External Variables!
      REAL*8 :: GasConc(HowManyEvolveGasChems)  
      REAL*8 :: AfterGasConc(HowManyEvolveGasChems)  
      REAL*8 :: TempGasConc(HowManyEvolveGasChems)   
      REAL*8 :: NumConc(10)   
      REAL*8 :: MassConc(10,HowManyAqChems+HowManyAqCations+HowManyAqAnions+HowManyOrgChems+HowManyAqOrgChems)   
      REAL*8 :: AfterNum(10)   
      REAL*8 :: AfterMass(10,HowManyAqChems+HowManyAqCations+HowManyAqAnions+HowManyOrgChems+HowManyAqOrgChems)

      !REAL*8 :: TermVel(HowManyEvolveGasChems)  !!MM DEBUG
      REAL*8 :: TermVel(10) 
      REAL*8 :: t_TotOA   
      REAL*8 :: PhotolysisRates(51)   
      ! ---------------------- end from ASP


      Mgas = GetM()
      !	Enter Fire Emissions, when over fire - else use tracers background concentrations
      zlower = 1!int(round(float(vals[6])*1000./float(dz))) ! this is number of the lower box
      zhigher = zlower+1 ! this is the upper INCLUSIVE box
      emisgridy = 1
      idown=(101-xextent)/2
      iup=(100+xextent)/2
      AfterGasConc = 0.0D0
      Transmissivity = 1.0
      tstep = 10 

      !doAerosols = 1

      ! calculate solar zenith angle
      iday = int(day)
      hour = (day - iday)*24.
      tv = iday-vernal_equinox ! days since vernal equinox  
      lambda = 360*(tv/365.24)! longitude of Earth from vernal equinox
      dels = asin(sin(earth_tilt/180*pi)*sin(lambda/180*pi))
      H = 360*((hour-12)/24.)
      solar_zenith_angle = acos(sin(latitude0/180*pi)*sin(dels)+cos(latitude0/180*pi)*cos(dels)*cos(H/180*pi))*180/pi	


      !cur_dist = 100-((plume_y-dy)/vg/dt)
      !cur_dist = 100+(((plume_y-dy)/vg/dt)*2)
      cur_dist = (((plume_y-dy)/vg/dt)*2)
      write(*,*)'In Tracer Physics:solar cur_dist = ',cur_dist
      !if (cur_dist.lt.plume_y-dy) then


      ! Adding Fire Emissions if over fire
       !if (nstep.gt.cur_dist.and.nstep.lt.100) then  
       if (nstep.gt.0.and.nstep.le.cur_dist) then 
         write(*,*) 'cur_dist = ',cur_dist,' plume_y = ',plume_y,nstep
         !write(*,*) 'rank L681 tracers_physics= ',rank !AD DEBUG
!AD DEBUG remove rank criterion as this seems to stay at 0?
               if (rank.eq.0) then !AD DEBUG
!         if (rank.eq.5) then! only read in plume in the first slab of the 12 x-subdomains)
	    !do iden = 100,101
	    do iden = 10,11
		do zden = 29,32
		     prestm = dble(pres(zden)*100.d0) +p(iden,emisgridy,zden)! in Pa  
	             temptm = dble(tabs(iden,emisgridy,zden)) ! in K

	             airdens = prestm*airmw/(gasc*temptm) ! in kg/m3
                     airdenstm = prestm*Av/(gasc*temptm)/1e6 ! in molec/cm3

		     ntrtrack = 1
                     write(*,*) ' In fire initialization box',iden,zden,nstep
		     do ntr=1, HowManyEvolveGasChems
                        write(*,*) tracername(ntrtrack), tracer(iden,emisgridy,zden,ntrtrack)*airmw/GasMolecularMass(ntrtrack)*1e9! kg/kg_air to ppb
               		tracer(iden,emisgridy,zden,ntrtrack) = tracer(iden,emisgridy,zden,ntrtrack)+tracerfire(iden,emisgridy,zden,ntrtrack)! kg/kg_air! 

                        write(*,*) tracername(ntrtrack), tracer(iden,emisgridy,zden,ntrtrack)*airmw/GasMolecularMass(ntrtrack)*1e9! kg/kg_air to ppb

			ntrtrack = ntrtrack+1
   
		     enddo 

                     do binnum=1,HowManyBins
                        ! Aerosol MassConc
                        ! Aqueous Phase Chems + Organic Phase Chems + Aq Organic Phase Chems
                      do ntr = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions+HowManyOrgChems+HowManyAqOrgChems!

	   		   tracer(iden,emisgridy,zden,ntrtrack) = tracer(iden,emisgridy,zden,ntrtrack) + tracerfire(iden,emisgridy,zden,ntrtrack)! in kg/kg_air
                       !if (iden.eq.12.and.zden.eq.30) then
                           write(*,*) tracername(ntrtrack), binnum,tracer(iden,emisgridy,zden,ntrtrack)*(airdenstm/Av*airmw)*1e12
                       !endif
		       ntrtrack = ntrtrack+1
                      enddo

                       ! NumConc
                      tracer(iden,emisgridy,zden,1476+binnum) = tracer(iden,emisgridy,zden,1476+binnum)+ tracerfire(iden,emisgridy,zden,1476+binnum)! in #/kg_air
                      !if (i.eq.12.and.l.eq.30) then
                        write(*,*) tracername(1476+binnum), binnum,tracer(iden,emisgridy,zden,1476+binnum)*(airdenstm/Av*airmw)*1e-3
                        write(*,*)''
                      !endif
                     enddo ! binnum

		     !end if ! do aerosol
!---------------------------     
	   	enddo ! zden
	    enddo ! iden
         endif ! rank eq 5
      endif ! nstep eq.

   do i=1,nx
   !write (*,*) "i,nx", i, nx  !AD DEBUG
      do j=1,ny
      !write (*,*) "j,ny", j, ny  !AD DEBUG
!         do l=1,nz
         do l=1,nzm  !AD DEBUG
         !write (*,*) "l,nzm", l, nzm !AD DEBUG

            prestm = dble(pres(l)*100.d0) + p(i,j,l) ! in Pa     

            temptm = dble(tabs(I,J,L)) ! in K
            if (l.eq.1)then
               boxvol = dble(dx*dy*dz*1.0d6) ! in cm3
            else
               boxvol = dble(dx*dy*(zi(l)-zi(l-1))*1.0d6) ! in cm3
            endif
            airdens = prestm*airmw/(gasc*temptm) ! in kg/m3
	    airdenstm = prestm*Av/(gasc*temptm)/1.0d6! in molec/cm3

            boxmass = dble(boxvol/1.0d6*airdens) ! in kg!

            wvmixsat = dble(qsatw(real(temptm),real(prestm)/100.)) ! kg water/kg air
            wvmixamb = dble(qv(i,j,l)) ! kg water/kg air
            rhtm = wvmixamb/wvmixsat ! relative humidity 0-1
            svp = (6.112*exp(6816*((1/273.15)-(1/temptm))+(5.1309*log(273.15/temptm))))*100.d0 ! in Pa
	    waterconc = (svp*airdenstm)*rhtm/(prestm) ! 
	    airmols = (pres(l)+p(i,j,l))*100.d0/(gasc*tabs(i,j,l)) ! in mol/m3 ideal gas law
	    boxmols = (boxvol/1e6)*airmols ! in mols

	    ntrtrack = 1
	    ! Gas Phase Chems

            do ntr=1,HowManyEvolveGasChems
		GasConc(ntr)  = tracer(i,j,l,ntr)*airmw/GasMolecularMass(ntr)*1e9! kg/kg_air to ppb

            	if (IsNaN(GasConc(ntr)))then
	    		GasConc(ntr) = 0.0
           	endif

            	if (GasConc(ntr).lt.0)then
	    		GasConc(ntr) = 0.0
           	endif

		ntrtrack = ntrtrack +1
                ntrtrack2 = ntrtrack +1
            enddo

            GasConc(1) = waterconc
!write(*,*) "Before Fire Loop tracers.f90 sum(tracer) L786", sum(tracer) !AD DEBUG
            !write(*,*), "L779 tracers_physics waterconc", waterconc
            !!!! Identifed as a fire grid box so running ASP 
 	    !write(*,*) "L781 tracers_physics rank, i, l, nstep",rank,i,l,nstep !AD DEBUG
!AD DEBUG remove rank criterion as this seems to stay at 0?
               if (i.eq.12.and.l.eq.30.and.nstep.lt.200) then !AD DEBUG
!              if (rank.eq.5.and.i.eq.12.and.l.eq.30.and.nstep.lt.200) then
                  write(*,*) 'Before ASP is ever called...' ,nstep
                  write(*,*) 'water = ',GasConc(1)
                  write(*,*) sum(MassConc(:,1)),MassConc(:,1)
                  write(*,*) 'T = ',temptm
                  write(*,*) 'P = ',prestm
                  write(*,*) 'sza = ',solar_zenith_angle
                  write(*,*) 'bf: API + OH (1381) = ',GasPhaseChemicalRates(1381)
                  write(*,*) 'bf: API + O3 (1402) = ',GasPhaseChemicalRates(1402)
                  write(*,*) 'bf: API + NO3 (1413) = ',GasPhaseChemicalRates(1413)
              end if


	
            !IF (GasConc(2).gt.130.and.nstep.gt.200) THEN
!            write(*,*) 'CO Conc gt 130 chk',GasConc(2)
            IF (GasConc(2).gt.130.and.nstep.gt.cur_dist) THEN	!AD DEBUG
		write(*,*) 'Fire Loop CO Conc gt 130 chk',GasConc(2)
!	     IF (.FALSE.) THEN 	 !AD DEBUG SKIP LOOP ENTIRELY
            	! Aerosol MassConc
            	! Aqueous Phase Chems + Org Phase + Aq Org Phase

		!if (doAerosols.gt.0) then
               do binnum=1,HowManyBins
                  do ntr = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions+HowManyOrgChems+HowManyAqOrgChems

			MassConc(binnum,ntr)=tracer(i,j,l,ntrtrack)*(airdenstm/Av*airmw)*1e12
			ntrtrack = ntrtrack+1
                  enddo
            	  NumConc(binnum) = tracer(i,j,l,1476+binnum)*(airdenstm/Av*airmw*1e-3)

               enddo
	       !write(*,*) "Fire Loop L824 sum(tracer(i,j,l,:))", sum(tracer(i,j,l,:)) !AD DEBUG
               !write(*,*) 'Fire Loop sum(tracer) L824', sum(tracer) !AD DEBUG
               do binnum=1,HowManyBins
                  if (sum(MassConc(binnum,:)).lt.1e-30) then
                     NumConc(binnum) = 0.0
                  end if
               enddo

               ! writing print statements
               t_TotOA = 0.0
!AD DEBUG remove rank criterion as this seems to stay at 0?
               if (i.eq.12.and.l.eq.30) then
!               if (rank.eq.5.and.i.eq.12.and.l.eq.30) then
                  write(*,*) 'Gas+Aerosol going into ASP...' ,nstep
                  write(*,*) 'water = ',GasConc(1)
                  write(*,*) sum(MassConc(:,1)),MassConc(:,1)
                  write(*,*) 'T = ',temptm
                  write(*,*) 'P = ',prestm
                  write(*,*) 'sza = ',solar_zenith_angle
                  write(*,*) 'Mgas = ',Mgas
                  write(*,*) ' SAMASP air dens = ',airdenstm

                  ntrtrack2 = 673


                  do ntr=1,HowManyEvolveGasChems
                     !write(*,*) GasPhaseChemicalNames(ntr),   GasConc(ntr)  , ntr  
                     if (ntr .EQ. findchem('CO', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   GasConc(ntr)  , ntr  
                     endif
                     if (ntr .EQ. findchem('OH', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   GasConc(ntr)  , ntr  
                     endif
                     if (ntr .EQ. findchem('NO', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   GasConc(ntr)  , ntr  
                     endif

                     if (ntr .EQ. findchem('NO3', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   GasConc(ntr)  , ntr  
                     endif
                     if (ntr .EQ. findchem('HO2', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   GasConc(ntr)  , ntr  
                     endif
                     ! HC5SOA
                     if (ntr .EQ. findchem('HC5', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   GasConc(ntr)  , ntr  
                     endif

                     if (ntr .EQ. findchem('HC5SOA', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   GasConc(ntr)  , ntr  
                     endif

                     !HC8SOA
                     if (ntr .EQ. findchem('HC8', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   GasConc(ntr)  , ntr  
                     endif
                     if (ntr .EQ. findchem('HC8SOA', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   GasConc(ntr)  , ntr  
                     endif

                     !OLTSOA
                     if (ntr .EQ. findchem('OLT', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   GasConc(ntr)  , ntr  
                     endif
                     if (ntr .EQ. findchem('OLTSOA', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   GasConc(ntr)  , ntr  
                     endif
                     
                     ! OLISOA
                     if (ntr .EQ. findchem('OLI', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   GasConc(ntr)  , ntr  
                     endif             
                      if (ntr .EQ. findchem('DIEN', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   GasConc(ntr)  , ntr  
                     endif         
                     if (ntr .EQ. findchem('OLISOA', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   GasConc(ntr)  , ntr  
                     endif

                     !TOLSOA
                     if (ntr .EQ. findchem('TOL', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   GasConc(ntr)  , ntr  
                     endif
                     if (ntr .EQ. findchem('TOLSOA', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   GasConc(ntr)  , ntr  
                     endif

                     !XLYSOA
                     if (ntr .EQ. findchem('XYM', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   GasConc(ntr)  , ntr  
                     endif
                     if (ntr .EQ. findchem('XYP', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   GasConc(ntr)  , ntr  
                     endif
                     if (ntr .EQ. findchem('XYO', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   GasConc(ntr)  , ntr  
                     endif
                     if (ntr .EQ. findchem('CSL', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   GasConc(ntr)  , ntr  
                     endif
                     if (ntr .EQ. findchem('XYLSOA', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   GasConc(ntr)  , ntr  
                     endif

                     ! ISOSOA

                     if (ntr .EQ. findchem('ISOP', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   GasConc(ntr)  , ntr  
                     endif
                     if (ntr .EQ. findchem('ISOSOA', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   GasConc(ntr)  , ntr  
                     endif

                     ! TERPSOA
                     if (ntr .EQ. findchem('API', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   GasConc(ntr)  , ntr  
                     endif
                     if (ntr .EQ. findchem('LIM', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   GasConc(ntr)  , ntr  
                     endif
                     if (ntr .EQ. findchem('TERPSOA', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   GasConc(ntr)  , ntr  
                     endif
                  enddo

                  !do ntr = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions+HowManyOrgChems+HowManyAqOrgChems                  
                  do ntr = 56,64

                        write(*,*) tracername(ntrtrack2-1), sum(MassConc(:,ntr)),MassConc(:,ntr)
			ntrtrack2 = ntrtrack2+1
                        t_TotOA = t_TotOA + sum(MassConc(:,ntr))
                  enddo
                  
                  write(*,*) 'NoCo ',sum(NumConc(:)),NumConc(:)
                  write(*,*) 'TotOA = ',t_TotOA
               endif

		do ntr=1,HowManyEvolveGasChems
		    !write(*,*) "Fire Loop L960 i,j,l:",i,j,l,ntr !AD DEBUG
                    !write(*,*) "Fire Loop L960 tracer(i,j,l,ntr):", tracer(i,j,l,ntr) !AD DEBUG
                    !write(*,*) "Fire Loop L960 GasConc(ntr):", GasConc(ntr) !AD DEBUG
            	    if (IsNaN(GasConc(ntr)))then
			write(*,*) 'ATTENTION!!!!', GasPhaseChemicalNames(ntr),tracer(i,j,l,ntr)
		    endif
	        enddo
!write(*,*) "Fire Loop L963 sum(tracer(i,j,l,:))", sum(tracer(i,j,l,:)) !AD DEBUG
!write(*,*) "Fire Loop tracers.f90 sum(tracer) L963", sum(tracer) !AD DEBUG
                !MJA 20190702 Add photolysis rate input

                write(*,*) "Before input PhotolysisRates"
                do irxn=1,51
!                     PhotolysisRates(irxn) = PhotoArr(i,j,kk,irxn)
                      PhotolysisRates(irxn) = PhotoArr(i,j,l,irxn) ! AD DEBUG
                enddo


                !write(*,*) "z level 'l':", l
                !write(*,*) "samvalj(:,l)",samvalj(:,l)
                !write(*,*) "PhotolysisRates going into SAM_wrapper", PhotolysisRates !AD DEBUG
                do irxn=1,51 !AD DEBUG
!                     PhotolysisRates(irxn) = PhotoArr(i,j,kk,irxn) !AD DEBUG
		     PhotolysisRates(irxn) = PhotoArr(i,j,l,irxn) !AD DEBUG above ... shouldn't be kk (not defined till L1180)?	
!                      PhotolysisRates(irxn) = 0 ! AD DEBUG
                enddo !AD DEBUG

                !calculate photo rates not provided by TUV
                PhotolysisRates(49) = 0.33*PhotolysisRates(22)
                PhotolysisRates(50) = 0.16*PhotolysisRates(22)
                PhotolysisRates(51) = 100*PhotolysisRates(25)

		!write(*,*) "PhotolysisRates passed to SAM_wrapper: ", PhotolysisRates !AD DEBUG
		!call flush(6) !AD DEBUG
!write(*,*) "Fire Loop L990 sum(tracer(i,j,l,:))", sum(tracer(i,j,l,:)) !MM DEBUG
!write(*,*) "Fire Loop tracers.f90 sum(tracer) L990", sum(tracer) !MM DEBUG
                write(*,*) "Before SAM_wrapper"
   	        call SAM_wrapper(size(GasConc), &
				size(GasConc), & 
                                solar_zenith_angle, &
                                GasConc, &
                                temptm, &!
				prestm, &   
                                airdens, &
				tstep, &
                                AfterGasConc,  &
                                MassConc, &!
				NumConc, &
				TermVel, &
				Transmissivity,i,l,PhotolysisRates,ExtCoeff, SingScat, Assym)
!write(*,*) "Fire Loop L1005 sum(tracer(i,j,l,:))", sum(tracer(i,j,l,:)) !MM DEBUG
!write(*,*) "Fire Loop tracers.f90 sum(tracer) L1005", sum(tracer) !MM DEBUG
                !write(*,*) "Fire Loop After SAM_wrapper" !AD DEBUG
		!write(*,*) "tracers.f90 NANCHECK sum(NumConc)", sum(NumConc)
                !write(*,*) "tracers.f90 NANCHECK sum(GasConc)", sum(GasConc)
                !write(*,*) "tracers.f90 NANCHECK sum(AfterGasConc)", sum(AfterGasConc)
                !write(*,*) "tracers.f90 NANCHECK sum(MassConc)", sum(MassConc)
            !MJA 20190702 Store Aerosol Optical Properties for grid
                ! start MM DEBUG BLOCK
                write(*,*) "tracers_physics:L1014 After SAM_wrapper" 
                !WRITE(*,*) "shape(ExtCoeff) = ", shape(ExtCoeff)
                !WRITE(*,*) "shape(SingScat) = ", shape(SingScat)
                !WRITE(*,*) "shape(Assym) = ", shape(Assym)
                IF ( ANY( ExtCoeff <= 0) ) then
                   WRITE(*,*) "ExtCoeff has a bad value <= 0"
                   WRITE(*,*) "ExtCoeff = ", ExtCoeff
                END IF
                
                IF ( ANY( SingScat < 0) .OR. ANY(SingScat > 1.0)) then
                   WRITE(*,*) "SingScat has a bad value < 0 OR > 1.0"
                   WRITE(*,*) "SingScat = ", SingScat
                END IF
                
                IF ( ANY( Assym < -1.0) .OR. ANY(Assym > 1.0)) then
                   WRITE(*,*) "Assym has a bad value < -1.0 OR > 1.0"
                   WRITE(*,*) "Assym = ", Assym
                END IF
                ! end MM DEBUG BLOCK
                  
                do nwvn = 1,7
                   ExtCoeffArr(i,j,l,nwvn) = ExtCoeff(nwvn) !mja 20210721
                   SingScatArr(i,j,l,nwvn) = SingScat(nwvn) !mja 20210721
                   AssymArr(i,j,l,nwvn) = Assym(nwvn) !mja 20210721
                enddo
                
                ! start MM DEBUG BLOCK
                write(*,*) "tracers_physics:L1038 After reassignments" 
                !WRITE(*,*) "shape(ExtCoeffArr) = ", shape(ExtCoeffArr)
                !WRITE(*,*) "shape(SingScatArr) = ", shape(SingScatArr)
                !WRITE(*,*) "shape(AssymArr) = ", shape(AssymArr)
                IF ( ANY( ExtCoeffArr <= 0) ) then
                   WRITE(*,*) "ExtCoeffArr has a bad value <= 0"
                   WRITE(*,*) "ExtCoeffArr = ", ExtCoeffArr
                END IF
                
                IF ( ANY( SingScatArr < 0) .OR. ANY(SingScatArr > 1.0)) then
                   WRITE(*,*) "SingScatArr has a bad value < 0 OR > 1.0"
                   WRITE(*,*) "SingScatArr = ", SingScatArr
                END IF
                
                IF ( ANY( AssymArr < -1.0) .OR. ANY(AssymArr > 1.0)) then
                   WRITE(*,*) "AssymArr has a bad value < -1.0 OR > 1.0"
                   WRITE(*,*) "AssymArr = ", AssymArr
                END IF
                ! end MM DEBUG BLOCK
                
                !write(*,*) "Fire Loop After Store Aerosol Optical Properties" !AD DEBUG

		ntrtrack =1
!write(*,*) "Fire Loop L1012 sum(tracer)", sum(tracer)
		! Gas Phase Species
            	do ntr=1,HowManyEvolveGasChems
            	   if (AfterGasConc(ntr).lt.0)then
   	    		AfterGasConc(ntr) = 0.0
           	   endif
            	   if (IsNaN(AfterGasConc(ntr)))then
   	    		AfterGasConc(ntr) = 0.0
           	   endif

		   !write(*,*) "Gas Phase Loop L1022 AfterGasConc(ntr)",ntr, AfterGasConc(ntr) !AD DEBUG
		    
		   tracer(i,j,l,ntr)=AfterGasConc(ntr)*GasMolecularMass(ntr)/airmw/1e9! ppb to kg/kg_air
                   !write(*,*) "Gas Phase Loop L1025 i,j,l:",i,j,l !AD DEBUG
		   !write(*,*) "Gas Phase Loop L1025 tracer(i,j,l,ntr):", tracer(i,j,l,ntr) !AD DEBUG 
		   !write(*,*) "Gas Phase Loop L1027 GasMolecularMass(ntr)",ntr, GasMolecularMass(ntr)
                    
		   ntrtrack  = ntrtrack+1
                   ntrtrack2 = ntrtrack+1
            	enddo

!write(*,*) "Fire Loop tracers.f90 sum(tracer) L1023", sum(tracer) !AD DEBUG
            	! Aqueous Phase Chems + Org Phase + Aq Org Phase
            	do binnum=1,HowManyBins
                   do ntr = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions+HowManyOrgChems+HowManyAqOrgChems
			tracer(i,j,l,ntrtrack)= MassConc(binnum,ntr)/(airdenstm/Av*airmw)/1e12
			ntrtrack = ntrtrack+1
            	   enddo

            	   ! NumConc
		   tracer(i,j,l,1476+binnum)= NumConc(binnum)/(airdenstm/Av*airmw*1e-3)

                enddo ! binnum
               !write(*,*) "Fire Loop tracers.f90 sum(tracer) L1035", sum(tracer) !AD DEBUG
               ! writing print statements
               t_TotOA = 0.0
!AD DEBUG remove rank criterion as this seems to stay at 0?
               if (i.eq.12.and.l.eq.30) then !AD DEBUG
!               if (rank.eq.5.and.i.eq.12.and.l.eq.30) then
                  write(*,*) 'Gas+Aerosol after ASP...' ,nstep
                  write(*,*) 'water = ',AfterGasConc(1)
                  write(*,*) sum(MassConc(:,1)),MassConc(:,1)
                  write(*,*) 'T = ',temptm
                  write(*,*) 'P = ',prestm
                  write(*,*) 'sza = ',solar_zenith_angle
                  write(*,*) 'GasPhaseChemicalRates(NO3 -> NO #7) = ',GasPhaseChemicalRates(7)
                  write(*,*) 'GasPhaseChemicalRates(NO3 -> NO2 #8) = ',GasPhaseChemicalRates(8)
                  ntrtrack2 = 673
                  do ntr=1,HowManyEvolveGasChems
                     !write(*,*) GasPhaseChemicalNames(ntr),   GasConc(ntr)  , ntr  
                     if (ntr .EQ. findchem('CO', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   AfterGasConc(ntr)  , ntr  
                     endif
                     if (ntr .EQ. findchem('OH', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   AfterGasConc(ntr)  , ntr  
                     endif
                     if (ntr .EQ. findchem('NO', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   AfterGasConc(ntr)  , ntr  
                     endif

                     if (ntr .EQ. findchem('NO3', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   AfterGasConc(ntr)  , ntr  
                     endif
                     if (ntr .EQ. findchem('HO2', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   AfterGasConc(ntr)  , ntr  
                     endif
                     ! HC5SOA
                     if (ntr .EQ. findchem('HC5', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   AfterGasConc(ntr)  , ntr  
                     endif

                     if (ntr .EQ. findchem('HC5SOA', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   AfterGasConc(ntr)  , ntr  
                     endif

                     !HC8SOA
                     if (ntr .EQ. findchem('HC8', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   AfterGasConc(ntr)  , ntr  
                     endif
                     if (ntr .EQ. findchem('HC8SOA', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   AfterGasConc(ntr)  , ntr  
                     endif

                     !OLTSOA
                     if (ntr .EQ. findchem('OLT', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   AfterGasConc(ntr)  , ntr  
                     endif
                     if (ntr .EQ. findchem('OLTSOA', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   AfterGasConc(ntr)  , ntr  
                     endif
                     
                     ! OLISOA
                     if (ntr .EQ. findchem('OLI', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   AfterGasConc(ntr)  , ntr  
                     endif             
                      if (ntr .EQ. findchem('DIEN', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   AfterGasConc(ntr)  , ntr  
                     endif         
                     if (ntr .EQ. findchem('OLISOA', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   AfterGasConc(ntr)  , ntr  
                     endif

                     !TOLSOA
                     if (ntr .EQ. findchem('TOL', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   AfterGasConc(ntr)  , ntr  
                     endif
                     if (ntr .EQ. findchem('TOLSOA', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   AfterGasConc(ntr)  , ntr  
                     endif

                     !XLYSOA
                     if (ntr .EQ. findchem('XYM', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   AfterGasConc(ntr)  , ntr  
                     endif
                     if (ntr .EQ. findchem('XYP', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   AfterGasConc(ntr)  , ntr  
                     endif
                     if (ntr .EQ. findchem('XYO', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   AfterGasConc(ntr)  , ntr  
                     endif
                     if (ntr .EQ. findchem('CSL', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   AfterGasConc(ntr)  , ntr  
                     endif
                     if (ntr .EQ. findchem('XYLSOA', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   AfterGasConc(ntr)  , ntr  
                     endif

                     ! ISOSOA

                     if (ntr .EQ. findchem('ISOP', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   AfterGasConc(ntr)  , ntr  
                     endif
                     if (ntr .EQ. findchem('ISOSOA', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   AfterGasConc(ntr)  , ntr  
                     endif

                     ! TERPSOA
                     if (ntr .EQ. findchem('API', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   AfterGasConc(ntr)  , ntr  
                     endif
                     if (ntr .EQ. findchem('LIM', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   AfterGasConc(ntr)  , ntr  
                     endif
                     if (ntr .EQ. findchem('TERPSOA', 0)) then
                        write(*,*) GasPhaseChemicalNames(ntr),   AfterGasConc(ntr)  , ntr  
                     endif

                  enddo

                  !do ntr = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions+HowManyOrgChems+HowManyAqOrgChems
                  do ntr = 56, 64
                       write(*,*) tracername(ntrtrack2-1), sum(MassConc(:,ntr)),MassConc(:,ntr)
                       t_TotOA = t_TotOA + sum(MassConc(:,ntr))
                       ntrtrack2 = ntrtrack2+1
                  enddo

                  write(*,*) 'NoCo ',sum(NumConc(:)),NumConc(:)
                  write(*,*) 'TotOA = ',t_TotOA
               endif
		!write(*,*) "At very end Fire Loop sum(tracer)", sum(tracer) !AD DEBUG
!                IF (sum(tracer) .LT. 0) THEN !AD DEBUG
!			STOP !AD DEBUG
!		ENDIF !AD DEBUG
            !!!! ENDING Identifed as a fire grid box so running ASP 		
            END IF ! GT GasConc(2)

         enddo
      !MJA 20190702 TODO: Call TUV to get photolysis rates for each column AFTER the main z loop is done
         !Call tuv
         !samtstart = hour - 7 !Convert from PST to UTC, kludge that will need to be fixed!
         samtstart = hour -7 !MM DEBUG
         !MM DEBUG: kk and nwvn loops added after samAOD, samSSA, and samASYM
         !dimensions modified to match those in TUV4SAMASP
         DO kk = 1,nzm
           DO nwvn = 1,7
             !samAOD(kk,nwvn) = 1e-6 !MM DEBUG clear sky test
             samAOD(kk,nwvn) = ExtCoeffArr(i,j,kk,nwvn)*dz/1.0e6 !In m   MJA 20220210 DEBUG
             !samAOD(kk,nwvn) = ExtCoeffArr(i,j,kk,nwvn)*dz !In m   MM DEBUG
             !write(*,*) "samAOD",samAOD(:,:) !AD DEBUG
             samSSA(kk,nwvn) = SingScatARR(i,j,kk,nwvn)
             !write(*,*) "samSSA",samSSA(:,:) !AD DEBUG
             samASYM(kk,nwvn) = AssymARR(i,j,kk,nwvn)
             !write(*,*) "samASYM",samASYM(:,:) !AD DEBUG
           ENDDO !nwvn
         ENDDO !kk

         ! start MM DEBUG BLOCK
         write(*,*) "tracers_physics:L1236 After reassignments, before tuv4samasp" 
         !WRITE(*,*) "shape(samAOD) = ", shape(samAOD)
         !WRITE(*,*) "shape(samSSA) = ", shape(samSSA)
         !WRITE(*,*) "shape(samASYM) = ", shape(samASYM)
         IF ( ANY( samAOD <= 0) ) then
            WRITE(*,*) "samAOD has a bad value <= 0"
            WRITE(*,*) "samAOD = ", samAOD
         END IF
                
         IF ( ANY( samSSA < 0) .OR. ANY(samSSA > 1.0)) then
            WRITE(*,*) "samSSA has a bad value < 0 OR > 1.0"
            WRITE(*,*) "samSSA = ", samSSA
         END IF
                
         IF ( ANY( samASYM < -1.0) .OR. ANY(samASYM > 1.0)) then
            WRITE(*,*) "samASYM has a bad value < -1.0 OR > 1.0"
            WRITE(*,*) "samASYM = ", samASYM
         END IF
         ! end MM DEBUG BLOCK
         
         write(*,*) "Before tuv4samasp", samtstart
         call tuv4samasp(samtstart, nzm,samAOD, samSSA, samASYM, samvalj)
         write(*,*) "After tuv4samasp"
         !Store PhotolysisRates for next time step
         do kk=1,nzm
            do irxn=1,51
!               PhotoArr(i,j,kk,irxn) = samvalj(kk,irxn)
	       PhotoArr(i,j,kk,irxn) = samvalj(irxn,kk) !AD DEBUG
               !write(*,*) "kk, irxn********************",kk,irxn !AD DEBUG
               !write(*,*) "samvalj(:,:)",samvalj(irxn,kk) !AD DEBUG
               !write(*,*) "PhotoArr(i,j,:,:)",PhotoArr(i,j,kk,irxn)!AD DEBUG
	       
            enddo !irxn=1,51
         enddo !kk=1,nzm
         write(*,*) "After store PhotolysisRates"
      enddo
   enddo ! end grid loop

write(*,*) "end tracers_physics loop"

!----------------------
 trphys = 0. ! Default tendency due to physics. You code should compute this to output statistics.
 end subroutine tracers_physics
!----------------------------------

 subroutine tracers_hbuf_init(namelist,deflist,unitlist,status,average_type,count,trcount)

! Initialize the list of tracers statistics variables written in statistics.f90

   character(*) namelist(*), deflist(*), unitlist(*)
   integer status(*),average_type(*),count,trcount
   integer ntr
!write(*,*) "L1169 tracers.f90"

   do ntr=1,ntracers
     count = count + 1
     trcount = trcount + 1
     !namelist(count) = trim(tracername(ntr))
     namelist(count) = tracername(ntr)
     !deflist(count) = trim(tracername(ntr))
     deflist(count) = tracername(ntr)
     unitlist(count) = trim(tracerunits(ntr))
     status(count) = 1
     average_type(count) = 0		
     count = count + 1
     trcount = trcount + 1
     namelist(count) = trim(tracername(ntr))//'FLX'
     deflist(count) = 'Total flux of '//trim(tracername(ntr))
     unitlist(count) = trim(tracerunits(ntr))//' kg/m2/s'
     status(count) = 1
     average_type(count) = 0
     count = count + 1
     trcount = trcount + 1
     namelist(count) = trim(tracername(ntr))//'FLXS'
     deflist(count) = 'SGS flux of '//trim(tracername(ntr))
     unitlist(count) = trim(tracerunits(ntr))//' kg/m2/s'
     status(count) = 1
     average_type(count) = 0
     count = count + 1
     trcount = trcount + 1
     namelist(count) = trim(tracername(ntr))//'ADV'
     deflist(count) = 'Tendency of '//trim(tracername(ntr)//'due to vertical advection')
     unitlist(count) = trim(tracerunits(ntr))//'/day'
     status(count) = 1
     average_type(count) = 0
     count = count + 1
     trcount = trcount + 1
     namelist(count) = trim(tracername(ntr))//'DIFF'
     deflist(count) = 'Tendency of '//trim(tracername(ntr)//'due to vertical SGS transport')
     unitlist(count) = trim(tracername(ntr))//'/day'
     status(count) = 1
     average_type(count) = 0
     count = count + 1
     trcount = trcount + 1
     namelist(count) = trim(tracername(ntr))//'PHYS'
     deflist(count) = 'Tendency of '//trim(tracername(ntr)//'due to physics')
     unitlist(count) = trim(tracername(ntr))//'/day'
     status(count) = 1
     average_type(count) = 0
   enddo

 end subroutine tracers_hbuf_init

end module tracers
