module params

use grid, only: nzm

implicit none

!   Constants:

real, parameter :: cp = 1004.             ! Specific heat of air, J/kg/K
real, parameter :: ggr = 9.81             ! Gravity acceleration, m/s2
real, parameter :: lcond = 2.5104e+06     ! Latent heat of condensation, J/kg
real, parameter :: lfus = 0.3336e+06      ! Latent heat of fusion, J/kg
real, parameter :: lsub = 2.8440e+06      ! Latent heat of sublimation, J/kg
real, parameter :: rv = 461.              ! Gas constant for water vapor, J/kg/K
real, parameter :: rgas = 287.            ! Gas constant for dry air, J/kg/K
real, parameter :: diffelq = 2.21e-05     ! Diffusivity of water vapor, m2/s
real, parameter :: therco = 2.40e-02      ! Thermal conductivity of air, J/m/s/K
real, parameter :: muelq = 1.717e-05      ! Dynamic viscosity of air

real, parameter :: fac_cond = lcond/cp 
real, parameter :: fac_fus = lfus/cp
real, parameter :: fac_sub = lsub/cp

real, parameter ::  pi = 3.141592653589793

!
! internally set parameters:

real   epsv     ! = (1-eps)/eps, where eps= Rv/Ra, or =0. if dosmoke=.true.
logical:: dosubsidence = .false.
real fcorz      ! Vertical Coriolis parameter
real coszrs
real salt_factor ! correction factor for water vapor saturation over sea-water

integer:: ncycle_max = 4  ! maximum number of subcycling within dt

!----------------------------------------------
! Parameters set by PARAMETERS namelist:
! Initialized to default values.
!----------------------------------------------
real:: plume_y = 0. ! width of plume in y direction
real:: plume_x = 0. ! width of plume in x direction
real:: ug = 0.        ! Velocity of the Domain's drift in x direction
real:: vg	= 0.  ! Velocity of the Domain's drift in y direction
real:: fcor = -999.   ! Coriolis parameter	
real:: longitude0 = 0.    ! latitude of the domain's center 
real:: latitude0  = 0.    ! longitude of the domain's center 
real:: nxco2 = 1         ! factor to modify co2 concentration
logical:: doradlat = .false.
logical:: doradlon = .false.
real:: massemis = 0. ! mass emission flux of biomass burning - CRL 
! Fire emission timing and location (tracers_physics). Emissions are added
! on steps fire_nstep_start+1 .. fire_nstep_start + 2*(plume_y-dy)/vg/dt,
! in the box of GLOBAL x indices fire_i1..fire_i2 and levels fire_k1..fire_k2
! (at y index 1). ASP is called after that, where CO exceeds 130 ppb.
integer:: fire_nstep_start = 1080
integer:: fire_i1 = 100, fire_i2 = 115
integer:: fire_k1 = 32, fire_k2 = 55
! ASP input-deck and output directories, relative to the run directory
! (each ends in '/'); passed to InitializeASP.
character(256):: asp_deck_dir = 'SRC_TUV/ASP/InorganicInputDecks/'
! Calendar year of the run; day0 and day are days since 1 January 00 UTC
! of that year. Used for TUV's date (Earth-Sun distance, solar declination)
integer:: year0 = 2009
! ASP (chemistry, aerosol) is run only where CO exceeds this (ppb), i.e. in
! the plume, after the fire emissions have stopped
real:: asp_co_threshold_ppb = 130.
character(256):: asp_output_dir = 'SRC_TUV/ASP/OutputFiles/'

real(8):: tabs_s =0.	! surface temperature,K
real:: delta_sst = 0.   ! amplitude of sin-pattern of sst about tabs_s (ocean_type=1)
real:: depth_slab_ocean = 2. ! thickness of the slab-ocean (m)
real:: Szero = 0.  ! mean ocean transport (W/m2)
real:: deltaS = 0. ! amplitude of linear variation of ocean transport (W/m2)
real:: timesimpleocean = 0. ! time to start simple ocean

real::   fluxt0 =0.  ! surface sensible flux, Km/s
real::   fluxq0 =0.  ! surface latent flux, m/s
real::   tau0   =0.  ! surface stress, m2/s2
real::   z0     =0.035	! roughness length
real::   soil_wetness =0.4! wetness coeff for soil (from 0 to 1.)
integer:: ocean_type =0 ! type of SST forcing
logical:: cem =.false.    ! flag for Cloud Ensemble Model
logical:: les =.false.    ! flag for Large-Eddy Simulation
logical:: ocean =.false.  ! flag indicating that surface is water
logical:: land =.false.   ! flag indicating that surface is land
logical:: sfc_flx_fxd =.false. ! surface sensible flux is fixed
logical:: sfc_tau_fxd =.false.! surface drag is fixed

real:: timelargescale =0. ! time to start large-scale forcing

! nudging boundaries (between z1 and z2, where z2 > z1): 
real:: nudging_uv_z1 =-1., nudging_uv_z2 = 1000000.
real:: nudging_t_z1 =-1., nudging_t_z2 = 1000000.
real:: nudging_q_z1 =-1., nudging_q_z2 = 1000000.
real:: tauls = 99999999.    ! nudging-to-large-scaler-profile time-scale
real:: tautqls = 99999999.! nudging-to-large-scaler-profile time-scale for scalars

logical:: dodamping = .false.
logical:: doupperbound = .false. 
logical:: docloud = .false. 
logical:: doprecip = .false.
logical:: dolongwave = .false. 
logical:: doshortwave = .false.
logical:: dosgs = .false.
logical:: docoriolis = .false. 
logical:: docoriolisz = .false. 
logical:: dofplane = .true.
logical:: dosurface = .false. 
logical:: dolargescale = .false. 
logical:: doradforcing = .false.
logical:: dosfcforcing = .false. 
logical:: doradsimple = .false. 
logical:: donudging_uv = .false. 
logical:: donudging_tq = .false.
logical:: donudging_t = .false. 
logical:: donudging_q = .false.
logical:: doensemble = .false. 
logical:: dowallx = .false. 
logical:: dowally = .false. 
logical:: docolumn = .false. 
logical:: docup = .false.
logical:: doperpetual = .false. 
logical:: doseasons = .false. 
logical:: doradhomo = .false. 
logical:: dosfchomo = .false.
logical:: dossthomo = .false. 
logical:: dodynamicocean = .false. 
logical:: dosolarconstant = .false.
logical:: dotracers = .false. 
logical:: dosmoke = .false. 
logical:: notracegases = .false.
logical:: doseawater = .false.

! Specify solar constant and zenith angle for perpetual insolation.
! Based onn Tompkins and Graig (1998)
! Note that if doperpetual=.true. and dosolarconstant=.false.
! the insolation will be set to the daily-averaged value on day0.
real:: solar_constant = 685. ! solar constant (in W/m2)
real:: zenith_angle = 51.7   ! zenith angle (in degrees)

integer:: nensemble =0   ! the number of subensemble set of perturbations
integer:: perturb_type  = 0 ! type of initial noise in setperturb()

! Initial bubble parameters. Activated when perturb_type = 2
  real:: bubble_x0 = 0.
  real:: bubble_y0 = 0.
  real:: bubble_z0 = 0.
  real:: bubble_radius_hor = 0.
  real:: bubble_radius_ver = 0.
  real:: bubble_dtemp = 0.
  real:: bubble_dq = 0.


end module params
