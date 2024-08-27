program crm

!       Main module.

use vars
use hbuffer
use microphysics
use sgs
use tracers
use movies, only: init_movies

implicit none
character *256 cmd
integer k, icyc, nn, nstatsteps
double precision cputime, oldtime, init_time, elapsed_time !bloss wallclocktime
double precision usrtime, systime
integer itmp1(1), oldstep

!-------------------------------------------------------------------
! determine the rank of the current task and of the neighbour's ranks

call task_init()
!------------------------------------------------------------------
! print time, version, etc

if(masterproc) call header()	
!------------------------------------------------------------------
! Initialize timing library.  2nd arg 0 means disable, 1 means enable

   call t_setoptionf (1, 0)
   call t_initializef ()

   call t_startf ('total')
   call t_startf ('initialize')
!------------------------------------------------------------------
! Get initial time of job

call t_stampf(init_time,usrtime,systime)
!------------------------------------------------------------------

call init()     ! initialize some statistics arrays
call setparm()	! set all parameters and constants

!------------------------------------------------------------------
! Initialize or restart from the save-dataset:

if(nrestart.eq.0) then
   day=day0 
   call setgrid() ! initialize vertical grid structure
   call setdata() ! initialize all variables
elseif(nrestart.eq.1) then
   call read_all()
   call setgrid() ! initialize vertical grid structure
   call diagnose()
   call sgs_init()
   call micro_init()  !initialize microphysics
elseif(nrestart.eq.2) then  ! branch run
   call read_all()
   call setgrid() ! initialize vertical grid structure
   call diagnose()
   call setparm() ! overwrite the parameters
   call sgs_init()
   call micro_init()  !initialize microphysics
   nstep = 0
   day0 = day
else
   print *,'Error: confused by value of NRESTART'
   call task_abort() 
endif

write(*,*) 'In main.f90 varying parameters: '
write(*,*) 'caseid = ',caseid
write(*,*) 'day0 =', day0
write(*,*) 'latitude0 = ',latitude0
write(*,*) 'longitude0 = ',longitude0


call init_movies()
call stat_2Dinit(1) ! argument of 1 means storage terms in stats are reset
call tracers_init() ! initialize tracers


call setforcing()



if(masterproc) call printout()
!------------------------------------------------------------------
!  Initialize statistics buffer:

call hbuf_init()

!------------------------------------------------------------------
nstatis = nstat/nstatfrq
nstat = nstatis * nstatfrq
nstatsteps = 0
call t_stopf ('initialize')
!------------------------------------------------------------------
!   Main time loop    
!------------------------------------------------------------------
call t_stampf(oldtime,usrtime,systime)
oldstep = nstep


do while(nstep.lt.nstop.and.nelapse.gt.0) 
        
  nstep = nstep + 1
  !write(*,*) "In main.f90..... nstep = ",nstep
  time = time + dt
  day = day0 + nstep*dt/86400.
  nelapse = nelapse - 1
!------------------------------------------------------------------
!  Check if the dynamical time step should be decreased 
!  to handle the cases when the flow being locally linearly unstable
!------------------------------------------------------------------

  ncycle = 1

  call kurant()

  total_water_before = total_water()
  total_water_evap = 0.
  total_water_prec = 0.
  total_water_ls = 0.

  do icyc=1,ncycle
     icycle = icyc
     dtn = dt/ncycle
     dt3(na) = dtn
     dtfactor = dtn/dt

     if(mod(nstep,nstatis).eq.0.and.icycle.eq.ncycle) then
        nstatsteps = nstatsteps + 1
        dostatis = .true.
        if(masterproc) print *,'Collecting statistics...'
     else
        dostatis = .false.
     endif
     !bloss:make special statistics flag for radiation,since it's only updated at icycle==1.
     dostatisrad = .false.
     if(mod(nstep,nstatis).eq.0.and.icycle.eq.1) dostatisrad = .true.

!---------------------------------------------
!  	the Adams-Bashforth scheme in time

     call abcoefs()
 
!---------------------------------------------
!  	initialize stuff: 
	
     call zero()

!-----------------------------------------------------------
!       Buoyancy term:
	     
     call buoyancy()

!------------------------------------------------------------

     total_water_ls =  total_water_ls - total_water()

!------------------------------------------------------------
!       Large-scale and surface forcing:

     call forcing()

!----------------------------------------------------------
!       Nadging to sounding:

     call nudging()

!----------------------------------------------------------
!   	spange-layer damping near the upper boundary:

     if(dodamping) call damping()

!----------------------------------------------------------

     total_water_ls =  total_water_ls + total_water()

!---------------------------------------------------------
!   Ice fall-out
   
      if(docloud) then
          call ice_fall()
      end if

!----------------------------------------------------------
!     Update scalar boundaries after large-scale processes:

     call boundaries(3)

!---------------------------------------------------------
!     Update boundaries for velocities:

      call boundaries(0)

!-----------------------------------------------
!     surface fluxes:

     if(dosurface) call surface()

!-----------------------------------------------------------
!  SGS physics:

     if (dosgs) call sgs_proc()

!----------------------------------------------------------
!     Fill boundaries for SGS diagnostic fields:

     call boundaries(4)
!-----------------------------------------------
!       advection of momentum:

     call advect_mom()

!----------------------------------------------------------
!	SGS effects on momentum:

     if(dosgs) call sgs_mom()

!-----------------------------------------------------------
!       Coriolis force:
	     
     if(docoriolis) call coriolis()
	 
!---------------------------------------------------------
!       compute rhs of the Poisson equation and solve it for pressure. 

     call pressure()

!---------------------------------------------------------
!       find velocity field at n+1/2 timestep needed for advection of scalars:
!  Note that at the end of the call, the velocities are in nondimensional form.
	 
     call adams()

!----------------------------------------------------------
!     Update boundaries for all prognostic scalar fields for advection:

     call boundaries(2)

!---------------------------------------------------------
!      advection of scalars :

     call advect_all_scalars()

!-----------------------------------------------------------
!    Convert velocity back from nondimensional form:

      call uvw()

!----------------------------------------------------------
!     Update boundaries for scalars to prepare for SGS effects:

     call boundaries(3)
   
!---------------------------------------------------------
!      SGS effects on scalars :

     if (dosgs) call sgs_scalars()

!-----------------------------------------------------------
!       Handle upper boundary for scalars

     if(doupperbound) call upperbound()

!-----------------------------------------------------------
!       Cloud condensation/evaporation and precipitation processes:

      if(docloud.or.dosmoke) call micro_proc()

!----------------------------------------------------------
!  Tracers' physics:

      call tracers_physics()
      !write(*,*) "IN MAIN.F90 END TRACERS_PHYSICS CALL" !AD DEBUG
!-----------------------------------------------------------
!	Radiation

      if(dolongwave.or.doshortwave) then 
        !write(*,*) "IN MAIN.F90 radiation" !AD DEBUG
	call radiation()     
      end if

!-----------------------------------------------------------
!    Compute diagnostic fields:
      !write(*,*) "IN MAIN.F90 before diagnose" !AD DEBUG
      call diagnose()

!----------------------------------------------------------

! Rotate the dynamic tendency arrays for Adams-bashforth scheme:
      write(*,*) "IN MAIN.F90 L294"
      nn=na
      na=nc
      nc=nb
      nb=nn
   end do ! icycle	
      write(*,*) "IN MAIN.F90 L300"
  total_water_after = total_water()
!----------------------------------------------------------
!  collect statistics, write save-file, etc.
write(*,*) "IN MAIN.F90 L304"
   call stepout(nstatsteps)
write(*,*) "IN MAIN.F90 L306 after stepout"  
!----------------------------------------------------------
!----------------------------------------------------------

   ! bloss:  only stop when stat and 2D fields are output.
   if ( (nstep.lt.nstop).and.(mod(nstep,nstat).eq.0) .AND. &
        ( (.NOT.save2Davg) .OR. (mod(nstep,nsave2D).eq.0) ) ) then

     if(masterproc) then
       call t_stampf(cputime,usrtime,systime)
       write(*,999) cputime-oldtime, float(nstep-oldstep)*dt/3600.
999    format('CPU TIME = ',f12.4, ' OVER ', f8.2,' MODEL HOURS')

       elapsed_time = cputime-init_time
       if ((elapsed_time+2.*(cputime-oldtime)).gt.60.*float(nelapsemin)) then
         nelapse=0 ! Job will stop when nelapse=0
         write(*,*) 'Job stopping -- nelapsemin reached'
       end if

       oldtime = cputime
       oldstep = nstep
     end if

      if(dompi) then
        if(masterproc) itmp1 = nelapse
        call task_bcast_integer(0,itmp1,1)
        nelapse = itmp1(1)
      end if
   end if

end do ! main loop

!----------------------------------------------------------
!----------------------------------------------------------

   call t_stopf('total')
   if(masterproc) call t_prf(rank)

if(masterproc) write(*,*) 'Finished with SAM, exiting...'
if(masterproc) then
	cmd ='mv ./OUT_3D/'//trim(case)//'_'//trim(case)//'_12.bin2D COMPLETE_RUNS_AEROSOL/'//trim(case)//'_aerosol_12.bin2D'
	write(*,*) cmd
	call system(cmd)
end if
call task_stop()

end program crm
