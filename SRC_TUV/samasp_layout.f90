module samasp_layout

! Layout of the SAM tracer array for the ASP chemistry and aerosol species.
!
! All counts come from ASP after InitializeASP has read the input decks, so
! a change to the ASP species lists or the number of size bins does not need
! edits to SAM code. Only ntracers (domain.f90, a compile-time constant) has
! to match; samasp_layout_init stops with the value to use if it does not.
!
! tracer(:,:,:,n) holds, for n =
!   itr_gas(j)       = j                        gas species j, j = 1..nGas
!   itr_mass(bin,j)  = iAero0 + (bin-1)*nSpec + j
!                                               aerosol species j in size bin
!                                               bin, with j =
!                      1 .. nAq                 aqueous species and ions
!                      nAq+1 .. nAq+nOrg        organic (hydrophobic) phase
!                      nAq+nOrg+1 .. nSpec      aqueous organic species
!   itr_num(bin)     = iNum0 + bin              particle number in bin
! and nTracersNeeded = iNum0 + nBins.

implicit none
private

public :: samasp_layout_init, itr_gas, itr_mass, itr_num, &
          aero_species_name, num_conc_name, samasp_gas_index

integer, public, protected :: nGas = 0     ! gas species (ASP evolving gases)
integer, public, protected :: nAq = 0      ! aqueous species + cations + anions
integer, public, protected :: nOrg = 0     ! organic-phase species
integer, public, protected :: nAqOrg = 0   ! aqueous organic species
integer, public, protected :: nSpec = 0    ! aerosol species per bin
integer, public, protected :: nBins = 0    ! aerosol size bins
integer, public, protected :: iAero0 = 0   ! offset of the aerosol block
integer, public, protected :: iNum0 = 0    ! offset of the number block
integer, public, protected :: nTracersNeeded = 0
integer, public, protected :: iH2Og = 0    ! gas index of H2O
integer, public, protected :: iCOg = 0     ! gas index of CO

contains

subroutine samasp_layout_init(ntracers, masterproc)
  use ModelParameters, only: HowManyBins
  use Chemistry, only: HowManyEvolveGasChems, HowManyAqChems, &
       HowManyAqCations, HowManyAqAnions, HowManyOrgChems, HowManyAqOrgChems
  integer, intent(in) :: ntracers
  logical, intent(in) :: masterproc

  nGas   = HowManyEvolveGasChems
  nAq    = HowManyAqChems + HowManyAqCations + HowManyAqAnions
  nOrg   = HowManyOrgChems
  nAqOrg = HowManyAqOrgChems
  nSpec  = nAq + nOrg + nAqOrg
  nBins  = HowManyBins

  iAero0 = nGas
  iNum0  = nGas + nBins*nSpec
  nTracersNeeded = iNum0 + nBins

  iH2Og = samasp_gas_index('H2O')
  iCOg  = samasp_gas_index('CO')

  if (masterproc) then
     write(*,'(a)') 'SAM-ASP tracer layout:'
     write(*,'(a,i6)') '  gas species              ', nGas
     write(*,'(a,i6)') '  aerosol size bins        ', nBins
     write(*,'(a,i6,a,3i5,a)') '  aerosol species per bin  ', nSpec, &
          ' (aq, org, aq-org =', nAq, nOrg, nAqOrg, ')'
     write(*,'(a,i6)') '  tracers needed           ', nTracersNeeded
     write(*,'(a,2i6)') '  gas index of H2O, CO     ', iH2Og, iCOg
  end if

  if (nBins > 99) call layout_error('more than 99 aerosol bins; tracer '// &
       'names use a 2-digit bin suffix')
  if (nTracersNeeded /= ntracers) then
     if (masterproc) write(*,'(a,i6,a,i6,a)') 'ERROR: the ASP input decks need ', &
          nTracersNeeded, ' tracers but ntracers = ', ntracers, &
          '. Set ntracers to that value in domain.f90 and rebuild.'
     call task_abort()
  end if
end subroutine samasp_layout_init

integer function itr_gas(j)
  integer, intent(in) :: j
  itr_gas = j
end function itr_gas

integer function itr_mass(bin, j)
  integer, intent(in) :: bin, j
  itr_mass = iAero0 + (bin-1)*nSpec + j
end function itr_mass

integer function itr_num(bin)
  integer, intent(in) :: bin
  itr_num = iNum0 + bin
end function itr_num

! Tracer name of aerosol species j in size bin `bin`, e.g. 'K+01', 'POA105',
! 'AqPOA110' (species name plus a 2-digit bin number).
character(len=20) function aero_species_name(bin, j)
  use Chemistry, only: AqPhaseChemicalNames, OrgPhaseChemicalNames, &
       AqOrgPhaseChemicalNames
  integer, intent(in) :: bin, j
  character(len=2) :: suffix
  write(suffix, '(i2.2)') bin
  if (j <= nAq) then
     aero_species_name = trim(AqPhaseChemicalNames(j))//suffix
  else if (j <= nAq + nOrg) then
     aero_species_name = trim(OrgPhaseChemicalNames(j-nAq))//suffix
  else
     aero_species_name = 'Aq'//trim(AqOrgPhaseChemicalNames(j-nAq-nOrg))//suffix
  end if
end function aero_species_name

character(len=20) function num_conc_name(bin)
  integer, intent(in) :: bin
  write(num_conc_name, '(a,i2.2)') 'NoCo', bin
end function num_conc_name

! Index of a gas species in ASP's evolving-gas list (= its tracer index);
! stops if the species is not in the mechanism.
integer function samasp_gas_index(name)
  use Chemistry, only: GasPhaseChemicalNames
  character(len=*), intent(in) :: name
  integer :: i
  samasp_gas_index = 0
  do i = 1, nGas
     if (trim(GasPhaseChemicalNames(i)) == trim(name)) then
        samasp_gas_index = i
        return
     end if
  end do
  call layout_error('gas species '//trim(name)//' is not among the '// &
       'evolving gases of the ASP mechanism')
end function samasp_gas_index

subroutine layout_error(msg)
  character(len=*), intent(in) :: msg
  write(*,'(a)') 'ERROR (samasp_layout): '//trim(msg)
  call task_abort()
end subroutine layout_error

end module samasp_layout
