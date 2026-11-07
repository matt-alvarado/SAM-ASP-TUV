module samasp_photo

! Photolysis rates passed from TUV to ASP.
!
! ASP's photolysis reactions (GasChemicalMechanism.in, "; 0.; N.") use rate
! column N of its photolysis-rate table (PhotoRates.in; HowManyPhotoRxns
! columns). TUV computes J-values for the reactions selected in its input
! file, identified by their labels (jlabel). PhotoMap.in, in the ASP
! input-deck directory, says which TUV reaction(s) feed each ASP column:
!
!     ASP column ; scale factor ; TUV reaction label
!
! Lines with the same column are summed. The run stops at start-up if a
! column is outside 1..HowManyPhotoRxns, if a column used by the mechanism
! has no line, or (at the first TUV call) if a label is not among TUV's
! selected reactions.

implicit none
private

public :: samasp_photo_init, samasp_photo_set_labels, samasp_photo_map, &
          photo_col_of_label

integer, public, protected :: nPhoto = 0       ! ASP photolysis rate columns
logical, public, protected :: photo_labels_set = .false.

integer :: nMap = 0
integer, allocatable :: map_col(:)              ! ASP column
integer, allocatable :: map_tuv(:)              ! index among TUV's selected reactions
real, allocatable :: map_scale(:)               ! default real, as the constants it replaces
character(len=50), allocatable :: map_label(:)  ! TUV label

contains

subroutine samasp_photo_init(deckdir, masterproc)
  use Chemistry, only: HowManyPhotoRxns, GasPhaseReactionRates, &
       GasPhaseReactionRateArraySize
  character(len=*), intent(in) :: deckdir
  logical, intent(in) :: masterproc

  integer, parameter :: maxmap = 500
  integer :: col(maxmap), u, ios, p1, p2, L, n, k
  real :: scl(maxmap)
  character(len=50) :: lab(maxmap)
  character(len=512) :: line, fname
  logical, allocatable :: used(:)

  nPhoto = HowManyPhotoRxns
  fname = trim(deckdir)//'PhotoMap.in'
  open(newunit=u, file=trim(fname), status='old', action='read', iostat=ios)
  if (ios /= 0) call photo_error('cannot open '//trim(fname))

  n = 0
  do
     read(u, '(a)', iostat=ios) line
     if (ios /= 0) exit
     line = adjustl(line)
     if (len_trim(line) == 0 .or. line(1:1) == '!') cycle
     p1 = index(line, ';')
     p2 = p1 + index(line(p1+1:), ';')
     if (p1 == 0 .or. p2 == p1) call photo_error('PhotoMap.in: expected '// &
          '"column ; scale ; TUV label" in line: '//trim(line))
     n = n + 1
     if (n > maxmap) call photo_error('PhotoMap.in: more than 500 lines')
     read(line(1:p1-1), *, iostat=ios) col(n)
     if (ios /= 0) call photo_error('PhotoMap.in: bad column in line: '//trim(line))
     read(line(p1+1:p2-1), *, iostat=ios) scl(n)
     if (ios /= 0) call photo_error('PhotoMap.in: bad scale in line: '//trim(line))
     lab(n) = adjustl(line(p2+1:))
     if (col(n) < 1 .or. col(n) > nPhoto) then
        write(line, '(a,i0,a,i0,a)') 'PhotoMap.in: column ', col(n), &
             ' is outside 1..', nPhoto, ' (the columns of PhotoRates.in)'
        call photo_error(trim(line))
     end if
  end do
  close(u)

  nMap = n
  allocate(map_col(n), map_tuv(n), map_scale(n), map_label(n))
  map_col = col(1:n); map_scale = scl(1:n); map_label = lab(1:n); map_tuv = 0

  ! every column used by a photolysis reaction of the mechanism needs a source
  allocate(used(nPhoto))
  used = .false.
  do L = 1, GasPhaseReactionRateArraySize(1)
     if (int(GasPhaseReactionRates(L,1)) == 1 .and. &
         int(GasPhaseReactionRates(L,3)) == 0) then
        k = int(GasPhaseReactionRates(L,4))
        if (k < 1 .or. k > nPhoto) then
           write(line, '(a,i0,a,i0)') 'the mechanism uses photolysis column ', &
                k, ' but PhotoRates.in has columns 1..', nPhoto
           call photo_error(trim(line))
        end if
        used(k) = .true.
     end if
  end do
  do k = 1, nPhoto
     if (used(k) .and. .not. any(map_col == k)) then
        write(line, '(a,i0,a)') 'photolysis column ', k, &
             ' is used by the mechanism but has no line in PhotoMap.in'
        call photo_error(trim(line))
     end if
  end do

  if (masterproc) write(*,'(a,i0,a,i0,a)') 'SAM-ASP photolysis: ', nMap, &
       ' PhotoMap.in entries for ', nPhoto, ' ASP rate columns'
end subroutine samasp_photo_init

! Match the PhotoMap.in labels to TUV's selected reactions (first TUV call).
subroutine samasp_photo_set_labels(nj, jlabel)
  integer, intent(in) :: nj
  character(len=*), intent(in) :: jlabel(:)
  integer :: k, ij
  do k = 1, nMap
     map_tuv(k) = 0
     do ij = 1, nj
        if (trim(adjustl(jlabel(ij))) == trim(map_label(k))) then
           map_tuv(k) = ij
           exit
        end if
     end do
     if (map_tuv(k) == 0) call photo_error('PhotoMap.in: TUV reaction "'// &
          trim(map_label(k))//'" is not among the reactions selected in '// &
          'the TUV input file')
  end do
  photo_labels_set = .true.
end subroutine samasp_photo_set_labels

! jasp(iz, col) = sum over PhotoMap.in lines for col of scale * jtuv(ij, iz),
! for iz = 1..nz. Contributions are added in the order of TUV's reactions.
subroutine samasp_photo_map(nj, nz, jtuv, jasp)
  integer, intent(in) :: nj, nz
  real*8, intent(in) :: jtuv(:,:)          ! (TUV reaction, level)
  real*8, intent(out) :: jasp(:,:)         ! (level, ASP column)
  integer :: ij, k
  jasp(1:nz,:) = 0.0d0
  do ij = 1, nj
     do k = 1, nMap
        if (map_tuv(k) /= ij) cycle
        if (map_scale(k) == 1.0) then
           jasp(1:nz, map_col(k)) = jasp(1:nz, map_col(k)) + jtuv(ij, 1:nz)
        else
           jasp(1:nz, map_col(k)) = jasp(1:nz, map_col(k)) + &
                map_scale(k)*jtuv(ij, 1:nz)
        end if
     end do
  end do
end subroutine samasp_photo_map

! ASP column fed by TUV reaction `label` with scale 1 (0 if none), e.g. for
! choosing J-value output fields by reaction.
integer function photo_col_of_label(label)
  character(len=*), intent(in) :: label
  integer :: k
  photo_col_of_label = 0
  do k = 1, nMap
     if (trim(map_label(k)) == trim(label) .and. map_scale(k) == 1.0) then
        photo_col_of_label = map_col(k)
        return
     end if
  end do
end function photo_col_of_label

subroutine photo_error(msg)
  character(len=*), intent(in) :: msg
  write(*,'(a)') 'ERROR (samasp_photo): '//trim(msg)
  call task_abort()
end subroutine photo_error

end module samasp_photo
