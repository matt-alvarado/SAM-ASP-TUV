module samasp_balance

! Load balancing of per-cell work (ASP chemistry/aerosol calls) across MPI
! ranks. Each rank has ncell work items (columns of bufin, nin values each).
! The items of all ranks, taken in global order (rank 0's first, then rank
! 1's, ...), are split into nproc contiguous blocks of nearly equal size;
! rank r computes block r with the routine `work`, and the nout results per
! item are returned to the rank the item came from, in its original order.
!
! The results are those of computing every item on its own rank only if the
! work for an item depends on that item's inputs alone (for ASP this is
! ensured by the cell-state reset, SetASPCellReset).

implicit none
private
public :: balanced_compute

abstract interface
   subroutine work_routine(nin, nout, n, in, out)
     integer, intent(in) :: nin, nout, n
     real*8, intent(in) :: in(nin, n)
     real*8, intent(out) :: out(nout, n)
   end subroutine work_routine
end interface

contains

subroutine balanced_compute(dompi, nin, nout, ncell, bufin, bufout, work)
  include 'mpif.h'
  logical, intent(in) :: dompi
  integer, intent(in) :: nin, nout, ncell
  real*8, intent(in) :: bufin(nin, ncell)
  real*8, intent(out) :: bufout(nout, ncell)
  procedure(work_routine) :: work

  integer :: nproc, me, ierr, r, total, myoff, mylo, myhi, nmine, a, b
  integer, allocatable :: counts(:), offs(:), lo(:), hi(:)
  integer, allocatable :: scnt(:), sdsp(:), rcnt(:), rdsp(:)
  real*8, allocatable :: myin(:,:), myout(:,:)

  if (.not. dompi) then
     if (ncell > 0) call work(nin, nout, ncell, bufin, bufout)
     return
  end if

  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, me, ierr)
  allocate(counts(0:nproc-1), offs(0:nproc-1), lo(0:nproc-1), hi(0:nproc-1), &
           scnt(0:nproc-1), sdsp(0:nproc-1), rcnt(0:nproc-1), rdsp(0:nproc-1))
  call MPI_Allgather(ncell, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

  total = sum(counts)
  if (total == 0) return
  offs(0) = 0
  do r = 1, nproc-1
     offs(r) = offs(r-1) + counts(r-1)
  end do
  ! block of global items [lo(r), hi(r)) computed by rank r
  do r = 0, nproc-1
     lo(r) = int((int(total,8)*r)/nproc)
     hi(r) = int((int(total,8)*(r+1))/nproc)
  end do
  myoff = offs(me); mylo = lo(me); myhi = hi(me)
  nmine = myhi - mylo

  ! my items [myoff, myoff+ncell) go to the ranks whose blocks overlap them;
  ! I receive from each rank the part of its items inside my block
  do r = 0, nproc-1
     a = max(myoff, lo(r)); b = min(myoff+ncell, hi(r))
     scnt(r) = max(0, b-a)
     a = max(offs(r), mylo); b = min(offs(r)+counts(r), myhi)
     rcnt(r) = max(0, b-a)
  end do
  sdsp(0) = 0; rdsp(0) = 0
  do r = 1, nproc-1
     sdsp(r) = sdsp(r-1) + scnt(r-1)
     rdsp(r) = rdsp(r-1) + rcnt(r-1)
  end do

  allocate(myin(nin, max(nmine,1)), myout(nout, max(nmine,1)))
  call MPI_Alltoallv(bufin, scnt*nin, sdsp*nin, MPI_DOUBLE_PRECISION, &
                     myin, rcnt*nin, rdsp*nin, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
  if (nmine > 0) call work(nin, nout, nmine, myin, myout)
  call MPI_Alltoallv(myout, rcnt*nout, rdsp*nout, MPI_DOUBLE_PRECISION, &
                     bufout, scnt*nout, sdsp*nout, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
end subroutine balanced_compute

end module samasp_balance
