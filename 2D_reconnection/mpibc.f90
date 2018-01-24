subroutine mpibc(U,ix,jx,myrank,npe)
!-----------------------------------------------------------------------
!     Boundary conditions for the parallel code (mainp)
!-----------------------------------------------------------------------
  implicit none
  include 'mpif.h'
  include 'param.h'
  integer, intent(in) :: ix, jx, myrank, npe
  real(8), intent(inout) :: U(ix,jx,var1)
!----------------------------------------------------------------------
  integer :: mmx, merr, mright, mleft
  real(8) :: bufsnd1(jx,var1), bufrcv1(jx,var1)
  real(8) :: bufsnd2(jx,var1), bufrcv2(jx,var1)
  integer :: mreq1(2), mreq2(2)
!----------------------------------------------------------------------

  mmx    = jx*var1
  mright = myrank+1
  mleft  = myrank-1
  if( myrank == npe-1 )  mright = mpi_proc_null
  if( myrank == 0     )  mleft  = mpi_proc_null

! nonblocking communication (mreq1)
  call mpi_irecv(bufrcv1,mmx,mpi_real8,mright,0,mpi_comm_world,mreq1(1),merr)
  bufsnd1(:,:) = U(2,:,:)
  call mpi_isend(bufsnd1,mmx,mpi_real8,mleft ,0,mpi_comm_world,mreq1(2),merr)

! nonblocking communication (mreq2)
  call mpi_irecv(bufrcv2,mmx,mpi_real8,mleft ,1,mpi_comm_world,mreq2(1),merr)
  bufsnd2(:,:) = U(ix-1,:,:)
  call mpi_isend(bufsnd2,mmx,mpi_real8,mright,1,mpi_comm_world,mreq2(2),merr)

! wait for completion (mreq1)
  call mpi_waitall(2,mreq1,mpi_statuses_ignore,merr)
! right boundary
  if( myrank /= npe-1 ) then
     U(ix,:,:) = bufrcv1(:,:)
  else
     U(ix,:,:) = U(ix-1,:,:)
  endif

! wait for completion (mreq2)
  call mpi_waitall(2,mreq2,mpi_statuses_ignore,merr)
! left boundary
  if( myrank /= 0 ) then
     U(1,:,:) = bufrcv2(:,:)
  else
     U(1,:,:) = U(2,:,:)
  endif

! top/bottom boundaries
  U(:, 1,:)  =  U(:,   2,:)
  U(:,jx,:)  =  U(:,jx-1,:)
  U(:, 1,my) = -U(:,   2,my)
  U(:,jx,my) = -U(:,jx-1,my)
  U(:, 1,by) = -U(:,   2,by) 
  U(:,jx,by) = -U(:,jx-1,by) 

  return
end subroutine mpibc
