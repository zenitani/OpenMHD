subroutine mpibc(U,ix,jx,myrank,npe)
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
  if( myrank == npe-1 )  mright = 0
  if( myrank == 0     )  mleft  = npe-1

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
  U(ix,:,:) = bufrcv1(:,:)

! wait for completion (mreq2)
  call mpi_waitall(2,mreq2,mpi_statuses_ignore,merr)
  U(1,:,:) = bufrcv2(:,:)

! top/bottom boundaries
  U(:,jx,:) = U(:,2,:)
  U(:,1,:)  = U(:,jx-1,:)

  return
end subroutine mpibc
