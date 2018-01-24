subroutine mpibc_vlvr_f(VL,VR,ix,jx,myrank,npe)
  implicit none
  include 'mpif.h'
  include 'param.h'
  integer, intent(in) :: ix, jx, myrank, npe
! left flux (VL) [input]
  real(8), intent(inout) :: VL(ix,jx,var1)
! right flux (VR) [input]
  real(8), intent(inout) :: VR(ix,jx,var1)
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
  bufsnd1(:,:) = VR(1,:,:)
  call mpi_isend(bufsnd1,mmx,mpi_real8,mleft ,0,mpi_comm_world,mreq1(2),merr)

! nonblocking communication (mreq2)
  call mpi_irecv(bufrcv2,mmx,mpi_real8,mleft ,1,mpi_comm_world,mreq2(1),merr)
  bufsnd2(:,:) = VL(ix-1,:,:)
  call mpi_isend(bufsnd2,mmx,mpi_real8,mright,1,mpi_comm_world,mreq2(2),merr)

! wait for completetion (mreq1)
  call mpi_waitall(2,mreq1,mpi_statuses_ignore,merr)
  VR(ix-1,:,:) = bufrcv1(:,:)

! wait for completetion (mreq2)
  call mpi_waitall(2,mreq2,mpi_statuses_ignore,merr)
  VL(1,:,:) = bufrcv2(:,:)

  return
end subroutine mpibc_vlvr_f
