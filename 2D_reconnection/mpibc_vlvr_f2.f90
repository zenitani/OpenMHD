subroutine mpibc_vlvr_f2(VL,VR,ix,jx,myrank,npe)
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
  if( myrank == npe-1 )  mright = mpi_proc_null
  if( myrank == 0     )  mleft  = mpi_proc_null

! nonblocking communication (mreq1)
  call mpi_irecv(bufrcv1,mmx,mpi_real8,mright,0,mpi_comm_world,mreq1(1),merr)
  bufsnd1(:,:) = VR(1,:,:)
  call mpi_isend(bufsnd1,mmx,mpi_real8,mleft ,0,mpi_comm_world,mreq1(2),merr)

! nonblocking communication (mreq2)
  call mpi_irecv(bufrcv2,mmx,mpi_real8,mleft ,1,mpi_comm_world,mreq2(1),merr)
  bufsnd2(:,:) = VL(ix-1,:,:)
  call mpi_isend(bufsnd2,mmx,mpi_real8,mright,1,mpi_comm_world,mreq2(2),merr)

! wait for completion (mreq1)
  call mpi_waitall(2,mreq1,mpi_statuses_ignore,merr)
! right boundary
  if( myrank /= npe-1 ) then
     VR(ix-1,:,:) = bufrcv1(:,:)
  else
     VR(ix-1,:, ro) =  VL(ix-1,:, ro)
     VR(ix-1,:, vx) = -VL(ix-1,:, vx)
     VR(ix-1,:, vy) =  VL(ix-1,:, vy)
     VR(ix-1,:, vz) = -VL(ix-1,:, vz)
     VR(ix-1,:, pr) =  VL(ix-1,:, pr)
     VR(ix-1,:, bx) =  VL(ix-1,:, bx)
     VR(ix-1,:, by) = -VL(ix-1,:, by)
     VR(ix-1,:, bz) =  VL(ix-1,:, bz)
     VR(ix-1,:, ps) = -VL(ix-1,:, ps)
  endif

! wait for completion (mreq2)
  call mpi_waitall(2,mreq2,mpi_statuses_ignore,merr)
! left boundary
  if( myrank /= 0 ) then
     VL(1,:,:) = bufrcv2(:,:)
  else
     VL(1,:, ro) =  VR(1,:, ro)
     VL(1,:, vx) = -VR(1,:, vx)
     VL(1,:, vy) =  VR(1,:, vy)
     VL(1,:, vz) = -VR(1,:, vz)
     VL(1,:, pr) =  VR(1,:, pr)
     VL(1,:, bx) =  VR(1,:, bx)
     VL(1,:, by) = -VR(1,:, by)
     VL(1,:, bz) =  VR(1,:, bz)
     VL(1,:, ps) = -VR(1,:, ps)
  endif

  return
end subroutine mpibc_vlvr_f2
