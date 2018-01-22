subroutine mpibc_vlvr_f(VL,VR,ix,jx,myrank,npe)
!-----------------------------------------------------------------------
!     2010/01/27  S. Zenitani   MPI ex for HLL
!-----------------------------------------------------------------------
  implicit none
  include 'mpif.h'
  include 'param.h'
  integer, intent(in) :: ix, jx, myrank, npe
! left flux (VL) [input]
  real(8) :: VL(ix,jx,var1)
! right flux (VR) [input]
  real(8) :: VR(ix,jx,var1)

  integer :: mstatus(mpi_status_size)
  integer :: mmx, merr, mright, mleft
  real(8) :: bufsnd(jx,var1), bufrcv(jx,var1)

  mmx = jx*var1

!----------------------------------------------------------------------
!  from PE(myrank) to PE(myrank-1) for new da(ix)
!----------------------------------------------------------------------

  mright = myrank+1
  mleft  = myrank-1
  if( myrank == npe-1 )  mright = 0
  if( myrank == 0     )  mleft  = npe-1

  bufsnd(:,:) = VR(1,:,:)
! call mpi_barrier(mpi_comm_world,merr)
  call mpi_sendrecv( &
       bufsnd,mmx,mpi_real8,mleft ,0, &
       bufrcv,mmx,mpi_real8,mright,0, &
       mpi_comm_world,mstatus,merr)
  VR(ix-1,:,:) = bufrcv(:,:)

!----------------------------------------------------------------------
!  from PE(myrank) to PE(myrank+1) for new da(1)
!----------------------------------------------------------------------

  bufsnd(:,:) = VL(ix-1,:,:)
! call mpi_barrier(mpi_comm_world,merr)
  call mpi_sendrecv( &
       bufsnd,mmx,mpi_real8,mright,1, &
       bufrcv,mmx,mpi_real8,mleft ,1, &
       mpi_comm_world,mstatus,merr)
  VL(1,:,:) = bufrcv(:,:)
!  call mpi_barrier(mpi_comm_world,merr)

  return
end subroutine mpibc_vlvr_f
