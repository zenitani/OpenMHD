subroutine mpibc_vlvr_g3(VL,VR,ix,jx,myrank,npe)
!-----------------------------------------------------------------------
!     2010/01/27  S. Zenitani   MPI ex for HLL
!-----------------------------------------------------------------------
  implicit none
  include 'mpif.h'
  include 'param.h'
  integer :: ix, jx, myrank, npe
! left flux (VL) [input]
  real(8) :: VL(ix,jx,var1)
! right flux (VR) [input]
  real(8) :: VR(ix,jx,var1)

  integer :: mstatus(mpi_status_size)
  integer :: mmx, merr, mpair
  real(8) :: bufsnd(ix,var1), bufrcv(ix,var1)

  mmx = ix*var1

  mpair = (npe/2) + npe - myrank - 1
  if (mpair.ge.npe)  mpair = mpair - npe

  bufsnd(:,:) = VL(ix:1:-1,jx-1,:)
  call mpi_barrier(mpi_comm_world,merr)
  call mpi_sendrecv( &
       bufsnd,mmx,mpi_double_precision,mpair,0, &
       bufrcv,mmx,mpi_double_precision,mpair,0, &
       mpi_comm_world,mstatus,merr)
  VL(:,1,:)  = bufrcv(:,:)
  VL(:,1,vx) = -VL(:,1,vx)
  VL(:,1,bx) = -VL(:,1,bx)
  VL(:,1,ps) = -VL(:,1,ps)

  bufsnd(:,:) = VR(ix:1:-1,1,:)
  call mpi_barrier(mpi_comm_world,merr)
  call mpi_sendrecv( &
       bufsnd,mmx,mpi_double_precision,mpair,1, &
       bufrcv,mmx,mpi_double_precision,mpair,1, &
       mpi_comm_world,mstatus,merr)
  VR(:,jx-1,:)  = bufrcv(:,:)
  VR(:,jx-1,vx) = -VR(:,jx-1,vx)
  VR(:,jx-1,bx) = -VR(:,jx-1,bx)
  VR(:,jx-1,ps) = -VR(:,jx-1,ps)

  return
end subroutine mpibc_vlvr_g3
