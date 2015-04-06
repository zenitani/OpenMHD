subroutine mpibc3(U,ix,jx,myrank,npe)
!-----------------------------------------------------------------------
!     Boundary conditions for the parallel code 3 (mainp3)
!-----------------------------------------------------------------------
  implicit none
  include 'mpif.h'
  include 'param.h'
  integer, intent(in) :: ix, jx, myrank, npe
  real(8) :: U(ix,jx,var1)
  integer :: mstatus(mpi_status_size)
  integer :: mmx, mmx2, merr, mright, mleft, mpair
  real(8) :: bufsnd(jx,var1),  bufrcv(jx,var1)
  real(8) :: bufsnd2(ix,var1), bufrcv2(ix,var1)

  mmx  = jx*var1
  mmx2 = ix*var1

!----------------------------------------------------------------------
!  leftward transfer :  PE(myrank-1)  <---  PE(myrank)
!----------------------------------------------------------------------

  mright = myrank+1
  mleft  = myrank-1
  if (myrank.eq.(npe-1))  mright  = 0
  if (myrank.eq.0      )  mleft   = (npe-1)

  bufsnd(:,:) = U(2,:,:)

  call mpi_barrier(mpi_comm_world,merr)
  call mpi_sendrecv( &
       bufsnd,mmx,mpi_double_precision,mleft ,0, &
       bufrcv,mmx,mpi_double_precision,mright,0, &
       mpi_comm_world,mstatus,merr)
  call mpi_barrier(mpi_comm_world,merr)
      
  U(ix,:,:) = bufrcv(:,:)

!----------------------------------------------------------------------
!  rightward transfer :  PE(myrank)  --->  PE(myrank+1)
!----------------------------------------------------------------------

  mright = myrank+1
  mleft  = myrank-1
  if (myrank.eq.(npe-1))  mright  = 0
  if (myrank.eq.0      )  mleft   = (npe-1)

  bufsnd(:,:) = U(ix-1,:,:)

  call mpi_barrier(mpi_comm_world,merr)
  call mpi_sendrecv( &
       bufsnd,mmx,mpi_double_precision,mright,1, &
       bufrcv,mmx,mpi_double_precision,mleft ,1, &
       mpi_comm_world,mstatus,merr)

  U(1,:,:) = bufrcv(:,:)

!----------------------------------------------------------------------
!  top/bottom walls :  PE(myrank)  <--->  PE(mpair)
!----------------------------------------------------------------------

  mpair = (npe/2) + npe - myrank - 1
  if (mpair.ge.npe)  mpair = mpair - npe

  bufsnd2(:,:) = U(ix:1:-1,jx-1,:)
  call mpi_barrier(mpi_comm_world,merr)
  call mpi_sendrecv( &
       bufsnd2,mmx2,mpi_double_precision,mpair,0, &
       bufrcv2,mmx2,mpi_double_precision,mpair,0, &
       mpi_comm_world,mstatus,merr)
  U(:,1,:) = bufrcv2(:,:)
  U(:,1,mx) = -U(:,1,mx)
  U(:,1,bx) = -U(:,1,bx)

  bufsnd2(:,:) = U(ix:1:-1,2,:)
  call mpi_barrier(mpi_comm_world,merr)
  call mpi_sendrecv( &
       bufsnd2,mmx2,mpi_double_precision,mpair,1, &
       bufrcv2,mmx2,mpi_double_precision,mpair,1, &
       mpi_comm_world,mstatus,merr)
  U(:,jx,:) = bufrcv2(:,:)
  U(:,jx,mx) = -U(:,jx,mx)
  U(:,jx,bx) = -U(:,jx,bx)

!----------------------------------------------------------------------
!  top/bottom walls
!----------------------------------------------------------------------
!
!  U(:, 1,:)  =  U(:,   2,:)
!  U(:,jx,:)  =  U(:,jx-1,:)

  return
end subroutine mpibc3
