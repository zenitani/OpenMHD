subroutine mpibc(U,ix,jx,myrank,npe)
  implicit none
  include 'mpif.h'
  include 'param.h'
  integer, intent(in) :: ix, jx, myrank, npe
  real(8) :: U(ix,jx,var1)
  integer :: mstatus(mpi_status_size)
  integer :: mmx, merr, mright, mleft
  real(8) :: bufsnd(jx,var1), bufrcv(jx,var1)

  mmx = jx*var1

!----------------------------------------------------------------------
!  leftward transfer :  PE(myrank-1)  <---  PE(myrank)
!----------------------------------------------------------------------

  mright = myrank+1
  mleft  = myrank-1
  if( myrank == npe-1 )  mright = 0
  if( myrank == 0     )  mleft  = npe-1

  bufsnd(:,:) = U(2,:,:)
! call mpi_barrier(mpi_comm_world,merr)
  call mpi_sendrecv( &
       bufsnd,mmx,mpi_real8,mleft ,0, &
       bufrcv,mmx,mpi_real8,mright,0, &
       mpi_comm_world,mstatus,merr)
  U(ix,:,:) = bufrcv(:,:)

!----------------------------------------------------------------------
!  rightward transfer :  PE(myrank)  --->  PE(myrank+1)
!----------------------------------------------------------------------

  bufsnd(:,:) = U(ix-1,:,:)
! call mpi_barrier(mpi_comm_world,merr)
  call mpi_sendrecv( &
       bufsnd,mmx,mpi_real8,mright,1, &
       bufrcv,mmx,mpi_real8,mleft ,1, &
       mpi_comm_world,mstatus,merr)
  U(1,:,:) = bufrcv(:,:)

!----------------------------------------------------------------------
!  top/bottom walls
!----------------------------------------------------------------------

  U(:,jx,:) = U(:,2,:)
  U(:,1,:)  = U(:,jx-1,:)

  return
end subroutine mpibc
