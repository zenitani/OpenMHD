subroutine mpibc(U,ix,jx,myrank,npe)
  implicit none
  include 'mpif.h'
  include 'param.h'
  integer :: ix, jx, myrank, npe
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
  if (myrank.eq.(npe-1))  mright  = 0
  if (myrank.eq.0      )  mleft   = npe-1

  bufsnd(:,:) = U(2,:,:)

  call mpi_barrier(mpi_comm_world,merr)
  call mpi_sendrecv( &
       bufsnd,mmx,mpi_double_precision,mleft ,0, &
       bufrcv,mmx,mpi_double_precision,mright,0, &
       mpi_comm_world,mstatus,merr)
      
  U(ix,:,:) = bufrcv(:,:)

!----------------------------------------------------------------------
!  rightward transfer :  PE(myrank)  --->  PE(myrank+1)
!----------------------------------------------------------------------

  mright = myrank+1
  mleft  = myrank-1
  if (myrank.eq.(npe-1))  mright  = 0
  if (myrank.eq.0      )  mleft   = npe-1

  bufsnd(:,:) = U(ix-1,:,:)

  call mpi_barrier(mpi_comm_world,merr)
  call mpi_sendrecv( &
       bufsnd,mmx,mpi_double_precision,mright,1, &
       bufrcv,mmx,mpi_double_precision,mleft ,1, &
       mpi_comm_world,mstatus,merr)

  U(1,:,:) = bufrcv(:,:)

!----------------------------------------------------------------------
!  top/bottom walls
!----------------------------------------------------------------------

  U(:,1, ro) = U(:,2, ro)
  U(:,1, mx) = U(:,2, mx)
  U(:,1, my) =-U(:,2, my)
  U(:,1, mz) = U(:,2, mz)
  U(:,1, en) = U(:,2, en)
  U(:,1, bx) = U(:,2, bx)
  U(:,1, by) =-U(:,2, by)
  U(:,1, bz) = U(:,2, bz)
  U(:,1, ps) = U(:,2, ps)

  U(:,jx,ro) = U(:,jx-1,ro)
  U(:,jx,mx) = U(:,jx-1,mx)
  U(:,jx,my) =-U(:,jx-1,my)
  U(:,jx,mz) = U(:,jx-1,mz)
  U(:,jx,en) = U(:,jx-1,en)
  U(:,jx,bx) = U(:,jx-1,bx)
  U(:,jx,by) =-U(:,jx-1,by)
  U(:,jx,bz) = U(:,jx-1,bz)
  U(:,jx,ps) = U(:,jx-1,ps)

  return
end subroutine mpibc
