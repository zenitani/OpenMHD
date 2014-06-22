subroutine mpibc(U,ix,jx,myrank,npe)
!-----------------------------------------------------------------------
!     Boundary conditions for the parallel code (mainp)
!-----------------------------------------------------------------------
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
  if (myrank.eq.(npe-1))  mright  = mpi_proc_null
  if (myrank.eq.0      )  mleft   = mpi_proc_null

  bufsnd(:,:) = U(2,:,:)

  call mpi_barrier(mpi_comm_world,merr)
  call mpi_sendrecv( &
       bufsnd,mmx,mpi_double_precision,mleft ,0, &
       bufrcv,mmx,mpi_double_precision,mright,0, &
       mpi_comm_world,mstatus,merr)
      
  if( myrank.ne.(npe-1) ) then
     U(ix,:,:) = bufrcv(:,:)
  endif

!----------------------------------------------------------------------
!  rightward transfer :  PE(myrank)  --->  PE(myrank+1)
!----------------------------------------------------------------------

  mright = myrank+1
  mleft  = myrank-1
  if (myrank.eq.(npe-1))  mright  = mpi_proc_null
  if (myrank.eq.0      )  mleft   = mpi_proc_null

  bufsnd(:,:) = U(ix-1,:,:)

  call mpi_barrier(mpi_comm_world,merr)
  call mpi_sendrecv( &
       bufsnd,mmx,mpi_double_precision,mright,1, &
       bufrcv,mmx,mpi_double_precision,mleft ,1, &
       mpi_comm_world,mstatus,merr)

  if( myrank.ne.0 ) then
     U(1,:,:) = bufrcv(:,:)
  endif

!----------------------------------------------------------------------
!  outflow BC
!----------------------------------------------------------------------

  if( myrank.eq.0 ) then
     U(1,:,:) = U(2,:,:)
  endif
  if( myrank.eq.(npe-1) ) then
     U(ix,:,:) = U(ix-1,:,:)
  endif

!----------------------------------------------------------------------
!  top/bottom walls
!----------------------------------------------------------------------

  U(:, 1,:)  =  U(:,   2,:)
  U(:,jx,:)  =  U(:,jx-1,:)
  U(:, 1,my) = -U(:,   2,my)
  U(:,jx,my) = -U(:,jx-1,my)
  U(:, 1,by) = -U(:,   2,by) 
  U(:,jx,by) = -U(:,jx-1,by) 

  return
end subroutine mpibc
