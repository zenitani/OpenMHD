subroutine mpibc2(U,ix,jx,myrank,npe)
!-----------------------------------------------------------------------
!     Boundary conditions for the parallel code 2 (mainp2)
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
!     U(1,:,:) = U(2,:,:)
     U(1,:, ro) =  U(2,:, ro)
     U(1,:, mx) = -U(2,:, mx)
     U(1,:, my) =  U(2,:, my)
     U(1,:, mz) = -U(2,:, mz)
     U(1,:, en) =  U(2,:, en)
     U(1,:, bx) =  U(2,:, bx)
     U(1,:, by) = -U(2,:, by)
     U(1,:, bz) =  U(2,:, bz)
     U(1,:, ps) = -U(2,:, ps)
  endif
  if( myrank.eq.(npe-1) ) then
!     U(ix,:,:) = U(ix-1,:,:)
     U(ix,:, ro) =  U(ix-1,:, ro)
     U(ix,:, mx) = -U(ix-1,:, mx)
     U(ix,:, my) =  U(ix-1,:, my)
     U(ix,:, mz) = -U(ix-1,:, mz)
     U(ix,:, en) =  U(ix-1,:, en)
     U(ix,:, bx) =  U(ix-1,:, bx)
     U(ix,:, by) = -U(ix-1,:, by)
     U(ix,:, bz) =  U(ix-1,:, bz)
     U(ix,:, ps) = -U(ix-1,:, ps)
  endif

!----------------------------------------------------------------------
!  top/bottom walls
!----------------------------------------------------------------------

! bottom
  U(:,1, ro) =  U(:,2, ro)
  U(:,1, mx) =  U(:,2, mx)
  U(:,1, my) = -U(:,2, my)
  U(:,1, mz) = -U(:,2, mz)
  U(:,1, en) =  U(:,2, en)
  U(:,1, bx) = -U(:,2, bx)
  U(:,1, by) =  U(:,2, by)
  U(:,1, bz) =  U(:,2, bz)
  U(:,1, ps) = -U(:,2, ps)

! top
  U(:,jx,:)  =  U(:,jx-1,:)
  U(:,jx,my) = -U(:,jx-1,my) 
  U(:,jx,by) = -U(:,jx-1,by)

  return
end subroutine mpibc2
