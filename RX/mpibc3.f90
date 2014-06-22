subroutine mpibc3(U,ix,jx,myrank,npe)
!-----------------------------------------------------------------------
!     Boundary conditions for the parallel code 3 (mainp3)
!-----------------------------------------------------------------------
  implicit none
  include 'mpif.h'
  include 'param.h'
  integer :: ix, jx, myrank, npe
  real(8) :: U(ix,jx,var1)
  integer :: mstatus(mpi_status_size)
  integer :: mmx, merr, mright, mleft
  integer :: mmx2, mpair
  real(8) :: bufsnd(jx,var1), bufrcv(jx,var1)
  real(8) :: bufsnd2(ix,var1), bufrcv2(ix,var1)

  mmx  = jx*var1
  mmx2 = ix*var1

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
!  bottom walls
!----------------------------------------------------------------------

  U(:,1, ro) =  U(:,2, ro)
  U(:,1, mx) =  U(:,2, mx)
  U(:,1, my) = -U(:,2, my)
  U(:,1, mz) = -U(:,2, mz)
  U(:,1, en) =  U(:,2, en)
  U(:,1, bx) = -U(:,2, bx)
  U(:,1, by) =  U(:,2, by)
  U(:,1, bz) =  U(:,2, bz)
  U(:,1, ps) = -U(:,2, ps)

!----------------------------------------------------------------------
!  top walls :  PE(myrank)  <--->  PE( npe-myrank-1 )
!----------------------------------------------------------------------

  mpair = npe - myrank - 1
!  write(6,*) npe, myrank, mpair
  bufsnd2(:,:) = U(ix:1:-1,jx-1,:)

  call mpi_barrier(mpi_comm_world,merr)
  call mpi_sendrecv( &
       bufsnd2,mmx2,mpi_double_precision,mpair,0, &
       bufrcv2,mmx2,mpi_double_precision,mpair,0, &
       mpi_comm_world,mstatus,merr)
      
  U(:,jx,:)  = bufrcv2(:,:)
!  U(:,jx,ro) =  U(:,jx,ro)
  U(:,jx,mx) = -U(:,jx,mx)
  U(:,jx,my) = -U(:,jx,my)
!  U(:,jx,mz) =  U(:,jx,mz) ! not sure this one
!  U(:,jx,en) =  U(:,jx,en)
!  U(:,jx,bx) =  U(:,jx,bx)
!  U(:,jx,by) =  U(:,jx,by)
!  U(:,jx,bz) =  U(:,jx,bz) ! not sure this one
  U(:,jx,ps) = -U(:,jx,ps)

  return
end subroutine mpibc3
