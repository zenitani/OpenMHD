subroutine mpibc_vlvr_f2(VL,VR,ix,jx,myrank,npe)
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

  mright=  myrank+1
  mleft =  myrank-1
  if (myrank.eq.(npe-1)) mright  = mpi_proc_null
  if (myrank.eq.0      ) mleft   = mpi_proc_null

  bufsnd(:,:)=VR(1,:,:)

  call mpi_barrier(mpi_comm_world,merr)
  call mpi_sendrecv( &
       bufsnd,mmx,mpi_double_precision,mleft ,0, &
       bufrcv,mmx,mpi_double_precision,mright,0, &
       mpi_comm_world,mstatus,merr)
      
  if (myrank.ne.(npe-1)) then
     VR(ix-1,:,:)=bufrcv(:,:)
  endif

!----------------------------------------------------------------------
!  from PE(myrank) to PE(myrank+1) for new da(1)
!----------------------------------------------------------------------

  mright=  myrank+1
  mleft =  myrank-1
  if (myrank.eq.(npe-1)) mright  = mpi_proc_null
  if (myrank.eq.0      ) mleft   = mpi_proc_null

  bufsnd(:,:)=VL(ix-1,:,:)

  call mpi_barrier(mpi_comm_world,merr)
  call mpi_sendrecv( &
       bufsnd,mmx,mpi_double_precision,mright,1, &
       bufrcv,mmx,mpi_double_precision,mleft ,1, &
       mpi_comm_world,mstatus,merr)

  if (myrank.ne.0) then
     VL(1,:,:)=bufrcv(:,:)
  endif
  if (myrank.eq.0) then
     VL(1,:, ro) =  VR(1,:, ro)
     VL(1,:, mx) = -VR(1,:, mx)
     VL(1,:, my) =  VR(1,:, my)
     VL(1,:, mz) = -VR(1,:, mz)
     VL(1,:, en) =  VR(1,:, en)
     VL(1,:, bx) =  VR(1,:, bx)
     VL(1,:, by) = -VR(1,:, by)
     VL(1,:, bz) =  VR(1,:, bz)
     VL(1,:, ps) = -VR(1,:, ps)
  endif
!  call mpi_barrier(mpi_comm_world,merr)

  return
end subroutine mpibc_vlvr_f2
