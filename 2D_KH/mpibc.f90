!-----------------------------------------------------------------------
!    MPI BC routines
!-----------------------------------------------------------------------
! Note that the parallel module automatically takes care of periodic BC.
! This file focuses on non-periodic BC operations.

subroutine mpibc_for_U(U,ix,jx)
  use parallel
  implicit none
  include 'param.h'
  integer, intent(in) :: ix, jx
  real(8), intent(inout) :: U(ix,jx,var1)

  ! south boundary
  if( ranks%south == mpi_proc_null ) then
     U(:,1,ro) =  U(:,2,ro)
     U(:,1,mx) =  U(:,2,mx)
     U(:,1,my) = -U(:,2,my)
     U(:,1,mz) =  U(:,2,mz)
     U(:,1,en) =  U(:,2,en)
     U(:,1,bx) =  U(:,2,bx)
     U(:,1,by) = -U(:,2,by)
     U(:,1,bz) =  U(:,2,bz)
     U(:,1,ps) =  U(:,2,ps)
  endif

  ! north boundary
  if( ranks%north == mpi_proc_null ) then
     U(:,jx,ro) =  U(:,jx-1,ro)
     U(:,jx,mx) =  U(:,jx-1,mx)
     U(:,jx,my) = -U(:,jx-1,my)
     U(:,jx,mz) =  U(:,jx-1,mz)
     U(:,jx,en) =  U(:,jx-1,en)
     U(:,jx,bx) =  U(:,jx-1,bx)
     U(:,jx,by) = -U(:,jx-1,by)
     U(:,jx,bz) =  U(:,jx-1,bz)
     U(:,jx,ps) =  U(:,jx-1,ps)
  endif

end subroutine mpibc_for_U


subroutine mpibc_for_G(VL,VR,ix,jx)
  use parallel
  implicit none
  include 'param.h'
  integer, intent(in) :: ix, jx
! left flux (VL) [input]
  real(8), intent(inout) :: VL(ix,jx,var1)
! right flux (VR) [input]
  real(8), intent(inout) :: VR(ix,jx,var1)

! south boundary
  if( ranks%south == mpi_proc_null ) then
     VL(:,1,ro) =  VR(:,1,ro)
     VL(:,1,vx) =  VR(:,1,vx)
     VL(:,1,vy) = -VR(:,1,vy)
     VL(:,1,vz) =  VR(:,1,vz)
     VL(:,1,en) =  VR(:,1,en)
     VL(:,1,bx) =  VR(:,1,bx)
     VL(:,1,by) = -VR(:,1,by)
     VL(:,1,bz) =  VR(:,1,bz)
     VL(:,1,ps) =  VR(:,1,ps)
  endif

  ! north boundary
  if( ranks%north == mpi_proc_null ) then
     VR(:,jx-1,ro) =  VL(:,jx-1,ro)
     VR(:,jx-1,vx) =  VL(:,jx-1,vx)
     VR(:,jx-1,vy) = -VL(:,jx-1,vy)
     VR(:,jx-1,vz) =  VL(:,jx-1,vz)
     VR(:,jx-1,en) =  VL(:,jx-1,en)
     VR(:,jx-1,bx) =  VL(:,jx-1,bx)
     VR(:,jx-1,by) = -VL(:,jx-1,by)
     VR(:,jx-1,bz) =  VL(:,jx-1,bz)
     VR(:,jx-1,ps) =  VL(:,jx-1,ps)
  endif

end subroutine mpibc_for_G
