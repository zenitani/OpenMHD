!-----------------------------------------------------------------------
!     Boundary conditions (parallel version)
!-----------------------------------------------------------------------
! Note that the parallel module automatically takes care of periodic BC.
! This file focuses on non-periodic BC operations.

! BC for U
subroutine mpibc_for_U(U,ix,jx)
  use parallel
  implicit none
  include 'param.h'
  integer, intent(in) :: ix, jx
  real(8), intent(inout) :: U(ix,jx,var1)

  ! west boundary
  if( ranks%west == mpi_proc_null ) then
     U(1,:,ro) =  U(2,:,ro)
     U(1,:,mx) = -U(2,:,mx)
     U(1,:,my) =  U(2,:,my)
     U(1,:,mz) = -U(2,:,mz)
     U(1,:,en) =  U(2,:,en)
     U(1,:,bx) =  U(2,:,bx)
     U(1,:,by) = -U(2,:,by)
     U(1,:,bz) =  U(2,:,bz)
     U(1,:,ps) = -U(2,:,ps)
  endif

  ! east boundary
  if( ranks%east == mpi_proc_null ) then
     U(ix,:,ro) =  U(ix-1,:,ro)
     U(ix,:,mx) = -U(ix-1,:,mx)
     U(ix,:,my) =  U(ix-1,:,my)
     U(ix,:,mz) = -U(ix-1,:,mz)
     U(ix,:,en) =  U(ix-1,:,en)
     U(ix,:,bx) =  U(ix-1,:,bx)
     U(ix,:,by) = -U(ix-1,:,by)
     U(ix,:,bz) =  U(ix-1,:,bz)
     U(ix,:,ps) = -U(ix-1,:,ps)
  endif

  ! south boundary
  if( ranks%south == mpi_proc_null ) then
     U(:,1,ro) =  U(:,2,ro)
     U(:,1,mx) =  U(:,2,mx)
     U(:,1,my) = -U(:,2,my)
     U(:,1,mz) = -U(:,2,mz)
     U(:,1,en) =  U(:,2,en)
     U(:,1,bx) = -U(:,2,bx)
     U(:,1,by) =  U(:,2,by)
     U(:,1,bz) =  U(:,2,bz)
     U(:,1,ps) = -U(:,2,ps)
  endif

  ! north boundary
  if( ranks%north == mpi_proc_null ) then
     U(:,jx,:)  =  U(:,jx-1,:)
     U(:,jx,my) = -U(:,jx-1,my)
     U(:,jx,by) = -U(:,jx-1,by)
  endif

end subroutine mpibc_for_U


! BC for F -- VL/VR in the X direction
subroutine mpibc_for_F( VL,VR,ix,jx )
  use parallel
  implicit none
  include 'param.h'
  integer, intent(in) :: ix, jx
  real(8), intent(inout) :: VL(ix,jx,var1), VR(ix,jx,var1) ! interpolated states

! west boundary
  if( ranks%west == mpi_proc_null ) then
     VL(1,:,ro) =  VR(1,:,ro)
     VL(1,:,vx) = -VR(1,:,vx)
     VL(1,:,vy) =  VR(1,:,vy)
     VL(1,:,vz) = -VR(1,:,vz)
     VL(1,:,pr) =  VR(1,:,pr)
     VL(1,:,bx) =  VR(1,:,bx)
     VL(1,:,by) = -VR(1,:,by)
     VL(1,:,bz) =  VR(1,:,bz)
     VL(1,:,ps) = -VR(1,:,ps)
  endif

! east boundary
  if( ranks%east == mpi_proc_null ) then
     VR(ix-1,:,ro) =  VL(ix-1,:,ro)
     VR(ix-1,:,vx) = -VL(ix-1,:,vx)
     VR(ix-1,:,vy) =  VL(ix-1,:,vy)
     VR(ix-1,:,vz) = -VL(ix-1,:,vz)
     VR(ix-1,:,pr) =  VL(ix-1,:,pr)
     VR(ix-1,:,bx) =  VL(ix-1,:,bx)
     VR(ix-1,:,by) = -VL(ix-1,:,by)
     VR(ix-1,:,bz) =  VL(ix-1,:,bz)
     VR(ix-1,:,ps) = -VL(ix-1,:,ps)
  endif

end subroutine mpibc_for_F


! BC for G -- VL/VR in the Y direction
subroutine mpibc_for_G( VL,VR,ix,jx )
  use parallel
  implicit none
  include 'param.h'
  integer, intent(in) :: ix, jx
  real(8), intent(inout) :: VL(ix,jx,var1), VR(ix,jx,var1) ! interpolated states

! south boundary
  if( ranks%south == mpi_proc_null ) then
     VL(:,1,ro) =  VR(:,1,ro)
     VL(:,1,vx) =  VR(:,1,vx)
     VL(:,1,vy) = -VR(:,1,vy)
     VL(:,1,vz) = -VR(:,1,vz)
     VL(:,1,pr) =  VR(:,1,pr)
     VL(:,1,bx) = -VR(:,1,bx)
     VL(:,1,by) =  VR(:,1,by)
     VL(:,1,bz) =  VR(:,1,bz)
     VL(:,1,ps) = -VR(:,1,ps)
  endif

! north boundary
  if( ranks%north == mpi_proc_null ) then
     VR(:,jx-1,ro) =  VL(:,jx-1,ro)
     VR(:,jx-1,vx) =  VL(:,jx-1,vx)
     VR(:,jx-1,vy) = -VL(:,jx-1,vy)
     VR(:,jx-1,vz) =  VL(:,jx-1,vz)
     VR(:,jx-1,pr) =  VL(:,jx-1,pr)
     VR(:,jx-1,bx) =  VL(:,jx-1,bx)
     VR(:,jx-1,by) = -VL(:,jx-1,by)
     VR(:,jx-1,bz) =  VL(:,jx-1,bz)
     VR(:,jx-1,ps) =  VL(:,jx-1,ps)
  endif

end subroutine mpibc_for_G
