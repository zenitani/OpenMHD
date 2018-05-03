!-----------------------------------------------------------------------
!     BC routines
!-----------------------------------------------------------------------

! BC for U
subroutine bc_for_U(U,ix,jx)
  implicit none
  include 'param.h'
  integer, intent(in) :: ix, jx
  real(8), intent(inout) :: U(ix,jx,var1)

! left/right
!  U(ix,:,:) = U(2,:,:)
!  U(1,:,:)  = U(ix-1,:,:)
  U(1,:,:)  = U(2,:,:)
  U(ix,:,:) = U(ix-1,:,:)

! top/bottom (lazy BC)
  U(:, 1,:)  =  U(:,   2,:)
  U(:,jx,:)  =  U(:,jx-1,:)
  U(:, 1,my) = -U(:,   2,my)
  U(:,jx,my) = -U(:,jx-1,my)
  U(:, 1,by) = -U(:,   2,by) 
  U(:,jx,by) = -U(:,jx-1,by) 

end subroutine bc_for_U


! BC for G -- VL/VR in the Y direction
subroutine bc_for_G( VL,VR,ix,jx )
  implicit none
  include 'param.h'
  integer, intent(in) :: ix, jx
  real(8), intent(inout) :: VL(ix,jx,var1), VR(ix,jx,var1) ! interpolated states

! bottom boundary
  VL(:,1,ro) =  VR(:,1,ro)
  VL(:,1,vx) =  VR(:,1,vx)
  VL(:,1,vy) = -VR(:,1,vy)
  VL(:,1,vz) =  VR(:,1,vz)
  VL(:,1,pr) =  VR(:,1,pr)
  VL(:,1,bx) =  VR(:,1,bx)
  VL(:,1,by) = -VR(:,1,by)
  VL(:,1,bz) =  VR(:,1,bz)
  VL(:,1,ps) =  VR(:,1,ps)

! top boundary
  VR(:,jx-1,ro) =  VL(:,jx-1,ro)
  VR(:,jx-1,vx) =  VL(:,jx-1,vx)
  VR(:,jx-1,vy) = -VL(:,jx-1,vy)
  VR(:,jx-1,vz) =  VL(:,jx-1,vz)
  VR(:,jx-1,pr) =  VL(:,jx-1,pr)
  VR(:,jx-1,bx) =  VL(:,jx-1,bx)
  VR(:,jx-1,by) = -VL(:,jx-1,by)
  VR(:,jx-1,bz) =  VL(:,jx-1,bz)
  VR(:,jx-1,ps) =  VL(:,jx-1,ps)

end subroutine bc_for_G
