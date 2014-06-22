subroutine bc(U,ix,jx)
  implicit none
  include 'param.h'
  integer :: ix, jx
  real(8) :: U(ix,jx,var1)
!----------------------------------------------------------------------

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

  return
end subroutine bc
