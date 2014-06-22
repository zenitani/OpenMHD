subroutine bc(U,ix,jx)
  implicit none
  include 'param.h'
  integer :: ix, jx
  real(8) :: U(ix,jx,var1)
!----------------------------------------------------------------------

! left/right
  U(ix,:,:) = U(2,:,:)
  U(1,:,:)  = U(ix-1,:,:)

! top/bottom
  U(:,jx,:) = U(:,2,:)
  U(:,1,:)  = U(:,jx-1,:)

  return
end subroutine bc
