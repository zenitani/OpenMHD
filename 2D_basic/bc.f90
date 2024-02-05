subroutine bc(U,ix,jx)
  implicit none
  include 'param.h'
  integer, intent(in) :: ix, jx
  real(8), intent(inout) :: U(ix,jx,var1)
!----------------------------------------------------------------------

! west/east
  U(ix,:,:) = U(2,:,:)
  U(1,:,:)  = U(ix-1,:,:)

! south/north
  U(:,jx,:) = U(:,2,:)
  U(:,1,:)  = U(:,jx-1,:)

end subroutine bc
