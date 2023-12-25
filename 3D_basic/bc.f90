subroutine bc(U,ix,jx,kx)
  implicit none
  include 'param.h'
  integer, intent(in) :: ix, jx, kx
  real(8), intent(inout) :: U(ix,jx,kx,var1)
!----------------------------------------------------------------------

! east/west
  U(ix,:,:,:) = U(   2,:,:,:)
  U( 1,:,:,:) = U(ix-1,:,:,:)

! north/south
  U(:,jx,:,:) = U(:,   2,:,:)
  U(:, 1,:,:) = U(:,jx-1,:,:)

! upstairs/downstairs
  U(:,:,kx,:) = U(:,:,   2,:)
  U(:,:, 1,:) = U(:,:,kx-1,:)

end subroutine bc
