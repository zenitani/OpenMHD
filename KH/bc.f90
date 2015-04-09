subroutine bc(U,ix,jx)
  implicit none
  include 'param.h'
  integer, intent(in) :: ix, jx
  real(8) :: U(ix,jx,var1)
!----------------------------------------------------------------------

! periodic BC
  U(ix,:,:) = U(   2,:,:)
  U( 1,:,:) = U(ix-1,:,:)

! top/bottom walls
  U(:,1, ro) = U(:,2, ro)
  U(:,1, mx) = U(:,2, mx)
  U(:,1, my) =-U(:,2, my)
  U(:,1, mz) = U(:,2, mz)
  U(:,1, en) = U(:,2, en)
  U(:,1, bx) = U(:,2, bx)
  U(:,1, by) =-U(:,2, by)
  U(:,1, bz) = U(:,2, bz)
  U(:,1, ps) = U(:,2, ps)

  U(:,jx,ro) = U(:,jx-1,ro)
  U(:,jx,mx) = U(:,jx-1,mx)
  U(:,jx,my) =-U(:,jx-1,my)
  U(:,jx,mz) = U(:,jx-1,mz)
  U(:,jx,en) = U(:,jx-1,en)
  U(:,jx,bx) = U(:,jx-1,bx)
  U(:,jx,by) =-U(:,jx-1,by)
  U(:,jx,bz) = U(:,jx-1,bz)
  U(:,jx,ps) = U(:,jx-1,ps)

  return
end subroutine bc
