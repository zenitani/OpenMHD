subroutine modelp(U,V,x,y,dx,ix,jx)
  use parallel
  implicit none
  include 'param.h'
  real(8), intent(out) :: U(ix,jx,var1)
  real(8), intent(out) :: V(ix,jx,var2)
  real(8), intent(out) :: x(ix), y(jx), dx
  integer, intent(in)  :: ix, jx
! ---------------------------------------------------
  real(8), parameter :: pi = 4.d0*atan(1.d0)
! x & y positions (Note: domain_y(2) is automatically calculated)
  real(8), parameter :: domain_x(2) = (/0.d0, 2*pi/)
  real(8), parameter :: domain_y(1) = (/0.d0/)
! ---------------------------------------------------
  integer :: i, j
  integer :: iix, jjx
  real(8) :: B2, v2, f1
  real(8) :: tmpx(cart2d%sizes(1)*(ix-2) + 2)
  real(8) :: tmpy(cart2d%sizes(2)*(jx-2) + 2)
! ---------------------------------------------------

! grid in the X direction
  iix = cart2d%sizes(1)*(ix-2) + 2
  dx = ( domain_x(2) - domain_x(1) ) / dble( iix-2 )
  tmpx(1)   = domain_x(1) - dx/2
! tmpx(iix) = domain_x(2) + dx/2
  do i=2,iix
     tmpx(i) = tmpx(i-1) + dx
  enddo
  do i=1,ix
     x(i) = tmpx(cart2d%coords(1)*(ix-2) + i)
  enddo

! grid in the Y direction
  jjx = cart2d%sizes(2)*(jx-2) + 2
  tmpy(1) = domain_y(1) - dx/2
  do j=2,jjx
     tmpy(j) = tmpy(j-1) + dx
  enddo
  do j=1,jx
     y(j) = tmpy(cart2d%coords(2)*(jx-2) + j)
  enddo
! ---------------------------------------------------

  do j=1,jx
  do i=1,ix

     U(i,j,ro) = gamma**2

     V(i,j,pr) = gamma
     V(i,j,vx) = -sin(y(j))
     V(i,j,vy) = sin(x(i))
     V(i,j,vz) = 0.d0

     U(i,j,bx) = -sin(y(j))
     U(i,j,by) = sin(2*x(i))
     U(i,j,bz) = 0.d0
     U(i,j,ps) = 0.d0

     f1 = 1.d0 / ( gamma - 1 )
     v2 = dot_product( V(i,j,vx:vz), V(i,j,vx:vz) )
     B2 = dot_product( U(i,j,bx:bz), U(i,j,bx:bz) )
     U(i,j,mx:mz) = U(i,j,ro) * V(i,j,vx:vz)
     U(i,j,en)    = 0.5d0 * ( U(i,j,ro)*v2 + B2 ) + f1*V(i,j,pr)

  enddo
  enddo

  return
end subroutine modelp
