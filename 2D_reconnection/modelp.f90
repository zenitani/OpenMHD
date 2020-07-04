subroutine modelp(U,V,x,y,dx,ix,jx)
!-----------------------------------------------------------------------
!     Initial configuration for the parallel version (mainp2)
!-----------------------------------------------------------------------
  use parallel
  implicit none
  include 'param.h'
  real(8), intent(out) :: U(ix,jx,var1)
  real(8), intent(out) :: V(ix,jx,var2)
  real(8), intent(out) :: x(ix), y(jx), dx
  integer :: ix, jx
! ---------------------------------------------------
! x & y positions (Note: domain_y(2) is automatically calculated)
  real(8), parameter :: domain_x(2) = (/0.d0, 200.d0/)
  real(8), parameter :: domain_y(1) = (/0.d0/)
! plasma beta in the upstream region
  real(8), parameter :: beta = 0.2d0
! ---------------------------------------------------
  integer :: i, j
  integer :: iix, jjx
  real(8) :: B2, v2, f1, r2, b1
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

! Here we use min(..., 25) to avoid NaN on some systems.
! Note that (cosh(25))**(-2) = 7.7E-22 is extremely small.
     U(i,j,ro) = 1.d0 + cosh(min( y(j), 25.d0 ))**(-2) / beta
!     U(i,j,ro) = 1.d0 + cosh(y(j))**(-2) / beta
     V(i,j,pr) = 0.5d0 * beta * U(i,j,ro)
     V(i,j,vx) = 0.d0
     V(i,j,vy) = 0.d0
     V(i,j,vz) = 0.d0

     U(i,j,bx) = tanh(y(j))
     U(i,j,by) = 0.d0
     U(i,j,bz) = 0.d0
     U(i,j,ps) = 0.d0

! -- initial perturbation ---
     b1 = 0.03d0
     r2 = x(i)**2 + y(j)**2
     U(i,j,bx) = U(i,j,bx) - b1 * y(j) * exp(-r2/4.d0)
     U(i,j,by) = U(i,j,by) + b1 * x(i) * exp(-r2/4.d0)
! -- initial perturbation ---

     f1 = 1.d0 / ( gamma - 1 )
     v2 = dot_product( V(i,j,vx:vz), V(i,j,vx:vz) )
     B2 = dot_product( U(i,j,bx:bz), U(i,j,bx:bz) )
     U(i,j,mx:mz) = U(i,j,ro) * V(i,j,vx:vz)
     U(i,j,en)    = 0.5d0 * ( U(i,j,ro)*v2 + B2 ) + f1*V(i,j,pr)

  enddo
  enddo
  

  return
end subroutine modelp
