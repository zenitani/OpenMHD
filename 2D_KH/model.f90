subroutine model(U,V,x,y,dx,ix,jx)
  implicit none
  include 'param.h'
  real(8), intent(out) :: U(ix,jx,var1)
  real(8), intent(out) :: V(ix,jx,var2)
  real(8), intent(out) :: x(ix), y(jx), dx
  integer, intent(in)  :: ix, jx
! ---------------------------------------------------
! model parameters
! x & y positions (Note: domain_y(2) is automatically calculated)
  real(8), parameter :: domain_x(2) = (/0.d0, 16.d0/)
  real(8), parameter :: domain_y(1) = (/-10.d0/)
  real(8), parameter :: pi2   = 8.d0 * atan(1.d0)
  real(8), parameter :: alpha = 0.2d0
! ---------------------------------------------------
  integer :: i, j
  real(8) :: B2, v2, f1, Lx
! ---------------------------------------------------

! grid in the X direction
  Lx = ( domain_x(2) - domain_x(1) )
  dx = Lx / dble( ix-2 )
  x(1)  = domain_x(1) - dx/2
! x(ix) = domain_x(2) + dx/2
  do i=2,ix
     x(i) = x(i-1) + dx
  enddo
  y(1) = domain_y(1) - dx/2
  do j=2,jx
     y(j) = y(j-1) + dx
  enddo
! ---------------------------------------------------

  do j=1,jx
  do i=1,ix

     U(i,j,ro) = 0.5d0 * ( (1.d0+alpha) + (1.d0-alpha)*tanh(y(j)) )
     V(i,j,pr) = 0.15d0
     V(i,j,vx) = -0.5d0 * tanh(y(j))
! Here we use min(..., 25) to avoid NaN on some systems.
! Note that (cosh(25))**(-2) = 7.7E-22 is extremely small.
     V(i,j,vy) = 0.1d0 * sin( pi2*x(i) / Lx ) / cosh(min( y(j), 25.d0 ))**2
!     V(i,j,vy) = 0.1d0 * sin( pi2*x(i) / Lx ) / cosh(y(j))**2
     V(i,j,vz) = 0.d0

     U(i,j,bx) = 0.d0
     U(i,j,by) = 0.d0
     U(i,j,bz) = 1.d0
     U(i,j,ps) = 0.d0

     f1 = 1.d0 / ( gamma - 1 )
     v2 = dot_product( V(i,j,vx:vz), V(i,j,vx:vz) )
     B2 = dot_product( U(i,j,bx:bz), U(i,j,bx:bz) )
     U(i,j,mx:mz) = U(i,j,ro) * V(i,j,vx:vz)
     U(i,j,en)    = 0.5d0 * ( U(i,j,ro)*v2 + B2 ) + f1*V(i,j,pr)

  enddo
  enddo

  return
end subroutine model
