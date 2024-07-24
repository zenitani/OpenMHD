subroutine model(U,V,x,y,z,dx,ix,jx,kx)
!-----------------------------------------------------------------------
!     Initial configuration for the serial code
!-----------------------------------------------------------------------
  implicit none
  include 'param.h'
  real(8), intent(out) :: U(ix,jx,kx,var1)
  real(8), intent(out) :: V(ix,jx,kx,var2)
  real(8), intent(out) :: x(ix), y(jx), z(kx), dx
  integer, intent(in)  :: ix, jx, kx
! ---------------------------------------------------
! x & y positions (Note: domain_y(2) is automatically calculated)
  real(8), parameter :: domain_x(2) = (/0.d0, 160.d0/)
  real(8), parameter :: domain_y(1) = (/0.d0/)
  real(8), parameter :: domain_z(1) = (/0.d0/)
! plasma beta in the upstream region
  real(8), parameter :: beta = 0.2d0
! ---------------------------------------------------
  integer :: i, j, k
  real(8) :: B2, v2, f1, r2, b1
! ---------------------------------------------------

  dx = ( domain_x(2) - domain_x(1) ) / dble( ix-2 )
  x(1)  = domain_x(1) - dx/2
! x(ix) = domain_x(2) + dx/2
  do i=2,ix
     x(i) = x(i-1) + dx
  enddo
  y(1) = domain_y(1) - dx/2
  do j=2,jx
     y(j) = y(j-1) + dx
  enddo
  z(1) = domain_z(1) - dx/2
  do k=2,kx
     z(k) = z(k-1) + dx
  enddo
! ---------------------------------------------------

  do k=1,kx
  do j=1,jx
  do i=1,ix

! Here we use min(..., 25) to avoid NaN on some systems.
! Note that (cosh(25))**(-2) = 7.7E-22 is extremely small.
!     U(i,j,k,ro) = 1.d0 + cosh(min( y(j), 25.d0 ))**(-2) / beta
     U(i,j,k,ro) = 1.d0 + cosh(y(j))**(-2) / beta
     V(i,j,k,pr) = 0.5d0 * beta * U(i,j,k,ro)
     V(i,j,k,vx) = 0.d0
     V(i,j,k,vy) = 0.d0
     V(i,j,k,vz) = 0.d0

     U(i,j,k,bx) = tanh(y(j))
     U(i,j,k,by) = 0.d0
     U(i,j,k,bz) = 0.d0
     U(i,j,k,ps) = 0.d0

! -- initial perturbation ---
     b1 = 0.03d0
     r2 = x(i)**2 + y(j)**2
     U(i,j,k,bx) = U(i,j,k,bx) - b1 * y(j) * exp(-r2/4.d0)
     U(i,j,k,by) = U(i,j,k,by) + b1 * x(i) * exp(-r2/4.d0)
! -- initial perturbation ---

     f1 = 1.d0 / ( gamma - 1 )
     v2 = dot_product( V(i,j,k,vx:vz), V(i,j,k,vx:vz) )
     B2 = dot_product( U(i,j,k,bx:bz), U(i,j,k,bx:bz) )
     U(i,j,k,mx:mz) = U(i,j,k,ro) * V(i,j,k,vx:vz)
     U(i,j,k,en)    = 0.5d0 * ( U(i,j,k,ro)*v2 + B2 ) + f1*V(i,j,k,pr)

  enddo
  enddo
  enddo
  

  return
end subroutine model
