subroutine modelp(U,V,x,y,z,dx,ix,jx,kx)
!-----------------------------------------------------------------------
!     Initial configuration for the parallel version (mainp2)
!-----------------------------------------------------------------------
  use parallel
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
  integer :: iix, jjx, kkx
  real(8) :: B2, v2, f1, r2, b1
  real(8) :: tmpx(cart3d%sizes(1)*(ix-2) + 2)
  real(8) :: tmpy(cart3d%sizes(2)*(jx-2) + 2)
  real(8) :: tmpz(cart3d%sizes(3)*(kx-2) + 2)
! ---------------------------------------------------

! grid in the X direction
  iix = cart3d%sizes(1)*(ix-2) + 2
  dx = ( domain_x(2) - domain_x(1) ) / dble( iix-2 )
  tmpx(1)   = domain_x(1) - dx/2
! tmpx(iix) = domain_x(2) + dx/2
  do i=2,iix
     tmpx(i) = tmpx(i-1) + dx
  enddo
  do i=1,ix
     x(i) = tmpx(cart3d%coords(1)*(ix-2) + i)
  enddo

! grid in the Y direction
  jjx = cart3d%sizes(2)*(jx-2) + 2
  tmpy(1) = domain_y(1) - dx/2
  do j=2,jjx
     tmpy(j) = tmpy(j-1) + dx
  enddo
  do j=1,jx
     y(j) = tmpy(cart3d%coords(2)*(jx-2) + j)
  enddo

! grid in the Z direction
  kkx = cart3d%sizes(3)*(kx-2) + 2
  tmpz(1) = domain_z(1) - dx/2
  do k=2,kkx
     tmpz(k) = tmpz(k-1) + dx
  enddo
  do k=1,kx
     z(k) = tmpz(cart3d%coords(3)*(kx-2) + k)
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
end subroutine modelp
