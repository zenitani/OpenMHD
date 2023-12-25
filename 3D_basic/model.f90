subroutine model(U,V,x,y,z,dx,ix,jx,kx)
  implicit none
  include 'param.h'
  real(8), intent(out) :: U(ix,jx,kx,var1)
  real(8), intent(out) :: V(ix,jx,kx,var2)
  real(8), intent(out) :: x(ix), y(jx), z(kx), dx
  integer, intent(in)  :: ix, jx, kx
! ---------------------------------------------------
! model parameters
  real(8), parameter :: pi = 4.d0*atan(1.d0)
! X, Y, Z ranges (Note: domain_y(2) and domain_z(2) are automatically calculated)
  real(8), parameter :: domain_x(2) = (/0.d0, 2*pi/)
  real(8), parameter :: domain_y(1) = (/0.d0/)
  real(8), parameter :: domain_z(1) = (/0.d0/)
! ---------------------------------------------------
  integer :: i, j, k
  real(8) :: B2, v2, f1
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

     U(i,j,k,ro) = gamma**2
     V(i,j,k,pr) = gamma

     V(i,j,k,vx) = -sin(y(j))
     V(i,j,k,vy) = sin(x(i))
     V(i,j,k,vz) = 0.d0
!     V(i,j,k,vx) = 0.d0
!     V(i,j,k,vy) = -sin(z(k))
!     V(i,j,k,vz) = sin(y(j))
!     V(i,j,k,vx) = sin(z(k))
!     V(i,j,k,vy) = 0.d0
!     V(i,j,k,vz) = -sin(x(i))

     U(i,j,k,bx) = -sin(y(j))
     U(i,j,k,by) = sin(2*x(i))
     U(i,j,k,bz) = 0.d0
!     U(i,j,k,bx) = 0.d0
!     U(i,j,k,by) = -sin(z(k))
!     U(i,j,k,bz) = sin(2*y(j))
!     U(i,j,k,by) = sin(2*z(k))
!     U(i,j,k,bz) = 0.d0
!     U(i,j,k,bx) = -sin(x(i))

     U(i,j,k,ps) = 0.d0

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
