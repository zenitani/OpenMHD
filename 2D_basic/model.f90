subroutine model(U,V,x,y,dx,ix,jx)
  implicit none
  include 'param.h'
  real(8), intent(out) :: U(ix,jx,var1)
  real(8), intent(out) :: V(ix,jx,var2)
  real(8), intent(out) :: x(ix), y(jx), dx
  integer, intent(in)  :: ix, jx
! ---------------------------------------------------
! model parameters
  real(8), parameter :: pi = 4.d0*atan(1.d0)
! x & y positions (Note: domain_y(2) is automatically calculated)
  real(8), parameter :: domain_x(2) = (/0.d0, 2*pi/)
  real(8), parameter :: domain_y(1) = (/0.d0/)
! ---------------------------------------------------
  integer :: i, j
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
end subroutine model
