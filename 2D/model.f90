subroutine model(U,V,x,y,dx,ix,jx)
  implicit none
  include 'param.h'
  integer, intent(in) :: ix, jx
  real(8) :: U(ix,jx,var1)
  real(8) :: V(ix,jx,var2)
  real(8) :: x(ix), y(jx), dx
  integer :: i, j, izero, jzero
  real(8) :: B2, v2, f1

! grid in the X direction : L_x = 2 pi
!  dx = 1.d0/dble(ix-2)
  dx = ( 8.d0*atan(1.d0) ) / dble(ix-2)
!  izero=ix/2
  izero=1
  x(izero)=-dx/2
  do i=izero+1,ix
     x(i) = x(i-1)+dx
  enddo
  do i=izero-1,1,-1
     x(i) = x(i+1)-dx
  enddo

!  jzero = jx/2
  jzero = 1
  y(jzero) = -dx/2
  do j=jzero+1,jx
     y(j) = y(j-1)+dx
  enddo
  do j=jzero-1,1,-1
     y(j) = y(j+1)-dx
  enddo

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
