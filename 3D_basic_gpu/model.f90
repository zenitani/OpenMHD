function mya3(xx,yy,zz)
  real(8), intent(in) :: xx, yy, zz
  real(8), parameter :: myr = 0.3d0
  real(8) :: mya3
  mya3 = max( 0.d0, myr - sqrt( yy**2 + min( (2*xx+zz-2)**2, (2*xx+zz)**2, (2*xx+zz+2)**2 )/5.d0 ) )
end function mya3

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
  real(8), parameter :: domain_x(2) = (/-0.5d0, 0.5d0/)
  real(8), parameter :: domain_y(1) = (/-0.5d0/)
  real(8), parameter :: domain_z(1) = (/-1.d0/)
! ---------------------------------------------------
  integer :: i, j, k
  real(8) :: B2, v2, f1
  real(8) :: fx, fz, dx2, xx, yy, zz
  real(8) :: mya3
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

  fx = - (0.001d0/dx) * (1.d0/sqrt(5.d0))
  fz = + (0.001d0/dx) * (2.d0/sqrt(5.d0))
  dx2 = dx/2

  do k=1,kx
  do j=1,jx
  do i=1,ix

     U(i,j,k,ro) = 1.d0
     V(i,j,k,pr) = 1.d0

     V(i,j,k,vx) = 1.d0
     V(i,j,k,vy) = 1.d0
     V(i,j,k,vz) = 2.d0

     xx = x(i); yy = y(j); zz = z(k)
     U(i,j,k,bx) =  fz * ( mya3(xx,yy+dx2,zz) - mya3(xx,yy-dx2,zz)  )
     U(i,j,k,by) =  fx * ( mya3(xx,yy,zz+dx2) - mya3(xx,yy,zz-dx2)  ) - fz * ( mya3(xx+dx2,yy,zz) - mya3(xx-dx2,yy,zz)  )
     U(i,j,k,bz) = -fx * ( mya3(xx,yy+dx2,zz) - mya3(xx,yy-dx2,zz)  )

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
