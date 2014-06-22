subroutine modelp(U,V,x,y,dx,ix,jx,myrank,npe)
  implicit none
  include 'param.h'
  real(8) :: U(ix,jx,var1)
  real(8) :: V(ix,jx,var2)
  real(8) :: x(ix), y(jx), dx
  integer :: ix, jx
  integer :: myrank, npe    ! for MPI
  real(8) :: tmpx(npe*(ix-2)+2)  ! for MPI
! ---------------------------------------------------
  real(8), parameter :: Lx    = 16.d0
  real(8), parameter :: alpha = 0.5d0
! ---------------------------------------------------
  real(8), parameter :: pi2   = 8.d0 * atan(1.d0)
  integer :: i, j, izero, jzero
  real(8) :: B2, v2, f1

! ---------------------------------------------------
! grid
  dx = Lx / dble( npe*(ix-2) )
  izero = 1
!  izero = (npe*(ix-2))/2+1
  tmpx(izero) = -dx/2
  do i=izero+1,npe*(ix-2)+2
     tmpx(i) = tmpx(i-1) + dx
  enddo
  do i=izero-1,1,-1
     tmpx(i) = tmpx(i+1) - dx
  enddo
  do i=1,ix
     x(i) = tmpx(myrank*(ix-2) + i)
  enddo

  jzero = jx/2
  y(jzero) = -dx/2
  do j=jzero+1,jx
     y(j) = y(j-1)+dx
  enddo
  do j=jzero-1,1,-1
     y(j) = y(j+1)-dx
  enddo
! ---------------------------------------------------

  do j=1,jx
  do i=1,ix

     U(i,j,ro) = 0.5d0 * ( (1.d0+alpha) + (1.d0-alpha)*tanh(y(j)) )
     V(i,j,pr) = 0.15d0
     V(i,j,vx) = -0.5d0 * tanh(y(j))
     V(i,j,vy) = 0.1d0 * sin( pi2*x(i) / Lx ) / cosh(y(j))**2
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
end subroutine modelp
