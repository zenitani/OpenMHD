subroutine modelp(U,V,x,y,dx,ix,jx,myrank,npe)
  implicit none
  include 'param.h'
  real(8) :: U(ix,jx,var1)
  real(8) :: V(ix,jx,var2)
  real(8) :: x(ix), y(jx), dx
  integer :: ix, jx
  integer :: myrank, npe    ! for MPI
  real(8) :: tmpx(npe*(ix-2)+2)  ! for MPI
  integer :: i, j, izero, jzero
  real(8) :: B2, v2, f1

! ---------------------------------------------------
! grid in the X direction : L_x = 2 pi
  dx = ( 8.d0*atan(1.d0) ) / dble( npe*(ix-2) )
!  izero = (npe*(ix-2))/2+1
  izero = 1
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

!  jzero = jx/2
  jzero = 1
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
