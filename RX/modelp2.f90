subroutine modelp2(U,V,x,y,dx,ix,jx,myrank,npe)
!-----------------------------------------------------------------------
!     Initial configuration for the parallel version (mainp2)
!-----------------------------------------------------------------------
  implicit none
  include 'param.h'
  real(8) :: U(ix,jx,var1)
  real(8) :: V(ix,jx,var2)
  real(8) :: x(ix), y(jx), dx
  integer :: ix, jx
  integer :: myrank, npe    ! for MPI
  real(8) :: tmpx(npe*(ix-2)+2)  ! for MPI
! ---------------------------------------------------
! system size in the X direction (L_x): L_y is automatically calculated
  real(8), parameter :: Lx   = 200.d0
! plasma beta in the upstream region
  real(8), parameter :: beta =   0.2d0
! ---------------------------------------------------
  integer :: i, j, izero, jzero
  real(8) :: B2, v2, f1, r2, b1

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

     U(i,j,ro) = 1.d0 + cosh(y(j))**(-2) / beta
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
end subroutine modelp2
