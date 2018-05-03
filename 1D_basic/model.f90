subroutine model(U,V,x,y,dx,ix,jx)
  implicit none
  include 'param.h'
  real(8) :: U(ix,jx,var1)
  real(8) :: V(ix,jx,var2)
  real(8) :: x(ix), y(jx), dx
  integer :: ix, jx
  integer :: i, j
  real(8) :: ro0, vx0, vy0, vz0, pr0, bx0, by0, bz0
  real(8) :: B2, v2, f1
! ---------------------------------------------------
! x locations
  real(8), parameter :: domain_x(2) = (/-0.5d0, 0.5d0/)
! ---------------------------------------------------

! 1D in the x direction
  dx = ( domain_x(2) - domain_x(1) ) / dble( ix-2 )
  x(1)  = domain_x(1) - dx/2
! x(ix) = domain_x(2) + dx/2
  do i=2,ix
     x(i) = x(i-1) + dx
  enddo
  y(1) = 0.d0

  do j=1,jx
  do i=1,ix
! -------- left --------------------
     if ( x(i) < 0.d0 ) then
! Dai & Woodward (1994)
        ro0 = 1.08d0
        pr0 = 0.95d0
        vx0 = 1.2d0
        vy0 = 0.01d0
        vz0 = 0.5d0
        bx0 = 2.d0 / sqrt(16*atan(1.d0))
        by0 = 3.6d0 / sqrt(16*atan(1.d0))
        bz0 = 2.d0 / sqrt(16*atan(1.d0))
! Brio & Wu (1988)
        ro0 = 1.d0
        pr0 = 1.d0
        vx0 = 0.d0
        vy0 = 0.d0
        vz0 = 0.d0
        bx0 = 0.75d0
        by0 = 1.d0
        bz0 = 0.d0
! Falle et al. (1998) - slow shock
!        ro0 = 1.368d0
!        pr0 = 1.769d0
!        vx0 = 0.269d0
!        vy0 = 1.d0
!        vz0 = 0.d0
!        bx0 = 1.d0
!        by0 = 0.d0
!        bz0 = 0.d0
! Falle et al. (1998) - rarefaction
!        ro0 = 1.d0
!        pr0 = 2.d0
!        vx0 = 0.d0
!        vy0 = 0.d0
!        vz0 = 0.d0
!        bx0 = 1.d0
!        by0 = 0.d0
!        bz0 = 0.d0
! -------- right --------------------
     else
! Dai & Woodward (1994)
        ro0 = 1.d0
        pr0 = 1.d0
        vx0 = 0.d0
        vy0 = 0.d0
        vz0 = 0.d0
        bx0 = 2.d0 / sqrt(16*atan(1.d0))
        by0 = 4.d0 / sqrt(16*atan(1.d0))
        bz0 = 2.d0 / sqrt(16*atan(1.d0))
! Brio & Wu (1988)
        ro0 = 0.125d0
        pr0 = 0.1d0
        vx0 = 0.d0
        vy0 = 0.d0
        vz0 = 0.d0
        bx0 = 0.75d0
        by0 =-1.d0
        bz0 = 0.d0
! Falle et al. (1998) - slow shock
!        ro0 = 1.d0
!        pr0 = 1.d0
!        vx0 = 0.d0
!        vy0 = 0.d0
!        vz0 = 0.d0
!        bx0 = 1.d0
!        by0 = 1.d0
!        bz0 = 0.d0
! Falle et al. (1998) - rarefaction
!        ro0 = 0.2d0
!        pr0 = 0.1368d0
!        vx0 = 1.186d0
!        vy0 = 2.967d0
!        vz0 = 0.d0
!        bx0 = 1.d0
!        by0 = 1.6405d0
!        bz0 = 0.d0
     endif

     V(i,j,vx) = vx0
     V(i,j,vy) = vy0
     V(i,j,vz) = vz0
     V(i,j,pr) = pr0
     U(i,j,ro) = ro0
     U(i,j,bx) = bx0
     U(i,j,by) = by0
     U(i,j,bz) = bz0
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
