subroutine model(U,x,y,dx,ix,jx)
  implicit none
  include 'param_rela.h'
  real(8) :: U(ix,jx,var1)
  real(8) :: V(var1)
  real(8) :: x(ix), y(jx), dx
  integer :: ix, jx
  integer :: i, j, izero
  real(8) :: u0, B2, bu
!-----------------------------------------------------------------------
! Initial configurations for a relativistic jet problem
!        Ref: S. Zenitani, M. Hesse, A. Klimas, ApJ, 712, 951 (2010)
!-----------------------------------------------------------------------
  real(8), parameter :: roL = 0.1d0
  real(8), parameter :: uxL = 0.d0, uyL = 0.d0, uzL = sqrt(48.d0)
  real(8), parameter :: prL = 10.d0, byL = 0.d0, bzL = 0.d0 ! model H1
!  real(8), parameter :: prL = 2.d0, byL = 28.d0, bzL = 0.d0 ! model M1
!  real(8), parameter :: prL = 2.d0, byL = 0.d0, bzL = 4.d0 ! model M2
  real(8), parameter :: roR = 1.d0
  real(8), parameter :: uxR = 0.d0, uyR = 0.d0, uzR = 0.d0
  real(8), parameter :: prR = 1.d0, byR = 0.d0, bzR = 0.d0
!-----------------------------------------------------------------------

! grid
  dx=0.4d0/real(ix-2)
  izero=ix/2
  x(izero)=-dx/2
  do i=izero+1,ix
     x(i) = x(i-1)+dx
  enddo
  do i=izero-1,1,-1
     x(i) = x(i+1)-dx
  enddo

! initial condition
  write(6,*) '[Left]'
  write(6,999) roL, prL, uzL, byL, bzL
  write(6,*) '[Right]'
  write(6,999) roR, prR, uzR, byR, bzR
999 format ('ro:', 1p, e10.3, ' pr:', 1p, e10.3, ' uz:', 1p, e10.3, &
         ' by:', 1p, e10.3, ' bz:', 1p, e10.3 )

  j = 1
  y(j) = 0.d0
  do i=1,ix
!    left
     if ( x(i) .lt. 0.d0 ) then
        V(ro) = roL
        V(pr) = prL
        V(ux) = uxL
        V(uy) = uyL
        V(uz) = uzL
        V(by) = byL
        V(bz) = bzL
!     right
     else
        V(ro) = roR
        V(pr) = prR
        V(ux) = uxR
        V(uy) = uyR
        V(uz) = uzR
        V(by) = byR
        V(bz) = bzR
     endif
     V(bx) = 0.d0

     B2 = dot_product( V(bx:bz), V(bx:bz) )
     bu = dot_product( V(bx:bz), V(ux:uz) )
     u0 = sqrt(1.d0 + dot_product( V(ux:uz), V(ux:uz) ))

     U(i,j,mx) = u0 * ( V(ro) + 4.d0 * V(pr) ) * V(ux) + (B2*V(ux) - bu*V(bx))/u0
     U(i,j,my) = u0 * ( V(ro) + 4.d0 * V(pr) ) * V(uy) + (B2*V(uy) - bu*V(by))/u0
     U(i,j,mz) = u0 * ( V(ro) + 4.d0 * V(pr) ) * V(uz) + (B2*V(uz) - bu*V(bz))/u0
     U(i,j,en) = u0*u0*(V(ro)+4*V(pr))-V(pr) + B2 - 0.5d0*(B2+bu**2)/(u0**2)
     U(i,j,de) = u0 * V(ro)
     U(i,j,bx) = V(bx)
     U(i,j,by) = V(by)
     U(i,j,bz) = V(bz)
     U(i,j,ps) = 0.d0
  enddo
     
  return
end subroutine model
