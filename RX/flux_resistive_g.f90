subroutine flux_resistive_g(G,U,VL,VR,EtS,dx,ix,jx)
!-----------------------------------------------------------------------
!     Resistive HLL solver in the Y direction
!-----------------------------------------------------------------------
!     2010/09/23  S. Zenitani  resistive HLL solver
!     2010/11/30  S. Zenitani  fixed bug in J_x
!     2015/07/29  S. Zenitani  if-statements ==> max/min functions
!     2016/09/06  S. Zenitani  resistive part
!-----------------------------------------------------------------------
  implicit none
  include 'param.h'
!-----------------------------------------------------------------------
  integer, intent(in) :: ix, jx
  real(8), intent(in) :: dx
! numerical flux (G) [inout]
  real(8), intent(inout) :: G(ix,jx,var1)
! left/right states (VL,VR) [input]
  real(8), intent(in) :: VL(ix,jx,var1), VR(ix,jx,var1)
! conserved variables (U) [input]
  real(8), intent(in) :: U(ix,jx,var1)
! resistivity at the cell surface (EtS) [input]
  real(8), intent(in) :: EtS(ix,jx)
!-----------------------------------------------------------------------
! electric current at the cell surface (these J's are local)
  real(8) :: JxS(ix,jx), JzS(ix,jx)
  integer :: i, j, is, ie, js, je
  real(8) :: f1

  is = min(2,jx); ie = max(1,ix-1)
  js = 1        ; je = jx-1

! surface current (Toth+ 2008, JCP)
  JxS = 0.d0
  JzS = 0.d0
  f1 = 1.d0 / dx
  do j=js,je
  do i=is,ie
     JxS(i,j) = f1*( U(i,j+1,bz)-U(i,j,bz) )
!     JyS(i,j) = -f1*0.25d0*( U(i+1,j+1,bz)+U(i+1,j,bz)-U(i-1,j+1,bz)-U(i-1,j,bz) )
     JzS(i,j) = f1*( 0.25d0*( U(i+1,j+1,by)+U(i+1,j,by)-U(i-1,j+1,by)-U(i-1,j,by) ) &
          - ( U(i,j+1,bx)-U(i,j,bx) ) )
  enddo
  enddo

! resistive fix to G
! Caution: J is surface value
!          B is taken from the left (VL) and the right states (VR)
  do j=js,je
  do i=is,ie
     G(i,j,en) = G(i,j,en) + 0.5d0 * EtS(i,j) * &
          ( JzS(i,j)*(VL(i,j,bx)+VR(i,j,bx)) - JxS(i,j)*(VL(i,j,bz)+VR(i,j,bz)) )
     G(i,j,bx) = G(i,j,bx) + EtS(i,j) * JzS(i,j)
     G(i,j,bz) = G(i,j,bz) - EtS(i,j) * JxS(i,j)
  enddo
  enddo
!-----------------------------------------------------------------------

  return
end subroutine flux_resistive_g
