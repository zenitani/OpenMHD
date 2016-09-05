subroutine flux_resistive_f(F,U,VL,VR,EtS,dx,ix,jx)
!-----------------------------------------------------------------------
!     Resistive part in the X direction
!-----------------------------------------------------------------------
!     2010/09/23  S. Zenitani  resistive HLL solver
!     2015/07/29  S. Zenitani  if-statements ==> max/min functions
!     2016/09/06  S. Zenitani  resistive part
!-----------------------------------------------------------------------
  implicit none
  include 'param.h'
!-----------------------------------------------------------------------
  integer, intent(in) :: ix, jx
  real(8), intent(in) :: dx
! numerical flux (F) [inout]
  real(8), intent(inout) :: F(ix,jx,var1)
! left/right states (VL,VR) [input]
  real(8), intent(in) :: VL(ix,jx,var1), VR(ix,jx,var1)
! conserved variables (U) [input]
  real(8), intent(in) :: U(ix,jx,var1)
! resistivity at the cell surface (EtS) [input]
  real(8), intent(in) :: EtS(ix,jx)
!-----------------------------------------------------------------------
! electric current at the cell surface (these J's are local)
  real(8) :: JyS(ix,jx), JzS(ix,jx)
  integer :: i, j, is, ie, js, je
  real(8) :: f1

  is = 1        ; ie = ix-1
  js = min(2,jx); je = max(1,jx-1)

! surface current (Toth+ 2008, JCP)
  JyS = 0.d0
  JzS = 0.d0
  f1 = 1.d0 / dx
  do j=js,je
  do i=is,ie
!     JxS(i,j) = f1*0.25d0*( U(i,j+1,bz)+U(i+1,j+1,bz)-U(i,j-1,bz)-U(i+1,j-1,bz) )
     JyS(i,j) = -f1*( U(i+1,j,bz)-U(i,j,bz) )
     JzS(i,j) = f1*( ( U(i+1,j,by)-U(i,j,by) ) &
                     - 0.25d0*( U(i,j+1,bx)+U(i+1,j+1,bx)-U(i,j-1,bx)-U(i+1,j-1,bx) ) )
  enddo
  enddo

! resistive fix to F
! Caution: J is surface value
!          B is taken from the left (VL) and the right states (VR)
  do j=js,je
  do i=is,ie
     F(i,j,en) = F(i,j,en) + 0.5d0 * EtS(i,j) * &
          ( JyS(i,j)*(VL(i,j,bz)+VR(i,j,bz)) - JzS(i,j)*(VL(i,j,by)+VR(i,j,by)) )
     F(i,j,by) = F(i,j,by) - EtS(i,j) * JzS(i,j)
     F(i,j,bz) = F(i,j,bz) + EtS(i,j) * JyS(i,j)
  enddo
  enddo
!-----------------------------------------------------------------------

  return
end subroutine flux_resistive_f
