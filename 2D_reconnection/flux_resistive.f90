subroutine flux_resistive(F,U,VL,VR,EtS,dx,ix,jx,dir)
!-----------------------------------------------------------------------
!     Resistive flux terms in the X/Y directions
!-----------------------------------------------------------------------
!     2010/09/23  S. Zenitani  resistive HLL solver
!     2015/07/29  S. Zenitani  if-statements ==> max/min functions
!     2016/09/06  S. Zenitani  resistive part
!     2016/10/02  S. Zenitani  X/Y directions, loop jamming
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
! resistivity ** at the cell surface ** (EtS) [input]
  real(8), intent(in) :: EtS(ix,jx)
! direction [input]: 1 (X), 2 (Y), 3 (Z)
  integer, intent(in) :: dir
!-----------------------------------------------------------------------
! electric currents at the cell surface [local]
  real(8) :: JxS, JyS, JzS
  integer :: i, j, is, ie, js, je
  real(8) :: f1

  is = min(2,ix); ie = max(1,ix-1)
  js = min(2,jx); je = max(1,jx-1)

  select case(dir)
  case(1)
     is = 1; ie = ix-1
  case(2)
     js = 1; je = jx-1
  case(3)
  endselect

  f1 = 1.d0 / dx


  select case(dir)
!-----------------------------------------------------------------------
  case(1)

     do j=js,je
     do i=is,ie
! electric current at the X-surface (Toth+ 2008, JCP)
! B is taken from the cell center (U)
!       JxS = f1*0.25d0*( U(i,j+1,bz)+U(i+1,j+1,bz)-U(i,j-1,bz)-U(i+1,j-1,bz) )
        JyS = -f1*( U(i+1,j,bz)-U(i,j,bz) )
        JzS = f1*( ( U(i+1,j,by)-U(i,j,by) ) &
             - 0.25d0*( U(i,j+1,bx)+U(i+1,j+1,bx)-U(i,j-1,bx)-U(i+1,j-1,bx) ) )
! resistive flux
! B is taken from the left (VL) and the right states (VR)
        F(i,j,en) = F(i,j,en) + 0.5d0 * EtS(i,j) * &
             (JyS*(VL(i,j,bz)+VR(i,j,bz))-JzS*(VL(i,j,by)+VR(i,j,by)))
        F(i,j,by) = F(i,j,by) - EtS(i,j) * JzS
        F(i,j,bz) = F(i,j,bz) + EtS(i,j) * JyS
     enddo
     enddo

!-----------------------------------------------------------------------
  case(2)

     do j=js,je
     do i=is,ie
! electric current at the Y-surface (Toth+ 2008, JCP)
! B is taken from the cell center (U)
        JxS = f1*( U(i,j+1,bz)-U(i,j,bz) )
!       JyS = -f1*0.25d0*( U(i+1,j+1,bz)+U(i+1,j,bz)-U(i-1,j+1,bz)-U(i-1,j,bz) )
        JzS = f1*( 0.25d0*(U(i+1,j+1,by)+U(i+1,j,by)-U(i-1,j+1,by)-U(i-1,j,by)) &
             - ( U(i,j+1,bx)-U(i,j,bx) ) )
! resistive flux
! B is taken from the left (VL) and the right states (VR)
        F(i,j,en) = F(i,j,en) + 0.5d0 * EtS(i,j) * &
             ( JzS*(VL(i,j,bx)+VR(i,j,bx)) - JxS*(VL(i,j,bz)+VR(i,j,bz)) )
        F(i,j,bx) = F(i,j,bx) + EtS(i,j) * JzS
        F(i,j,bz) = F(i,j,bz) - EtS(i,j) * JxS
     enddo
     enddo

  endselect

  return
end subroutine flux_resistive
