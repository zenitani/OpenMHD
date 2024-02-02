function eta(xx,yy)
  real(8), intent(in) :: xx, yy
  real(8) :: eta, eta0, eta01
!-----------------------------------------------------------------------
! Resistivity
  real(8), parameter :: Rm1 = 60.d0   ! localized - same as Rm1 in set_dt_resistive
  real(8), parameter :: Rm0 = 1000.d0 ! background
!-----------------------------------------------------------------------
  eta0  = 1.d0 / Rm0
  eta01 = ( Rm0 - Rm1 ) / ( Rm0 * Rm1 )
!  eta = eta0 + eta01 * ( cosh( sqrt( xx**2+yy**2 )/2.d0 ) )**(-2)
  eta = eta0 + eta01 * ( cosh(min( sqrt( xx**2+yy**2 ), 25.d0 )) )**(-2)
! Here we use min(..., 25) to avoid an overflow problem on JSS/JAXA.
! Note that (cosh(25))**(-2) = 7.7E-22 is much smaller than
! the machine epsilon 10**(-16) in the 8-bit floating-point format.
end function eta

subroutine set_dt_resistive(dt,dx,cfl)
  real(8), intent(inout) :: dt
  real(8), intent(in) :: dx, cfl
!-----------------------------------------------------------------------
  real(8), parameter :: Rm1 = 60.d0   ! localized - same as Rm1 in eta function
!-----------------------------------------------------------------------
  dt = min( dt, 0.5d0*cfl*Rm1*(dx**2) )
end subroutine set_dt_resistive
  
subroutine flux_resistive(F,U,x,y,dx,ix,jx,dir)
!-----------------------------------------------------------------------
!     Resistive flux terms in the X/Y directions
!-----------------------------------------------------------------------
!     2010/09/23  S. Zenitani  resistive HLL solver
!     2016/09/06  S. Zenitani  resistive part
!     2016/10/02  S. Zenitani  X/Y directions, loop jamming
!     2024/02/02  S. Zenitani  eta(X,Y) function
!-----------------------------------------------------------------------
  implicit none
  include 'param.h'
!-----------------------------------------------------------------------
  integer, intent(in) :: ix, jx
! numerical flux (F) [inout]
  real(8), intent(inout) :: F(ix,jx,var1)
! conserved variables (U) [input]
  real(8), intent(in) :: U(ix,jx,var1)
! spatial location [input]
  real(8), intent(in) :: x(ix), y(jx), dx
! direction [input]: 1 (X), 2 (Y), 3 (Z)
  integer, intent(in) :: dir
!-----------------------------------------------------------------------
! electric currents at the cell surface [local]
  real(8) :: JxS, JyS, JzS
  integer :: i, j, is, ie, js, je
  real(8) :: eta, eta_dx, odx

  select case(dir)
  case(1)
     is = 1; ie = ix-1
     js = 2; je = jx-1
  case(2)
     is = 2; ie = ix-1
     js = 1; je = jx-1
  case(3)
  end select

  odx = 1.d0 / dx

  select case(dir)
!-----------------------------------------------------------------------
  case(1)

!$omp parallel do private(i,j,eta_dx,JxS,JyS,JzS)
     do j=js,je
     do i=is,ie
! electric current at the X-surface (Toth+ 2008, JCP)
! B is taken from the cell center (U)
        eta_dx = eta(x(i)+0.5d0*dx,y(j)) * odx
!       JxS = 0.25d0*eta_dx*( U(i,j+1,bz)+U(i+1,j+1,bz)-U(i,j-1,bz)-U(i+1,j-1,bz) )
        JyS = -eta_dx*( U(i+1,j,bz)-U(i,j,bz) )
        JzS =  eta_dx*( ( U(i+1,j,by)-U(i,j,by) ) &
             - 0.25d0*( U(i,j+1,bx)+U(i+1,j+1,bx)-U(i,j-1,bx)-U(i+1,j-1,bx) ) )
! resistive flux
        F(i,j,en) = F(i,j,en) + 0.5d0 * &
             (JyS*(U(i,j,bz)+U(i+1,j,bz))-JzS*(U(i,j,by)+U(i+1,j,by)))
!             (JyS*(VL(i,j,bz)+VR(i,j,bz))-JzS*(VL(i,j,by)+VR(i,j,by)))
        F(i,j,by) = F(i,j,by) - JzS
        F(i,j,bz) = F(i,j,bz) + JyS
     enddo
     enddo
!$omp end parallel do

!-----------------------------------------------------------------------
  case(2)

!$omp parallel do private(i,j,eta_dx,JxS,JyS,JzS)
     do j=js,je
     do i=is,ie
! electric current at the Y-surface (Toth+ 2008, JCP)
! B is taken from the cell center (U)
        eta_dx = eta(x(i),y(j)+0.5d0*dx) * odx
        JxS = eta_dx*( U(i,j+1,bz)-U(i,j,bz) )
!       JyS = -0.25d0*eta_dx*( U(i+1,j+1,bz)+U(i+1,j,bz)-U(i-1,j+1,bz)-U(i-1,j,bz) )
        JzS = eta_dx*( 0.25d0*(U(i+1,j+1,by)+U(i+1,j,by)-U(i-1,j+1,by)-U(i-1,j,by)) &
             - ( U(i,j+1,bx)-U(i,j,bx) ) )
! resistive flux
        F(i,j,en) = F(i,j,en) + 0.5d0 * &
             ( JzS*(U(i,j,bx)+U(i,j+1,bx)) - JxS*(U(i,j,bz)+U(i,j+1,bz)) )
!             ( JzS*(VL(i,j,bx)+VR(i,j,bx)) - JxS*(VL(i,j,bz)+VR(i,j,bz)) )
        F(i,j,bx) = F(i,j,bx) + JzS
        F(i,j,bz) = F(i,j,bz) - JxS
     enddo
     enddo
!$omp end parallel do

  end select

  return
end subroutine flux_resistive
