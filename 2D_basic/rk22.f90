subroutine rk22(wk,wk1,wF,wG,dt,dx,ix,jx)
!-----------------------------------------------------------------------
!     2/2 step of TVD Runge=Kutta method
!-----------------------------------------------------------------------
  implicit none
  include 'param.h'
!-----------------------------------------------------------------------
  integer, intent(in) :: ix, jx
  real(8) :: wk(ix,jx,var1)
  real(8), intent(in) :: wk1(ix,jx,var1)
  real(8), intent(in) :: wF(ix,jx,var1), wG(ix,jx,var1)
  real(8), intent(in) :: dt, dx
!-----------------------------------------------------------------------
  integer :: i, j, k
  real(8) :: dtx

  dtx = dt/dx
!$omp parallel do private(i,j,k)
  do k=1,var1
     do j=2,jx-1
        do i=2,ix-1
           wk(i,j,k) = 0.5d0*( wk(i,j,k)+wk1(i,j,k) + &
                dtx * ( wF(i-1,j,k) - wF(i,j,k) + wG(i,j-1,k) - wG(i,j,k) ) &
           )
        enddo
     enddo
  enddo
!$omp end parallel do

  return
end subroutine rk22
