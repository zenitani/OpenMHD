subroutine step1(wk,wk1,wF,wG,dt,dx,ix,jx)
!-----------------------------------------------------------------------
!     1/2 step of the standard Runge-Kutta method
!-----------------------------------------------------------------------
  implicit none
  include 'param.h'
!-----------------------------------------------------------------------
  integer, intent(in) :: ix, jx
  real(8) :: wk(ix,jx,var1), wk1(ix,jx,var1)
  real(8), intent(in) :: wF(ix,jx,var1), wG(ix,jx,var1)
  real(8), intent(in) :: dt, dx
!-----------------------------------------------------------------------
  integer :: i, j, k
  real(8) :: dtx

  dtx = 0.5d0*dt/dx
!  dtx = dt
!  write(6,*) dt, dx, dtx
  do k=1,var1
     do j=2,jx-1
        do i=2,ix-1
           wk1(i,j,k) = wk(i,j,k) + &
                dtx * ( wF(i-1,j,k) - wF(i,j,k) + wG(i,j-1,k) - wG(i,j,k) )
        enddo
     enddo
  enddo

  return
end subroutine step1
