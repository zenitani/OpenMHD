!-----------------------------------------------------------------------
!    Runge-Kutta routines
!-----------------------------------------------------------------------
! This file contains TVD-RK2 routines and standard RK2 routines
! Note that these routines assume dx == dy.


! 1/2 step of TVD Runge=Kutta method
subroutine rk_tvd21(wk,wk1,wF,wG,dt,dx,ix,jx)
  implicit none
  include 'param.h'
!-----------------------------------------------------------------------
  integer, intent(in) :: ix, jx
  real(8), intent(in) :: wk(ix,jx,var1)
  real(8), intent(out):: wk1(ix,jx,var1)
  real(8), intent(in) :: wF(ix,jx,var1), wG(ix,jx,var1)
  real(8), intent(in) :: dt, dx
!-----------------------------------------------------------------------
  integer :: i, j, k
  real(8) :: dtx

  dtx = dt/dx
!  dtx = dt
!  write(6,*) dt, dx, dtx
! Here, we employ collapse(2) in order to parallelize the outermost loop +1
!$omp parallel do private(i,j,k) collapse(2)
  do k=1,var1
     do j=2,jx-1
        do i=2,ix-1
           wk1(i,j,k) = wk(i,j,k) + &
                dtx * ( wF(i-1,j,k) - wF(i,j,k) + wG(i,j-1,k) - wG(i,j,k) )
        enddo
     enddo
  enddo
!$omp end parallel do

  return
end subroutine rk_tvd21


! 2/2 step of TVD Runge=Kutta method
subroutine rk_tvd22(wk,wk1,wF,wG,dt,dx,ix,jx)
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
!$omp parallel do private(i,j,k) collapse(2)
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
end subroutine rk_tvd22


! 1/2 step of the standard Runge=Kutta method
subroutine rk_std21(wk,wk1,wF,wG,dt,dx,ix,jx)
  implicit none
  include 'param.h'
!-----------------------------------------------------------------------
  integer, intent(in) :: ix, jx
  real(8), intent(in) :: wk(ix,jx,var1)
  real(8), intent(out):: wk1(ix,jx,var1)
  real(8), intent(in) :: wF(ix,jx,var1), wG(ix,jx,var1)
  real(8), intent(in) :: dt, dx
!-----------------------------------------------------------------------
  integer :: i, j, k
  real(8) :: dtx

  dtx = 0.5d0*dt/dx
!  dtx = dt
!  write(6,*) dt, dx, dtx
! Here, we employ collapse(2) in order to parallelize the outermost loop +1
!$omp parallel do private(i,j,k) collapse(2)
  do k=1,var1
     do j=2,jx-1
        do i=2,ix-1
           wk1(i,j,k) = wk(i,j,k) + &
                dtx * ( wF(i-1,j,k) - wF(i,j,k) + wG(i,j-1,k) - wG(i,j,k) )
        enddo
     enddo
  enddo
!$omp end parallel do

  return
end subroutine rk_std21


! 2/2 step of the standard Runge=Kutta method
subroutine rk_std22(wk,wF,wG,dt,dx,ix,jx)
  implicit none
  include 'param.h'
!-----------------------------------------------------------------------
  integer, intent(in) :: ix, jx
  real(8), intent(inout) :: wk(ix,jx,var1)
  real(8), intent(in) :: wF(ix,jx,var1), wG(ix,jx,var1)
  real(8), intent(in) :: dt, dx
!-----------------------------------------------------------------------
  integer :: i, j, k
  real(8) :: dtx

  dtx = dt/dx
!$omp parallel do private(i,j,k) collapse(2)
  do k=1,var1
     do j=2,jx-1
        do i=2,ix-1
           wk(i,j,k) = wk(i,j,k) + &
                dtx * ( wF(i-1,j,k) - wF(i,j,k) + wG(i,j-1,k) - wG(i,j,k) )
        enddo
     enddo
  enddo
!$omp end parallel do
  
  return
end subroutine rk_std22

