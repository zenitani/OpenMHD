subroutine glm_ss2(U,ch,dt,ix,jx)
!-----------------------------------------------------------------------
!     GLM source term solver
!-----------------------------------------------------------------------
! *** caution ***
! The third argument is changed from dt to (1/2)dt.
! Previously, we use
!   call glm_ss(U,ch,0.5d0*dt,ix,jx)
! but now we use
!   call glm_ss2(U,ch,dt,ix,jx)
!-----------------------------------------------------------------------
  implicit none
  include 'param.h'
!-----------------------------------------------------------------------
! The ratio of the hyperbolic and parabolic effects
!   c_r = ( c_p**2 / c_h ) ~ 0.18
! Ref: Dedner et al. (2002)  p.657, p.661
  real(8), parameter :: cr = 0.2d0
!-----------------------------------------------------------------------
  integer, intent(in) :: ix, jx
  real(8), intent(inout) :: U(ix,jx,var1)  ! conserved variables (U)
  real(8), intent(in) :: ch, dt  ! div cleaning wave speed and the timestep
!-----------------------------------------------------------------------
  real(8) :: f1 = 1.d0
  real(8) :: tmp(ix,jx)

  f1 = exp( -0.5d0*dt*ch/cr )
!$omp parallel
!$omp workshare
  tmp(:,:) = f1 * U(:,:,ps)
!$omp end workshare
!$omp barrier
!$omp workshare
  U(:,:,ps) = tmp(:,:)
!$omp end workshare
!$omp end parallel

  return
end subroutine glm_ss2
