subroutine glm_ss2(U,ch,dt,ix,jx,kx)
!-----------------------------------------------------------------------
!     GLM source term solver
!-----------------------------------------------------------------------
  implicit none
  include 'param.h'
!-----------------------------------------------------------------------
! The ratio of the hyperbolic and parabolic effects
!   c_r = ( c_p**2 / c_h ) ~ 0.18
! Ref: Dedner et al. (2002)  p.657, p.661
  real(8), parameter :: cr = 0.2d0
!-----------------------------------------------------------------------
  integer, intent(in) :: ix, jx, kx
  real(8), intent(inout) :: U(ix,jx,kx,var1)  ! conserved variables (U)
  real(8), intent(in) :: ch, dt  ! div cleaning wave speed and the timestep
!-----------------------------------------------------------------------
  real(8) :: f1 = 1.d0
  real(8) :: tmp(ix,jx,kx)

  f1 = exp( -0.5d0*dt*ch/cr )
!$omp parallel
!$omp workshare
  tmp(:,:,:) = f1 * U(:,:,:,ps)
!$omp end workshare
!$omp barrier
!$omp workshare
  U(:,:,:,ps) = tmp(:,:,:)
!$omp end workshare
!$omp end parallel

  return
end subroutine glm_ss2
