subroutine glm_ss(U,ch,dt,ix,jx)
!-----------------------------------------------------------------------
!     GLM source term solver
!-----------------------------------------------------------------------
! *** caution ***
! A half timestep may be used for the 2nd order time-marching:
!   call glm_ss(U,ch,0.5d0*dt,ix,jx)
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

  f1 = exp( -dt*ch/cr )
!$omp parallel workshare
  U(:,:,ps) = f1 * U(:,:,ps)
!$omp end parallel workshare

  return
end subroutine glm_ss
