subroutine set_dt2(Rm,dt,dx,cfl)
!-----------------------------------------------------------------------
!     Timestep for the diffusion term (resistivity) : cfl < 1.0
!     For simplicity, we reuse the CFL factor.
!-----------------------------------------------------------------------
  implicit none
  real(8) :: Rm, dt, dx, cfl

  if( dt .gt. 0.5d0*cfl*Rm*(dx**2) ) then
     dt = 0.5d0*cfl*Rm*(dx**2)
  endif

  return
end subroutine set_dt2
