program main
!-----------------------------------------------------------------------
!     OpenMHD  Reconnection solver (serial version)
!-----------------------------------------------------------------------
!     2010/09/25  S. Zenitani  HLL reconnection code
!-----------------------------------------------------------------------
  implicit none
  include 'param.h'
  integer, parameter :: ix = 2000 + 2
  integer, parameter :: jx =  500 + 2
  integer, parameter :: loop_max = 1000000
  real(8), parameter :: tend  = 251.0d0
  real(8), parameter :: dtout =  25.d0 ! output interval
  real(8), parameter :: cfl   = 0.35d0 ! time step
! Slope limiter  (0: flat, 1: minmod, 2: MC, 3: van Leer, 4: Koren)
  integer, parameter :: lm_type   = 1
! Numerical flux (0: LLF, 1: HLL, 2: HLLC, 3: HLLD)
  integer, parameter :: flux_type = 3
! Time marching  (0: TVD RK2, 1: RK2)
  integer, parameter :: time_type = 0
!-----------------------------------------------------------------------
! See also model.f90
!-----------------------------------------------------------------------
  integer :: n_loop,n_output
  real(8) :: t, dt, t_output
  real(8) :: ch
  character*256 :: filename
!-----------------------------------------------------------------------
  real(8) :: x(ix), y(jx), dx
  real(8) :: U(ix,jx,var1)  ! conserved variables (U)
  real(8) :: U1(ix,jx,var1) ! conserved variables: medium state (U*)
  real(8) :: V(ix,jx,var2)  ! primitive variables (V)
  real(8) :: VL(ix,jx,var1), VR(ix,jx,var1) ! interpolated states
  real(8) :: F(ix,jx,var1), G(ix,jx,var1)   ! numerical flux (F,G)
!-----------------------------------------------------------------------

  t    =  0.d0
  dt   =  0.d0
  call model(U,V,x,y,dx,ix,jx)
  call bc_for_U(U,ix,jx)
  call set_dt(U,V,ch,dt,dx,cfl,ix,jx)
  call set_dt_resistive(dt,dx,cfl)
  t_output = -dt/3.d0
  n_output =  0

  if ( dt > dtout ) then
     write(6,*) 'error: ', dt, '>', dtout
     stop
  endif
  write(6,*) '[Params]'
  write(6,998) dt, dtout, ix, jx
  write(6,999) lm_type, flux_type, time_type
998 format (' dt:', 1p, e10.3, ' dtout:', 1p, e10.3, ' grids:', i5, i5 )
999 format (' limiter: ', i1, '  flux: ', i1, '  time-marching: ', i1 )
  write(6,*) '== start =='

!-----------------------------------------------------------------------
  do n_loop=1,loop_max

     write(6,*) ' t = ', t
!    Recovering primitive variables
!     write(6,*) 'U --> V'
     call u2v(U,V,ix,jx)
!   -----------------
!    [ output ]
     if ( t >= t_output ) then
        write(6,*) 'data output   t = ', t
        write(filename,990) n_output
990     format ('data/field-',i5.5,'.dat')
        call fileio_output(filename,ix,jx,t,x,y,U,V)
        n_output = n_output + 1
        t_output = t_output + dtout
     endif
!    [ end? ]
     if ( t >= tend )  exit
     if ( n_loop >= loop_max ) then
        write(6,*) 'max loop'
        exit
     endif
!   -----------------
!    CFL condition
     call set_dt(U,V,ch,dt,dx,cfl,ix,jx)
     call set_dt_resistive(dt,dx,cfl)
!    GLM solver for the first half timestep
!    This should be done after set_dt()
     call glm_ss2(U,ch,dt,ix,jx)

!    Slope limiters on primitive variables
!     write(6,*) 'V --> VL, VR (F)'
     call limiter(V(:,:,vx),VL(:,:,vx),VR(:,:,vx),ix,jx,1,lm_type)
     call limiter(V(:,:,vy),VL(:,:,vy),VR(:,:,vy),ix,jx,1,lm_type)
     call limiter(V(:,:,vz),VL(:,:,vz),VR(:,:,vz),ix,jx,1,lm_type)
     call limiter(V(:,:,pr),VL(:,:,pr),VR(:,:,pr),ix,jx,1,lm_type)
     call limiter(U(:,:,ro),VL(:,:,ro),VR(:,:,ro),ix,jx,1,lm_type)
     call limiter(U(:,:,bx),VL(:,:,bx),VR(:,:,bx),ix,jx,1,lm_type)
     call limiter(U(:,:,by),VL(:,:,by),VR(:,:,by),ix,jx,1,lm_type)
     call limiter(U(:,:,bz),VL(:,:,bz),VR(:,:,bz),ix,jx,1,lm_type)
     call limiter(U(:,:,ps),VL(:,:,ps),VR(:,:,ps),ix,jx,1,lm_type)
!    fix flux bc (F)
     call bc_for_F(VL,VR,ix,jx)
!    Numerical flux in the X direction (F)
!     write(6,*) 'VL, VR --> F'
     call flux_solver(F,VL,VR,ch,ix,jx,1,flux_type)
     call flux_resistive(F,U,x,y,dx,ix,jx,1)

!    Slope limiters on primitive variables
!     write(6,*) 'V --> VL, VR (G)'
     call limiter(V(:,:,vx),VL(:,:,vx),VR(:,:,vx),ix,jx,2,lm_type)
     call limiter(V(:,:,vy),VL(:,:,vy),VR(:,:,vy),ix,jx,2,lm_type)
     call limiter(V(:,:,vz),VL(:,:,vz),VR(:,:,vz),ix,jx,2,lm_type)
     call limiter(V(:,:,pr),VL(:,:,pr),VR(:,:,pr),ix,jx,2,lm_type)
     call limiter(U(:,:,ro),VL(:,:,ro),VR(:,:,ro),ix,jx,2,lm_type)
     call limiter(U(:,:,bx),VL(:,:,bx),VR(:,:,bx),ix,jx,2,lm_type)
     call limiter(U(:,:,by),VL(:,:,by),VR(:,:,by),ix,jx,2,lm_type)
     call limiter(U(:,:,bz),VL(:,:,bz),VR(:,:,bz),ix,jx,2,lm_type)
     call limiter(U(:,:,ps),VL(:,:,ps),VR(:,:,ps),ix,jx,2,lm_type)
!    fix flux bc (G)
     call bc_for_G(VL,VR,ix,jx)
!    Numerical flux in the Y direction (G)
!     write(6,*) 'VL, VR --> G'
     call flux_solver(G,VL,VR,ch,ix,jx,2,flux_type)
     call flux_resistive(G,U,x,y,dx,ix,jx,2)

     if( time_type == 0 ) then
!       write(6,*) 'U* = U + (dt/dx) (F-F)'
        call rk_tvd21(U,U1,F,G,dt,dx,ix,jx)
     elseif( time_type == 1 ) then
!       write(6,*) 'U*(n+1/2) = U + (0.5 dt/dx) (F-F)'
        call rk_std21(U,U1,F,G,dt,dx,ix,jx)
     endif
!    boundary condition
     call bc_for_U(U1,ix,jx)

!     write(6,*) 'U* --> V'
     call u2v(U1,V,ix,jx)
!    Slope limiters on primitive variables
!     write(6,*) 'V --> VL, VR (F)'
     call limiter(V(:,:,vx),VL(:,:,vx),VR(:,:,vx),ix,jx,1,lm_type)
     call limiter(V(:,:,vy),VL(:,:,vy),VR(:,:,vy),ix,jx,1,lm_type)
     call limiter(V(:,:,vz),VL(:,:,vz),VR(:,:,vz),ix,jx,1,lm_type)
     call limiter(V(:,:,pr),VL(:,:,pr),VR(:,:,pr),ix,jx,1,lm_type)
     call limiter(U1(:,:,ro),VL(:,:,ro),VR(:,:,ro),ix,jx,1,lm_type)
     call limiter(U1(:,:,bx),VL(:,:,bx),VR(:,:,bx),ix,jx,1,lm_type)
     call limiter(U1(:,:,by),VL(:,:,by),VR(:,:,by),ix,jx,1,lm_type)
     call limiter(U1(:,:,bz),VL(:,:,bz),VR(:,:,bz),ix,jx,1,lm_type)
     call limiter(U1(:,:,ps),VL(:,:,ps),VR(:,:,ps),ix,jx,1,lm_type)
!    fix flux bc (F)
     call bc_for_F(VL,VR,ix,jx)
!    Numerical flux in the X direction (F)
!     write(6,*) 'VL, VR --> F'
     call flux_solver(F,VL,VR,ch,ix,jx,1,flux_type)
     call flux_resistive(F,U1,x,y,dx,ix,jx,1)

!    Slope limiters on primitive variables
!     write(6,*) 'V --> VL, VR (G)'
     call limiter(V(:,:,vx),VL(:,:,vx),VR(:,:,vx),ix,jx,2,lm_type)
     call limiter(V(:,:,vy),VL(:,:,vy),VR(:,:,vy),ix,jx,2,lm_type)
     call limiter(V(:,:,vz),VL(:,:,vz),VR(:,:,vz),ix,jx,2,lm_type)
     call limiter(V(:,:,pr),VL(:,:,pr),VR(:,:,pr),ix,jx,2,lm_type)
     call limiter(U1(:,:,ro),VL(:,:,ro),VR(:,:,ro),ix,jx,2,lm_type)
     call limiter(U1(:,:,bx),VL(:,:,bx),VR(:,:,bx),ix,jx,2,lm_type)
     call limiter(U1(:,:,by),VL(:,:,by),VR(:,:,by),ix,jx,2,lm_type)
     call limiter(U1(:,:,bz),VL(:,:,bz),VR(:,:,bz),ix,jx,2,lm_type)
     call limiter(U1(:,:,ps),VL(:,:,ps),VR(:,:,ps),ix,jx,2,lm_type)
!    fix flux bc (G)
     call bc_for_G(VL,VR,ix,jx)
!    Numerical flux in the Y direction (G)
!     write(6,*) 'VL, VR --> G'
     call flux_solver(G,VL,VR,ch,ix,jx,2,flux_type)
     call flux_resistive(G,U1,x,y,dx,ix,jx,2)

     if( time_type == 0 ) then
!       write(6,*) 'U_new = 0.5( U_old + U* + F dt )'
        call rk_tvd22(U,U1,F,G,dt,dx,ix,jx)
     elseif( time_type == 1 ) then
!       write(6,*) 'U_new = U + (dt/dx) (F-F) (n+1/2)'
        call rk_std22(U,F,G,dt,dx,ix,jx)
     endif

!    boundary condition
     call bc_for_U(U,ix,jx)

!    GLM solver for the second half timestep
     call glm_ss2(U,ch,dt,ix,jx)

     t=t+dt

  enddo

  write(6,*) '== end =='
end program main
!-----------------------------------------------------------------------
