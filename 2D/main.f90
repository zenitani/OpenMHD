program main
!-----------------------------------------------------------------------
!     OpenMHD  Riemann solver
!-----------------------------------------------------------------------
!     2010/09/18  S. Zenitani  2D LLF/HLL/HLLC-G/HLLD solver
!     2010/09/25  S. Zenitani  Hyperbolic divergence cleaning
!-----------------------------------------------------------------------
  implicit none
  include 'param.h'
  integer, parameter :: version = 20100928   ! version number
  integer, parameter :: ix = 202
  integer, parameter :: jx = 202
  integer, parameter :: loop_max = 30000
  real(8), parameter :: tend  = 4.0d0
  real(8), parameter :: dtout = 0.1d0 ! output interval
  real(8), parameter :: cfl   = 0.4d0 ! time step
! Slope limiter  (0: flat, 1: minmod, 2: MC, 3: van Leer, 4: Koren)
  integer, parameter :: lm_type   = 1
! Numerical flux (0: LLF, 1: HLL, 2: HLLC, 3: HLLD)
  integer, parameter :: flux_type = 3
! Time marching  (0: TVD RK2, 1: RK2)
  integer, parameter :: time_type = 0
!-----------------------------------------------------------------------
  real(8) :: x(ix), y(jx), dx
  real(8) :: U(ix,jx,var1)  ! conserved variables (U)
  real(8) :: U1(ix,jx,var1) ! conserved variables: medium state (U*)
  real(8) :: V(ix,jx,var2)  ! primitive variables (V)
  real(8) :: VL(ix,jx,var1), VR(ix,jx,var1) ! interpolated states
  real(8) :: F(ix,jx,var1), G(ix,jx,var1)   ! numerical flux (F,G)
!-----------------------------------------------------------------------
  integer :: k
  integer :: n_output
  real(8) :: t, dt, t_output
  real(8) :: ch
  character*256 :: filename
!-----------------------------------------------------------------------

  t    =  0.d0
  dt   =  0.d0
  call model(U,V,x,y,dx,ix,jx)
  call bc(U,ix,jx)
  call set_dt(U,V,ch,dt,dx,cfl,ix,jx)
  t_output = -dt/3.d0
  n_output =  0

  if ( dt .gt. dtout ) then
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
  do k=1,loop_max

     write(6,*) ' t = ', t
!    Recovering primitive variables
!     write(6,*) 'U --> V'
     call u2v(U,V,ix,jx)
!   -----------------  
!    [ output ]
     if ( t .ge. t_output ) then
        write(6,*) 'data output   t = ', t
        write(filename,990) n_output
990     format ('data/field-',i5.5,'.dat')
        call output(filename,ix,jx,t,x,y,U,V)
        n_output = n_output + 1
        t_output = t_output + dtout
     endif
!    [ end? ]
     if ( t .ge. tend )  exit
     if ( k .eq. loop_max ) then
        write(6,*) 'max loop'
        exit
     endif
!   -----------------  
!    CFL condition
     call set_dt(U,V,ch,dt,dx,cfl,ix,jx)
!    GLM solver for the first half timestep
!    This should be done after set_dt()
     call glm_ss(U,ch,0.5d0*dt,ix,jx)

!    Slope limiters on primitive variables
!     write(6,*) 'V --> VL, VR (F)'
     call limiter_f(V(1,1,vx),VL(1,1,vx),VR(1,1,vx),ix,jx,lm_type)
     call limiter_f(V(1,1,vy),VL(1,1,vy),VR(1,1,vy),ix,jx,lm_type)
     call limiter_f(V(1,1,vz),VL(1,1,vz),VR(1,1,vz),ix,jx,lm_type)
     call limiter_f(V(1,1,pr),VL(1,1,pr),VR(1,1,pr),ix,jx,lm_type)
     call limiter_f(U(1,1,ro),VL(1,1,ro),VR(1,1,ro),ix,jx,lm_type)
     call limiter_f(U(1,1,bx),VL(1,1,bx),VR(1,1,bx),ix,jx,lm_type)
     call limiter_f(U(1,1,by),VL(1,1,by),VR(1,1,by),ix,jx,lm_type)
     call limiter_f(U(1,1,bz),VL(1,1,bz),VR(1,1,bz),ix,jx,lm_type)
     call limiter_f(U(1,1,ps),VL(1,1,ps),VR(1,1,ps),ix,jx,lm_type)
!    fix VL/VR for periodic bc (F)
     VR(ix-1,:,:) = VR(1,:,:)
     VL(1,:,:)    = VL(ix-1,:,:)
!    Numerical flux in the X direction (F)
!     write(6,*) 'VL, VR --> F'
     if( flux_type .eq. 0 )then
        call llf_f(F,VL,VR,ix,jx)
     elseif( flux_type .eq. 1 )then
        call hll_f(F,VL,VR,ix,jx)
     elseif( flux_type .eq. 2 )then
        call hllc_f(F,VL,VR,ix,jx)
     elseif( flux_type .eq. 3 )then
        call hlld_f(F,VL,VR,ix,jx)
     endif
     call glm_f(F,VL,VR,ch,ix,jx)

!    Slope limiters on primitive variables
!     write(6,*) 'V --> VL, VR (G)'
     call limiter_g(V(1,1,vx),VL(1,1,vx),VR(1,1,vx),ix,jx,lm_type)
     call limiter_g(V(1,1,vy),VL(1,1,vy),VR(1,1,vy),ix,jx,lm_type)
     call limiter_g(V(1,1,vz),VL(1,1,vz),VR(1,1,vz),ix,jx,lm_type)
     call limiter_g(V(1,1,pr),VL(1,1,pr),VR(1,1,pr),ix,jx,lm_type)
     call limiter_g(U(1,1,ro),VL(1,1,ro),VR(1,1,ro),ix,jx,lm_type)
     call limiter_g(U(1,1,bx),VL(1,1,bx),VR(1,1,bx),ix,jx,lm_type)
     call limiter_g(U(1,1,by),VL(1,1,by),VR(1,1,by),ix,jx,lm_type)
     call limiter_g(U(1,1,bz),VL(1,1,bz),VR(1,1,bz),ix,jx,lm_type)
     call limiter_g(U(1,1,ps),VL(1,1,ps),VR(1,1,ps),ix,jx,lm_type)
!    fix VL/VR for periodic bc (G)
     VL(:,   1,:) = VL(:,jx-1,:)
     VR(:,jx-1,:) = VR(:,1,:)
!    Numerical flux in the Y direction (G)
!     write(6,*) 'VL, VR --> G'
     if( flux_type .eq. 0 )then
        call llf_g(G,VL,VR,ix,jx)
     elseif( flux_type .eq. 1 )then
        call hll_g(G,VL,VR,ix,jx)
     elseif( flux_type .eq. 2 )then
        call hllc_g(G,VL,VR,ix,jx)
     elseif( flux_type .eq. 3 )then
        call hlld_g(G,VL,VR,ix,jx)
     endif
     call glm_g(G,VL,VR,ch,ix,jx)

     if( time_type .eq. 0 ) then
!       write(6,*) 'U* = U + (dt/dx) (F-F)'
        call rk21(U,U1,F,G,dt,dx,ix,jx)
     elseif( time_type .eq. 1 ) then
!       write(6,*) 'U*(n+1/2) = U + (0.5 dt/dx) (F-F)'
        call step1(U,U1,F,G,dt,dx,ix,jx)
     endif
!    boundary condition
     call bc(U1,ix,jx)
!     write(6,*) 'U* --> V'
     call u2v(U1,V,ix,jx)

!    Slope limiters on primitive variables
!     write(6,*) 'V --> VL, VR (F)'
     call limiter_f(V(1,1,vx),VL(1,1,vx),VR(1,1,vx),ix,jx,lm_type)
     call limiter_f(V(1,1,vy),VL(1,1,vy),VR(1,1,vy),ix,jx,lm_type)
     call limiter_f(V(1,1,vz),VL(1,1,vz),VR(1,1,vz),ix,jx,lm_type)
     call limiter_f(V(1,1,pr),VL(1,1,pr),VR(1,1,pr),ix,jx,lm_type)
     call limiter_f(U1(1,1,ro),VL(1,1,ro),VR(1,1,ro),ix,jx,lm_type)
     call limiter_f(U1(1,1,bx),VL(1,1,bx),VR(1,1,bx),ix,jx,lm_type)
     call limiter_f(U1(1,1,by),VL(1,1,by),VR(1,1,by),ix,jx,lm_type)
     call limiter_f(U1(1,1,bz),VL(1,1,bz),VR(1,1,bz),ix,jx,lm_type)
     call limiter_f(U1(1,1,ps),VL(1,1,ps),VR(1,1,ps),ix,jx,lm_type)
!    fix VL/VR for periodic bc (F)
     VR(ix-1,:,:) = VR(1,:,:)
     VL(1,:,:)    = VL(ix-1,:,:)
!    Numerical flux in the X direction (F)
!     write(6,*) 'VL, VR --> F'
     if( flux_type .eq. 0 )then
        call llf_f(F,VL,VR,ix,jx)
     elseif( flux_type .eq. 1 )then
        call hll_f(F,VL,VR,ix,jx)
     elseif( flux_type .eq. 2 )then
        call hllc_f(F,VL,VR,ix,jx)
     elseif( flux_type .eq. 3 )then
        call hlld_f(F,VL,VR,ix,jx)
     endif        
     call glm_f(F,VL,VR,ch,ix,jx)

!    Slope limiters on primitive variables
!     write(6,*) 'V --> VL, VR (G)'
     call limiter_g(V(1,1,vx),VL(1,1,vx),VR(1,1,vx),ix,jx,lm_type)
     call limiter_g(V(1,1,vy),VL(1,1,vy),VR(1,1,vy),ix,jx,lm_type)
     call limiter_g(V(1,1,vz),VL(1,1,vz),VR(1,1,vz),ix,jx,lm_type)
     call limiter_g(V(1,1,pr),VL(1,1,pr),VR(1,1,pr),ix,jx,lm_type)
     call limiter_g(U1(1,1,ro),VL(1,1,ro),VR(1,1,ro),ix,jx,lm_type)
     call limiter_g(U1(1,1,bx),VL(1,1,bx),VR(1,1,bx),ix,jx,lm_type)
     call limiter_g(U1(1,1,by),VL(1,1,by),VR(1,1,by),ix,jx,lm_type)
     call limiter_g(U1(1,1,bz),VL(1,1,bz),VR(1,1,bz),ix,jx,lm_type)
     call limiter_g(U1(1,1,ps),VL(1,1,ps),VR(1,1,ps),ix,jx,lm_type)
!    fix VL/VR for periodic bc (G)
     VL(:,   1,:) = VL(:,jx-1,:)
     VR(:,jx-1,:) = VR(:,1,:)
!    Numerical flux in the Y direction (G)
!     write(6,*) 'VL, VR --> G'
     if( flux_type .eq. 0 )then
        call llf_g(G,VL,VR,ix,jx)
     elseif( flux_type .eq. 1 )then
        call hll_g(G,VL,VR,ix,jx)
     elseif( flux_type .eq. 2 )then
        call hllc_g(G,VL,VR,ix,jx)
     elseif( flux_type .eq. 3 )then
        call hlld_g(G,VL,VR,ix,jx)
     endif
     call glm_g(G,VL,VR,ch,ix,jx)

     if( time_type .eq. 0 ) then
!       write(6,*) 'U_new = 0.5( U_old + U* + F dt )'
        call rk22(U,U1,F,G,dt,dx,ix,jx)
     elseif( time_type .eq. 1 ) then
!       write(6,*) 'U_new = U + (dt/dx) (F-F) (n+1/2)'
        call step2(U,F,G,dt,dx,ix,jx)
     endif

!    GLM solver for the second half timestep
     call glm_ss(U,ch,0.5d0*dt,ix,jx)

!    boundary condition
     call bc(U,ix,jx)
     t=t+dt
  enddo

  write(6,*) '== end =='
end program main
!-----------------------------------------------------------------------
