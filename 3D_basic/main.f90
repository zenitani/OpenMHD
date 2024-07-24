program main
!-----------------------------------------------------------------------
!     OpenMHD  Riemann solver
!-----------------------------------------------------------------------
!     2010/09/18  S. Zenitani  2D LLF/HLL/HLLC-G/HLLD solver
!     2023/12/25  S. Zenitani  3-D version
!-----------------------------------------------------------------------
  implicit none
  include 'param.h'
  integer, parameter :: ix = 100 + 2
  integer, parameter :: jx = 100 + 2
  integer, parameter :: kx = 200 + 2
  integer, parameter :: loop_max = 3000000
  real(8), parameter :: tend  = 1.01d0
  real(8), parameter :: dtout = 0.1d0 ! output interval
  real(8), parameter :: cfl   = 0.3d0 ! time step
! Slope limiter  (0: flat, 1: minmod, 2: MC, 3: van Leer, 4: Koren)
  integer, parameter :: lm_type   = 1
! Numerical flux (0: LLF, 1: HLL, 2: HLLC, 3: HLLD)
  integer, parameter :: flux_type = 3
! Time marching  (0: TVD RK2, 1: RK2)
  integer, parameter :: time_type = 0
!-----------------------------------------------------------------------
  real(8) :: x(ix), y(jx), z(kx), dx
  real(8) :: U(ix,jx,kx,var1)  ! conserved variables (U)
  real(8) :: U1(ix,jx,kx,var1) ! conserved variables: medium state (U*)
  real(8) :: V(ix,jx,kx,var2)  ! primitive variables (V)
  real(8) :: VL(ix,jx,kx,var1), VR(ix,jx,kx,var1) ! interpolated states
  real(8) :: F(ix,jx,kx,var1), G(ix,jx,kx,var1), H(ix,jx,kx,var1) ! numerical flux (F,G,H)
!-----------------------------------------------------------------------
  integer :: n_loop,n_output
  real(8) :: t, dt, t_output
  real(8) :: ch
  character*256 :: filename
!-----------------------------------------------------------------------

  t    =  0.d0
  dt   =  0.d0
  call model(U,V,x,y,z,dx,ix,jx,kx)
  call bc(U,ix,jx,kx)
  call set_dt(U,V,ch,dt,dx,cfl,ix,jx,kx)
  t_output = -dt/3.d0
  n_output =  0

  if ( dt > dtout ) then
     write(6,*) 'error: ', dt, '>', dtout
     stop
  endif
  write(6,*) '[Params]'
  write(6,998) dt, dtout, ix, jx, kx
  write(6,999) lm_type, flux_type, time_type
998 format (' dt:', 1p, e10.3, ' dtout:', 1p, e10.3, ' grids:', i5, i5, i5 )
999 format (' limiter: ', i1, '  flux: ', i1, '  time-marching: ', i1 )
  write(6,*) '== start =='

!-----------------------------------------------------------------------
  do n_loop=1,loop_max

     write(6,*) ' t = ', t
!    Recovering primitive variables
!     write(6,*) 'U --> V'
     call u2v(U,V,ix,jx,kx)
!   -----------------
!    [ output ]
     if ( t >= t_output ) then
        write(6,*) 'data output   t = ', t
        write(filename,990) n_output
990     format ('data/field-',i5.5,'.dat')
        call fileio_output(filename,ix,jx,kx,t,x,y,z,U,V)
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
     call set_dt(U,V,ch,dt,dx,cfl,ix,jx,kx)
!    GLM solver for the first half timestep
!    This should be done after set_dt()
     call glm_ss2(U,ch,dt,ix,jx,kx)

!    Slope limiters on primitive variables
!     write(6,*) 'V --> VL, VR (F)'
     call limiter(V(:,:,:,vx),VL(:,:,:,vx),VR(:,:,:,vx),ix,jx,kx,1,lm_type)
     call limiter(V(:,:,:,vy),VL(:,:,:,vy),VR(:,:,:,vy),ix,jx,kx,1,lm_type)
     call limiter(V(:,:,:,vz),VL(:,:,:,vz),VR(:,:,:,vz),ix,jx,kx,1,lm_type)
     call limiter(V(:,:,:,pr),VL(:,:,:,pr),VR(:,:,:,pr),ix,jx,kx,1,lm_type)
     call limiter(U(:,:,:,ro),VL(:,:,:,ro),VR(:,:,:,ro),ix,jx,kx,1,lm_type)
     call limiter(U(:,:,:,bx),VL(:,:,:,bx),VR(:,:,:,bx),ix,jx,kx,1,lm_type)
     call limiter(U(:,:,:,by),VL(:,:,:,by),VR(:,:,:,by),ix,jx,kx,1,lm_type)
     call limiter(U(:,:,:,bz),VL(:,:,:,bz),VR(:,:,:,bz),ix,jx,kx,1,lm_type)
     call limiter(U(:,:,:,ps),VL(:,:,:,ps),VR(:,:,:,ps),ix,jx,kx,1,lm_type)
!    fix VL/VR for periodic bc (F)
     VR(ix-1,:,:,:) = VR(   1,:,:,:)
     VL(   1,:,:,:) = VL(ix-1,:,:,:)
!    *** education ***
!    To disable hyperbolic div cleaning, pass zero to flux_solver
!    ch = 0.d0
!    *** education ***
!    Numerical flux in the X direction (F)
!     write(6,*) 'VL, VR --> F'
     call flux_solver(F,VL,VR,ch,ix,jx,kx,1,flux_type)

!    Slope limiters on primitive variables
!     write(6,*) 'V --> VL, VR (G)'
     call limiter(V(:,:,:,vx),VL(:,:,:,vx),VR(:,:,:,vx),ix,jx,kx,2,lm_type)
     call limiter(V(:,:,:,vy),VL(:,:,:,vy),VR(:,:,:,vy),ix,jx,kx,2,lm_type)
     call limiter(V(:,:,:,vz),VL(:,:,:,vz),VR(:,:,:,vz),ix,jx,kx,2,lm_type)
     call limiter(V(:,:,:,pr),VL(:,:,:,pr),VR(:,:,:,pr),ix,jx,kx,2,lm_type)
     call limiter(U(:,:,:,ro),VL(:,:,:,ro),VR(:,:,:,ro),ix,jx,kx,2,lm_type)
     call limiter(U(:,:,:,bx),VL(:,:,:,bx),VR(:,:,:,bx),ix,jx,kx,2,lm_type)
     call limiter(U(:,:,:,by),VL(:,:,:,by),VR(:,:,:,by),ix,jx,kx,2,lm_type)
     call limiter(U(:,:,:,bz),VL(:,:,:,bz),VR(:,:,:,bz),ix,jx,kx,2,lm_type)
     call limiter(U(:,:,:,ps),VL(:,:,:,ps),VR(:,:,:,ps),ix,jx,kx,2,lm_type)
!    fix VL/VR for periodic bc (G)
     VL(:,   1,:,:) = VL(:,jx-1,:,:)
     VR(:,jx-1,:,:) = VR(:,   1,:,:)
!    Numerical flux in the Y direction (G)
!     write(6,*) 'VL, VR --> G'
     call flux_solver(G,VL,VR,ch,ix,jx,kx,2,flux_type)

!    Slope limiters on primitive variables
!     write(6,*) 'V --> VL, VR (H)'
     call limiter(V(:,:,:,vx),VL(:,:,:,vx),VR(:,:,:,vx),ix,jx,kx,3,lm_type)
     call limiter(V(:,:,:,vy),VL(:,:,:,vy),VR(:,:,:,vy),ix,jx,kx,3,lm_type)
     call limiter(V(:,:,:,vz),VL(:,:,:,vz),VR(:,:,:,vz),ix,jx,kx,3,lm_type)
     call limiter(V(:,:,:,pr),VL(:,:,:,pr),VR(:,:,:,pr),ix,jx,kx,3,lm_type)
     call limiter(U(:,:,:,ro),VL(:,:,:,ro),VR(:,:,:,ro),ix,jx,kx,3,lm_type)
     call limiter(U(:,:,:,bx),VL(:,:,:,bx),VR(:,:,:,bx),ix,jx,kx,3,lm_type)
     call limiter(U(:,:,:,by),VL(:,:,:,by),VR(:,:,:,by),ix,jx,kx,3,lm_type)
     call limiter(U(:,:,:,bz),VL(:,:,:,bz),VR(:,:,:,bz),ix,jx,kx,3,lm_type)
     call limiter(U(:,:,:,ps),VL(:,:,:,ps),VR(:,:,:,ps),ix,jx,kx,3,lm_type)
!    fix VL/VR for periodic bc (H)
     VL(:,:,   1,:) = VL(:,:,kx-1,:)
     VR(:,:,kx-1,:) = VR(:,:,   1,:)
!    Numerical flux in the Z direction (H)
!     write(6,*) 'VL, VR --> H'
     call flux_solver(H,VL,VR,ch,ix,jx,kx,3,flux_type)

     if( time_type == 0 ) then
!       write(6,*) 'U* = U + (dt/dx) (F-F)'
        call rk_tvd21(U,U1,F,G,H,dt,dx,ix,jx,kx)
     elseif( time_type == 1 ) then
!       write(6,*) 'U*(n+1/2) = U + (0.5 dt/dx) (F-F)'
        call rk_std21(U,U1,F,G,H,dt,dx,ix,jx,kx)
     endif
!    boundary condition
     call bc(U1,ix,jx,kx)
!     write(6,*) 'U* --> V'
     call u2v(U1,V,ix,jx,kx)

!    Slope limiters on primitive variables
!     write(6,*) 'V --> VL, VR (F)'
     call limiter(V(:,:,:,vx),VL(:,:,:,vx),VR(:,:,:,vx),ix,jx,kx,1,lm_type)
     call limiter(V(:,:,:,vy),VL(:,:,:,vy),VR(:,:,:,vy),ix,jx,kx,1,lm_type)
     call limiter(V(:,:,:,vz),VL(:,:,:,vz),VR(:,:,:,vz),ix,jx,kx,1,lm_type)
     call limiter(V(:,:,:,pr),VL(:,:,:,pr),VR(:,:,:,pr),ix,jx,kx,1,lm_type)
     call limiter(U1(:,:,:,ro),VL(:,:,:,ro),VR(:,:,:,ro),ix,jx,kx,1,lm_type)
     call limiter(U1(:,:,:,bx),VL(:,:,:,bx),VR(:,:,:,bx),ix,jx,kx,1,lm_type)
     call limiter(U1(:,:,:,by),VL(:,:,:,by),VR(:,:,:,by),ix,jx,kx,1,lm_type)
     call limiter(U1(:,:,:,bz),VL(:,:,:,bz),VR(:,:,:,bz),ix,jx,kx,1,lm_type)
     call limiter(U1(:,:,:,ps),VL(:,:,:,ps),VR(:,:,:,ps),ix,jx,kx,1,lm_type)
!    fix VL/VR for periodic bc (F)
     VR(ix-1,:,:,:) = VR(   1,:,:,:)
     VL(   1,:,:,:) = VL(ix-1,:,:,:)
!    Numerical flux in the X direction (F)
!     write(6,*) 'VL, VR --> F'
     call flux_solver(F,VL,VR,ch,ix,jx,kx,1,flux_type)

!    Slope limiters on primitive variables
!     write(6,*) 'V --> VL, VR (G)'
     call limiter(V(:,:,:,vx),VL(:,:,:,vx),VR(:,:,:,vx),ix,jx,kx,2,lm_type)
     call limiter(V(:,:,:,vy),VL(:,:,:,vy),VR(:,:,:,vy),ix,jx,kx,2,lm_type)
     call limiter(V(:,:,:,vz),VL(:,:,:,vz),VR(:,:,:,vz),ix,jx,kx,2,lm_type)
     call limiter(V(:,:,:,pr),VL(:,:,:,pr),VR(:,:,:,pr),ix,jx,kx,2,lm_type)
     call limiter(U1(:,:,:,ro),VL(:,:,:,ro),VR(:,:,:,ro),ix,jx,kx,2,lm_type)
     call limiter(U1(:,:,:,bx),VL(:,:,:,bx),VR(:,:,:,bx),ix,jx,kx,2,lm_type)
     call limiter(U1(:,:,:,by),VL(:,:,:,by),VR(:,:,:,by),ix,jx,kx,2,lm_type)
     call limiter(U1(:,:,:,bz),VL(:,:,:,bz),VR(:,:,:,bz),ix,jx,kx,2,lm_type)
     call limiter(U1(:,:,:,ps),VL(:,:,:,ps),VR(:,:,:,ps),ix,jx,kx,2,lm_type)
!    fix VL/VR for periodic bc (G)
     VL(:,   1,:,:) = VL(:,jx-1,:,:)
     VR(:,jx-1,:,:) = VR(:,   1,:,:)
!    Numerical flux in the Y direction (G)
!     write(6,*) 'VL, VR --> G'
     call flux_solver(G,VL,VR,ch,ix,jx,kx,2,flux_type)

!    Slope limiters on primitive variables
!     write(6,*) 'V --> VL, VR (H)'
     call limiter(V(:,:,:,vx),VL(:,:,:,vx),VR(:,:,:,vx),ix,jx,kx,3,lm_type)
     call limiter(V(:,:,:,vy),VL(:,:,:,vy),VR(:,:,:,vy),ix,jx,kx,3,lm_type)
     call limiter(V(:,:,:,vz),VL(:,:,:,vz),VR(:,:,:,vz),ix,jx,kx,3,lm_type)
     call limiter(V(:,:,:,pr),VL(:,:,:,pr),VR(:,:,:,pr),ix,jx,kx,3,lm_type)
     call limiter(U1(:,:,:,ro),VL(:,:,:,ro),VR(:,:,:,ro),ix,jx,kx,3,lm_type)
     call limiter(U1(:,:,:,bx),VL(:,:,:,bx),VR(:,:,:,bx),ix,jx,kx,3,lm_type)
     call limiter(U1(:,:,:,by),VL(:,:,:,by),VR(:,:,:,by),ix,jx,kx,3,lm_type)
     call limiter(U1(:,:,:,bz),VL(:,:,:,bz),VR(:,:,:,bz),ix,jx,kx,3,lm_type)
     call limiter(U1(:,:,:,ps),VL(:,:,:,ps),VR(:,:,:,ps),ix,jx,kx,3,lm_type)
!    fix VL/VR for periodic bc (H)
     VL(:,:,   1,:) = VL(:,:,kx-1,:)
     VR(:,:,kx-1,:) = VR(:,:,   1,:)
!    Numerical flux in the Z direction (H)
!     write(6,*) 'VL, VR --> H'
     call flux_solver(H,VL,VR,ch,ix,jx,kx,3,flux_type)

     if( time_type == 0 ) then
!       write(6,*) 'U_new = 0.5( U_old + U* + F dt )'
        call rk_tvd22(U,U1,F,G,H,dt,dx,ix,jx,kx)
     elseif( time_type == 1 ) then
!       write(6,*) 'U_new = U + (dt/dx) (F-F) (n+1/2)'
        call rk_std22(U,F,G,H,dt,dx,ix,jx,kx)
     endif

!    boundary condition
     call bc(U,ix,jx,kx)

!    GLM solver for the second half timestep
     call glm_ss2(U,ch,dt,ix,jx,kx)

     t=t+dt
  enddo

  write(6,*) '== end =='
end program main
!-----------------------------------------------------------------------
