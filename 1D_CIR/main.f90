program main
!-----------------------------------------------------------------------
!     A corotating interaction region (CIR) problem
!        Ref: K. Tsubouchi, J. Geophys. Res., 114, A02101 (2009)
!-----------------------------------------------------------------------
!     2011/09/27  S. Zenitani  LLF/HLLC-G/HLLD solver
!-----------------------------------------------------------------------
  implicit none
  include 'param.h'
  integer, parameter :: ix = 10000 + 2
  integer, parameter :: jx = 1
  integer, parameter :: loop_max = 1000000
  real(8), parameter :: tend  = 600.d0
  real(8), parameter :: dtout = 100.d0 ! output interval
!  real(8), parameter :: cfl   = 0.4d0 ! time step
  real(8), parameter :: cfl   = 0.8d0 ! time step
! Slope limiter  (0: flat, 1: minmod, 2: MC, 3: van Leer, 4: Koren)
  integer, parameter :: lm_type   = 1
! Numerical flux (0: LLF, 1: HLL, 2: HLLC, 3: HLLD)
  integer, parameter :: flux_type = 3
! Time marching  (0: TVD RK2, 1: RK2)
  integer, parameter :: time_type = 0
!-----------------------------------------------------------------------
  real(8) :: x(ix), y(jx), dx
  real(8) :: U(ix,jx,var1)  ! conserved variables (U)
  real(8) :: U1(ix,jx,var1) ! another conserved variables (U*) for RK
  real(8) :: V(ix,jx,var2)  ! primitive variables (V)
  real(8) :: VL(ix,jx,var1), VR(ix,jx,var1) ! interpolated states
  real(8) :: F(ix,jx,var1)  ! numerical flux (F)
!-----------------------------------------------------------------------
  integer :: i
  integer :: n_loop,n_output
  real(8) :: t, dt, t_output, dtx
  real(8), parameter :: zero = 0.d0 ! 0 to disable hyperbolic div cleaning
  real(8) :: vmax
  character*256 :: filename
!-----------------------------------------------------------------------

  dt   =  0.d0
  t    =  0.d0
  call model(U,V,x,y,dx,ix,jx)
! --- BC ---
  U(ix,:,:) = U(ix-1,:,:)
  U(1,:,:)  = U(2,:,:)
! --- BC ---
  call set_dt(U,V,vmax,dt,dx,cfl,ix,jx)
  t_output = -dt/3.d0
  n_output =  0

  if ( dt > dtout ) then
     write(6,*) 'error: ', dt, '>', dtout
  endif
  write(6,*) '[Params]'
  write(6,999) lm_type, flux_type, time_type
  write(6,998) dt, dtout, ix, jx
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
990     format ('data/x-',i5.5,'.dat')
        call output(filename,ix,jx,t,x,y,U,V)
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
     call set_dt(U,V,vmax,dt,dx,cfl,ix,jx)
     dtx = dt/dx
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
!    Numerical flux in the X direction (F)
!     write(6,*) 'VL, VR --> F'
     call flux_solver(F,VL,VR,zero,ix,jx,1,flux_type)

     if( time_type == 0 ) then
!       write(6,*) 'U* = U + (dt/dx) (F-F)'
        do i=2,ix-1
           U1(i,1,:) = U(i,1,:) + dtx*( F(i-1,1,:) - F(i,1,:) )
        enddo
     elseif( time_type == 1 ) then
!       write(6,*) 'U*(n+1/2) = U + (0.5 dt/dx) (F-F)'
        do i=2,ix-1
           U1(i,1,:) = U(i,1,:) + 0.5d0*dtx*( F(i-1,1,:) - F(i,1,:) )
        enddo
     endif
! --- BC ---
     U1(ix,:,:) = U1(ix-1,:,:)
     U1(1,:,:)  = U1(2,:,:)
! --- BC ---
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
!    Numerical flux in the X direction (F)
!     write(6,*) 'VL, VR --> F'
     call flux_solver(F,VL,VR,zero,ix,jx,1,flux_type)

     if( time_type == 0 ) then
!       write(6,*) 'U_new = 0.5( U_old + U* + F dt )'
        do i=2,ix-1
           U(i,1,:) = 0.5d0*( U(i,1,:)+U1(i,1,:) + dtx*( F(i-1,1,:)-F(i,1,:) ) )
        enddo
     elseif( time_type == 1 ) then
!       write(6,*) 'U_new = U + (dt/dx) (F-F) (n+1/2)'
        do i=2,ix-1
           U(i,1,:) = U(i,1,:) + dtx * ( F(i-1,1,:) - F(i,1,:) )
        enddo
     endif

! --- BC ---
     U(ix,:,:) = U(ix-1,:,:)
     U(1,:,:)  = U(2,:,:)
! --- BC ---
     t=t+dt
  enddo

  write(6,*) '== end =='
end program main
!-----------------------------------------------------------------------
