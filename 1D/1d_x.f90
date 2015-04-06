program main
!-----------------------------------------------------------------------
!     Open MHD  Riemann solver
!-----------------------------------------------------------------------
!     2010/05/12  S. Zenitani  LLF/HLLC-G/HLLD solver
!-----------------------------------------------------------------------
  implicit none
  include 'param.h'
  integer, parameter :: version = 20100917   ! version number
  integer, parameter :: ix = 302
  integer, parameter :: jx = 1
  integer, parameter :: loop_max = 30000
  real(8), parameter :: tend  = 0.2d0
  real(8), parameter :: dtout = 0.1d0 ! output interval
  real(8), parameter :: cfl   = 0.8d0 ! time step
! Slope limiter  (0: flat, 1: minmod, 2: MC, 3: van Leer, 4: Koren)
  integer, parameter :: lm_type   = 3
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
  real(8) :: F(ix,jx,var1)  ! numerical flux (F,G)
!-----------------------------------------------------------------------
  integer :: i, k
  integer :: n_output
  real(8) :: t, dt, t_output, dtx
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

  if ( dt .gt. dtout ) then
     write(6,*) 'error: ', dt, '>', dtout
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
990     format ('data/x-',i5.5,'.dat')
        call output(filename,ix,jx,t,x,y,U,V)
        n_output = n_output + 1
        t_output = t_output + dtout
     endif
!    [ end? ]
     if ( k .eq. loop_max ) then
        write(6,*) 'max loop'
        exit
     endif
     if ( t .ge. tend )  exit
!   -----------------  
!    CFL condition
     call set_dt(U,V,vmax,dt,dx,cfl,ix,jx)
     dtx = dt/dx
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

     if( time_type .eq. 0 ) then
!       write(6,*) 'U* = U + (dt/dx) (F-F)'
        do i=2,ix-1
           U1(i,1,:) = U(i,1,:) + dtx*( F(i-1,1,:) - F(i,1,:) )
        enddo
     elseif( time_type .eq. 1 ) then
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
     call limiter_f(V(1,1,vx),VL(1,1,vx),VR(1,1,vx),ix,jx,lm_type)
     call limiter_f(V(1,1,vy),VL(1,1,vy),VR(1,1,vy),ix,jx,lm_type)
     call limiter_f(V(1,1,vz),VL(1,1,vz),VR(1,1,vz),ix,jx,lm_type)
     call limiter_f(V(1,1,pr),VL(1,1,pr),VR(1,1,pr),ix,jx,lm_type)
     call limiter_f(U1(1,1,ro),VL(1,1,ro),VR(1,1,ro),ix,jx,lm_type)
     call limiter_f(U1(1,1,bx),VL(1,1,bx),VR(1,1,bx),ix,jx,lm_type)
     call limiter_f(U1(1,1,by),VL(1,1,by),VR(1,1,by),ix,jx,lm_type)
     call limiter_f(U1(1,1,bz),VL(1,1,bz),VR(1,1,bz),ix,jx,lm_type)
     call limiter_f(U1(1,1,ps),VL(1,1,ps),VR(1,1,ps),ix,jx,lm_type)
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

     if( time_type .eq. 0 ) then
!       write(6,*) 'U_new = 0.5( U_old + U* + F dt )'
        do i=2,ix-1
           U(i,1,:) = 0.5d0*( U(i,1,:)+U1(i,1,:) + dtx*( F(i-1,1,:)-F(i,1,:) ) )
        enddo
     elseif( time_type .eq. 1 ) then
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
