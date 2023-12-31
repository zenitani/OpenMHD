program main
!-----------------------------------------------------------------------
!     One-dimensional relativistic magnetohydrodynamic code
!     without the normal magnetic field (B_x=0)
!        Ref: S. Zenitani, M. Hesse, A. Klimas, ApJ, 712, 951 (2010)
!-----------------------------------------------------------------------
!     2009/08/01  S. Zenitani  1D HLL/HLLD solver
!     2010/04/26  S. Zenitani  a bug fixed (sound speed in HLL flux)
!     2017/07/20  S. Zenitani  changed primitive variables: v --> u
!     2017/07/21  S. Zenitani  included in OpenMHD
!-----------------------------------------------------------------------
  implicit none
  include 'param_rela.h'
  integer, parameter :: ix = 6400 + 2
  integer, parameter :: jx = 1
  integer, parameter :: loop_max = 30000
  real(8), parameter :: cfl   = 0.8d0 / 2  ! Two Riemann fans (R+L)
  real(8), parameter :: tend  = 0.4d0
  real(8), parameter :: dtout = 0.1d0 ! output interval
! Slope limiter  (0: flat, 1: minmod, 2: MC, 3: van Leer, 4: Koren)
  integer, parameter :: lm_type   = 2
! Time marching  (0: TVD RK2, 1: RK2)
  integer, parameter :: time_type = 0
  logical, parameter :: hlld = .true. ! HLLD or HLL
  integer :: i
  integer :: n_loop,n_output
  real(8) :: t, dt, t_output, dtx
  character*256 :: filename
!-----------------------------------------------------------------------
  real(8) :: x(ix),y(jx),dx
  real(8) :: U(ix,jx,var1)  ! conserved variables (U)
  real(8) :: U1(ix,jx,var1) ! another conserved variables (U*) for RK
  real(8) :: V(ix,jx,var2)  ! primitive variables (V)
  real(8) :: VL(ix,jx,var1), VR(ix,jx,var1) ! interpolated states
  real(8) :: F(ix,jx,var1)  ! numerical flux (F,G)
!-----------------------------------------------------------------------

  call model(U,x,y,dx,ix,jx)
  dt   =  cfl * dx
  dtx  =  dt/dx
  t    =  0.d0
  t_output = -dt/3.d0
  n_output =  0

  if ( dt > dtout ) then
     write(6,*) 'error: ', dt, '>', dtout
     stop
  endif
  write(6,*) '[Params]'
  write(6,998) dt, dtout, ix
998 format ('dt:', 1p, e10.3, ' dtout:', 1p, e10.3, ' grids:', i5 )
  write(6,*) '== start =='

!-----------------------------------------------------------------------
  do n_loop=1,loop_max

     write(6,*) ' t = ', t
!    U ==> V (1)
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
     if ( t >= tend ) then
        goto 1000
     endif
     if ( n_loop >= loop_max ) then
        write(6,*) 'max loop'
        goto 1000
     endif
!   -----------------
!    V ==> VR, VL
     call limiter(V(:,:,ux),VL(:,:,ux),VR(:,:,ux),ix,jx,1,lm_type)
     call limiter(V(:,:,uy),VL(:,:,uy),VR(:,:,uy),ix,jx,1,lm_type)
     call limiter(V(:,:,uz),VL(:,:,uz),VR(:,:,uz),ix,jx,1,lm_type)
     call limiter(V(:,:,pr),VL(:,:,pr),VR(:,:,pr),ix,jx,1,lm_type)
     call limiter(V(:,:,ro),VL(:,:,ro),VR(:,:,ro),ix,jx,1,lm_type)
     call limiter(U(:,:,bx),VL(:,:,bx),VR(:,:,bx),ix,jx,1,lm_type)
     call limiter(U(:,:,by),VL(:,:,by),VR(:,:,by),ix,jx,1,lm_type)
     call limiter(U(:,:,bz),VL(:,:,bz),VR(:,:,bz),ix,jx,1,lm_type)
     call limiter(U(:,:,ps),VL(:,:,ps),VR(:,:,ps),ix,jx,1,lm_type)
!    VR, VL ==> Riemann Flux
     call flux_solver(F,VL,VR,ix,jx,1,hlld)

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

!    U* ==> V (1)
!     write(6,*) 'U* --> V'
     call u2v(U1,V,ix,jx)
!    V ==> VR, VL
     call limiter(V(:,:,ux),VL(:,:,ux),VR(:,:,ux),ix,jx,1,lm_type)
     call limiter(V(:,:,uy),VL(:,:,uy),VR(:,:,uy),ix,jx,1,lm_type)
     call limiter(V(:,:,uz),VL(:,:,uz),VR(:,:,uz),ix,jx,1,lm_type)
     call limiter(V(:,:,pr),VL(:,:,pr),VR(:,:,pr),ix,jx,1,lm_type)
     call limiter(V(:,:,ro),VL(:,:,ro),VR(:,:,ro),ix,jx,1,lm_type)
     call limiter(U1(:,:,bx),VL(:,:,bx),VR(:,:,bx),ix,jx,1,lm_type)
     call limiter(U1(:,:,by),VL(:,:,by),VR(:,:,by),ix,jx,1,lm_type)
     call limiter(U1(:,:,bz),VL(:,:,bz),VR(:,:,bz),ix,jx,1,lm_type)
     call limiter(U1(:,:,ps),VL(:,:,ps),VR(:,:,ps),ix,jx,1,lm_type)
!    VR, VL ==> Riemann Flux
     call flux_solver(F,VL,VR,ix,jx,1,hlld)
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

1000 continue

  write(6,*) '== end =='
end program main
!-----------------------------------------------------------------------
