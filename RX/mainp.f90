program main
!-----------------------------------------------------------------------
!     Open MHD  Reconnection solver (parallel version)
!-----------------------------------------------------------------------
!     2010/09/25  S. Zenitani  HLL reconnection code
!     2010/09/29  S. Zenitani  MPI version
!     2015/04/05  S. Zenitani  MPI-IO
!-----------------------------------------------------------------------
  implicit none
  include 'mpif.h' ! for MPI
  include 'param.h'
  integer, parameter :: version = 20150405   ! version number
  integer, parameter :: ix = 152
  integer, parameter :: jx = 602
  integer, parameter :: loop_max = 30000
  real(8), parameter :: tend  = 120.0d0
  real(8), parameter :: dtout = 5.d0 ! output interval
  real(8), parameter :: cfl   = 0.4d0 ! time step
! Slope limiter  (0: flat, 1: minmod, 2: MC, 3: van Leer, 4: Koren)
  integer, parameter :: lm_type   = 1
! Numerical flux (1: HLL, 3: HLLD)
  integer, parameter :: flux_type = 3
! Time marching  (0: TVD RK2, 1: RK2)
  integer, parameter :: time_type = 0
! Resistivity
  real(8), parameter :: Rm1 = 60.d0, Rm0 = 500.d0
!-----------------------------------------------------------------------
! See also modelp.f90
!-----------------------------------------------------------------------
  integer :: k
  integer :: n_output
  real(8) :: t, dt, t_output, dtg
  real(8) :: ch, chg
  character*256 :: filename
  integer :: merr, myrank, npe          ! for MPI
!-----------------------------------------------------------------------
  real(8) :: x(ix), y(jx), dx
  real(8) :: U(ix,jx,var1)  ! conserved variables (U)
  real(8) :: U1(ix,jx,var1) ! conserved variables: medium state (U*)
  real(8) :: V(ix,jx,var2)  ! primitive variables (V)
  real(8) :: VL(ix,jx,var1), VR(ix,jx,var1) ! interpolated states
  real(8) :: F(ix,jx,var1), G(ix,jx,var1)   ! numerical flux (F,G)
  real(8) :: E(ix,jx),EF(ix,jx), EG(ix,jx)  ! resistivity for U, F, G
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! for MPI
  call mpi_init(merr)
  call mpi_comm_size(mpi_comm_world,npe   ,merr)
  call mpi_comm_rank(mpi_comm_world,myrank,merr)
!-----------------------------------------------------------------------

  t    =  0.d0
  dt   =  0.d0
  call modelp(U,V,x,y,dx,ix,jx,myrank,npe)
  call set_eta(E,EF,EG,x,y,dx,Rm1,Rm0,ix,jx)
  call mpibc(U,ix,jx,myrank,npe)
  call set_dt(U,V,ch,dt,dx,cfl,ix,jx)
  call set_dt2(Rm1,dt,dx,cfl)

  call mpi_allreduce(dt,dtg,1,mpi_double_precision, &
       mpi_min,mpi_comm_world,merr)
  call mpi_allreduce(ch,chg,1,mpi_double_precision, &
       mpi_max,mpi_comm_world,merr)
  call mpi_barrier(mpi_comm_world,merr)
  dt = dtg
  ch = chg
  t_output = -dt/3.d0
  n_output =  0

  if ( dt .gt. dtout ) then
     write(6,*) 'error: ', dt, '>', dtout
     stop
  endif
  if( myrank.eq.0 ) then
     write(6,*) '[Params]'
     write(6,*) 'Code version: ', version, '  Core # : ', npe
     write(6,*) 'Reynolds #  : ', Rm1, ' and ', Rm0
     write(6,998) dt, dtout, npe*(ix-2)+2, ix, jx
     write(6,999) lm_type, flux_type, time_type
998  format (' dt: ',e10.3,' dtout: ',e10.3,' grids:',i5,' (',i5,') x ',i5 )
999  format (' limiter: ', i1, '  flux: ', i1, '  time-marching: ', i1 )
     write(6,*) '== start =='
  endif

!-----------------------------------------------------------------------
  do k=1,loop_max

     if( myrank.eq.0 ) then
        write(6,*) ' t = ', t
     endif
!    Recovering primitive variables
!     write(6,*) 'U --> V'
     call u2v(U,V,ix,jx)
!   -----------------  
!    [ output ]
     if ( t .ge. t_output ) then
!        write(6,*) 'data output   t = ', t
!        write(filename,990) myrank, n_output
!        call output(filename,ix,jx,t,x,y,U,V)
        if( myrank.eq.0 )  write(6,*) 'writing data ...   t = ', t
        write(filename,980) n_output
        call mpioutput(filename,ix,jx,t,x,y,U,V,myrank,npe)
        n_output = n_output + 1
        t_output = t_output + dtout
     endif
!    [ end? ]
     if ( t .ge. tend )  exit
     if ( k .eq. loop_max ) then
        if( myrank.eq.0 )  write(6,*) 'max loop'
        exit
     endif
!   -----------------  
!    CFL condition
     call set_dt(U,V,ch,dt,dx,cfl,ix,jx)
     call set_dt2(Rm1,dt,dx,cfl)
     call mpi_allreduce(dt,dtg,1,mpi_double_precision, &
          mpi_min,mpi_comm_world,merr)
     call mpi_allreduce(ch,chg,1,mpi_double_precision, &
          mpi_max,mpi_comm_world,merr)
     dt = dtg
     ch = chg

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
!     write(6,*) 'fix VL/VR at MPI boundary'
     call mpibc_vlvr_f(VL,VR,ix,jx,myrank,npe)
!    Numerical flux in the X direction (F)
!     write(6,*) 'VL, VR --> F'
     if( flux_type .eq. 1 )then
        call hll_resistive_f(F,U,VL,VR,EF,dx,ix,jx)
     elseif( flux_type .eq. 3 )then
        call hlld_resistive_f(F,U,VL,VR,EF,dx,ix,jx)
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
!    fix flux bc (G)
     call bc_vlvr_g(VL,VR,ix,jx)
!    Numerical flux in the Y direction (G)
!     write(6,*) 'VL, VR --> G'
     if( flux_type .eq. 1 )then
        call hll_resistive_g(G,U,VL,VR,EG,dx,ix,jx)
     elseif( flux_type .eq. 3 )then
        call hlld_resistive_g(G,U,VL,VR,EG,dx,ix,jx)
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
     call mpibc(U1,ix,jx,myrank,npe)
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
!     write(6,*) 'fix VL/VR at MPI boundary'
     call mpibc_vlvr_f(VL,VR,ix,jx,myrank,npe)
!    Numerical flux in the X direction (F)
!     write(6,*) 'VL, VR --> F'
     if( flux_type .eq. 1 )then
        call hll_resistive_f(F,U1,VL,VR,EF,dx,ix,jx)
     elseif( flux_type .eq. 3 )then
        call hlld_resistive_f(F,U1,VL,VR,EF,dx,ix,jx)
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
!    fix flux bc (G)
     call bc_vlvr_g(VL,VR,ix,jx)
!    Numerical flux in the Y direction (G)
!     write(6,*) 'VL, VR --> G'
     if( flux_type .eq. 1 )then
        call hll_resistive_g(G,U1,VL,VR,EG,dx,ix,jx)
     elseif( flux_type .eq. 3 )then
        call hlld_resistive_g(G,U1,VL,VR,EG,dx,ix,jx)
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
     call mpibc(U,ix,jx,myrank,npe)
     t=t+dt
  enddo
!-----------------------------------------------------------------------

  call mpi_finalize(merr)
  if( myrank.eq.0 )  write(6,*) '== end =='

980 format ('data/field-',i5.5,'.dat')
!990 format ('data/field-',i3.3,'-',i5.5,'.dat')

end program main
!-----------------------------------------------------------------------
