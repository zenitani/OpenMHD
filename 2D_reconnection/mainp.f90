program main
!-----------------------------------------------------------------------
!     OpenMHD  Reconnection solver (parallel version)
!-----------------------------------------------------------------------
!     2010/09/25  S. Zenitani  HLL reconnection code
!     2010/09/29  S. Zenitani  MPI version
!     2015/04/05  S. Zenitani  MPI-IO
!-----------------------------------------------------------------------
  use parallel
  implicit none
  include 'param.h'
  integer, parameter :: ix = 150 + 2  ! 150 (cells per core) x 8 (cores) = 1200
  integer, parameter :: jx = 600 + 2
  integer, parameter :: mpi_nums(2)       = (/8, 1/)  ! MPI numbers
  logical, parameter :: bc_periodicity(2) = (/.false., .false./)
  integer, parameter :: loop_max = 30000
  real(8), parameter :: tend  = 200.0d0
  real(8), parameter :: dtout =  10.0d0 ! output interval
  real(8), parameter :: cfl   =   0.4d0 ! time step
  integer, parameter :: n_start = 0     ! If non-zero, load previous data file
! Slope limiter  (0: flat, 1: minmod, 2: MC, 3: van Leer, 4: Koren)
  integer, parameter :: lm_type   = 1
! Numerical flux (0: LLF, 1: HLL, 2: HLLC, 3: HLLD)
  integer, parameter :: flux_type = 3
! Time marching  (0: TVD RK2, 1: RK2)
  integer, parameter :: time_type = 0
! File I/O  (0: Standard, 1: MPI-IO)
  integer, parameter :: io_type   = 1
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
  integer :: merr, myrank, mreq(2)   ! for MPI
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
!  call mpi_init(merr)
!  call mpi_comm_size(mpi_comm_world,npe   ,merr)
!  call mpi_comm_rank(mpi_comm_world,myrank,merr)
  call parallel_init(mpi_nums,bc_periodicity,ix,jx)
  myrank = ranks%myrank
!-----------------------------------------------------------------------

  t    =  0.d0
  dt   =  0.d0
  call modelp(U,V,x,y,dx,ix,jx)
  call set_eta(E,EF,EG,x,y,dx,Rm1,Rm0,ix,jx)
  call parallel_exchange(U,ix,jx,1)
  call parallel_exchange(U,ix,jx,2)
  call mpibc_for_U(U,ix,jx)
  call set_dt(U,V,ch,dt,dx,cfl,ix,jx)
  call set_dt2(Rm1,dt,dx,cfl)

!  call mpi_iallreduce(mpi_in_place,ch,1,mpi_real8,mpi_max,cart2d%comm,mreq(1),merr)
!  call mpi_iallreduce(mpi_in_place,dt,1,mpi_real8,mpi_min,cart2d%comm,mreq(2),merr)
  call mpi_iallreduce(ch,chg,1,mpi_real8,mpi_max,cart2d%comm,mreq(1),merr)
  call mpi_iallreduce(dt,dtg,1,mpi_real8,mpi_min,cart2d%comm,mreq(2),merr)
  call mpi_waitall(2,mreq,mpi_statuses_ignore,merr)
  dt = dtg; ch = chg
  call mpi_barrier(cart2d%comm,merr)

  if ( dt > dtout ) then
     write(6,*) 'error: ', dt, '>', dtout
     stop
  endif
  if( myrank == 0 ) then
     write(6,*) '[Params]'
     write(6,997) version, ranks%size, cart2d%sizes
     write(6,998) dt, dtout, cart2d%sizes(1)*(ix-2)+2, ix, cart2d%sizes(2)*(jx-2)+2, jx
     write(6,999) lm_type, flux_type, time_type
997  format ('Code version: ', i8, '  MPI node # : ', i5,' (',i4,' x ',i4,')' )
998  format (' dt: ',e10.3,' dtout: ',e10.3,' grids:',i6,' (',i5,') x ',i6,' (',i5,') ')
999  format (' limiter: ', i1, '  flux: ', i1, '  time-marching: ', i1 )
     write(6,*) 'Reynolds #  : ', Rm1, ' and ', Rm0
     write(6,*) '== start =='
  endif
  call mpi_barrier(cart2d%comm,merr)

  if ( n_start /= 0 ) then
     if ( io_type == 0 ) then
        write(6,*) 'reading data ...   rank = ', myrank
        write(filename,990) myrank, n_start
        call input(filename,ix,jx,t,x,y,U)
     else
        if( myrank == 0 )  write(6,*) 'reading data ...'
        write(filename,980) n_start
        call mpiinput(filename,ix,jx,t,x,y,U)
     endif
  endif
  t_output = t - dt/3.d0
  n_output = n_start

!-----------------------------------------------------------------------
  do k=1,loop_max

     if( myrank == 0 ) then
        write(6,*) ' t = ', t
     endif
!    Recovering primitive variables
!     write(6,*) 'U --> V'
     call u2v(U,V,ix,jx)
!   -----------------  
!    [ output ]
     if ( t >= t_output ) then
        if (( k > 1 ).or.( n_start == 0 )) then
           if ( io_type == 0 ) then
              write(6,*) 'writing data ...   t = ', t, ' rank = ', myrank
              write(filename,990) myrank, n_output
              call output(filename,ix,jx,t,x,y,U,V)
           else
              if( myrank == 0 )  write(6,*) 'writing data ...   t = ', t
              write(filename,980) n_output
              call mpioutput(filename,ix,jx,t,x,y,U,V)
           endif
        endif
        n_output = n_output + 1
        t_output = t_output + dtout
     endif
!    [ end? ]
     if ( t >= tend )  exit
     if ( k >= loop_max ) then
        if( myrank == 0 )  write(6,*) 'max loop'
        exit
     endif
!   -----------------  
!    CFL condition
     call set_dt(U,V,ch,dt,dx,cfl,ix,jx)
     call mpi_iallreduce(ch,chg,1,mpi_real8,mpi_max,cart2d%comm,mreq(1),merr)
!     call mpi_iallreduce(mpi_in_place,ch,1,mpi_real8,mpi_max,cart2d%comm,mreq(1),merr)
     call set_dt2(Rm1,dt,dx,cfl)
     call mpi_iallreduce(dt,dtg,1,mpi_real8,mpi_min,cart2d%comm,mreq(2),merr)
!     call mpi_iallreduce(mpi_in_place,dt,1,mpi_real8,mpi_min,cart2d%comm,mreq(2),merr)
     call mpi_waitall(2,mreq,mpi_statuses_ignore,merr)
     dt = dtg; ch = chg

!    GLM solver for the first half timestep
!    This should be done after set_dt()
     call glm_ss(U,ch,0.5d0*dt,ix,jx)

!    Slope limiters on primitive variables
!     write(6,*) 'V --> VL, VR (F)'
     call limiter(V(1,1,vx),VL(1,1,vx),VR(1,1,vx),ix,jx,1,lm_type)
     call limiter(V(1,1,vy),VL(1,1,vy),VR(1,1,vy),ix,jx,1,lm_type)
     call limiter(V(1,1,vz),VL(1,1,vz),VR(1,1,vz),ix,jx,1,lm_type)
     call limiter(V(1,1,pr),VL(1,1,pr),VR(1,1,pr),ix,jx,1,lm_type)
     call limiter(U(1,1,ro),VL(1,1,ro),VR(1,1,ro),ix,jx,1,lm_type)
     call limiter(U(1,1,bx),VL(1,1,bx),VR(1,1,bx),ix,jx,1,lm_type)
     call limiter(U(1,1,by),VL(1,1,by),VR(1,1,by),ix,jx,1,lm_type)
     call limiter(U(1,1,bz),VL(1,1,bz),VR(1,1,bz),ix,jx,1,lm_type)
     call limiter(U(1,1,ps),VL(1,1,ps),VR(1,1,ps),ix,jx,1,lm_type)
!    fix VL/VR for periodic bc (F)
!     write(6,*) 'fix VL/VR at MPI boundary'
     call parallel_exchange2(VL,VR,ix,jx,1)
!    Numerical flux in the X direction (F)
!     write(6,*) 'VL, VR --> F'
     call flux_solver(F,VL,VR,ix,jx,1,flux_type)
     call flux_glm(F,VL,VR,ch,ix,jx,1)
     call flux_resistive(F,U,VL,VR,EF,dx,ix,jx,1)

!    Slope limiters on primitive variables
!     write(6,*) 'V --> VL, VR (G)'
     call limiter(V(1,1,vx),VL(1,1,vx),VR(1,1,vx),ix,jx,2,lm_type)
     call limiter(V(1,1,vy),VL(1,1,vy),VR(1,1,vy),ix,jx,2,lm_type)
     call limiter(V(1,1,vz),VL(1,1,vz),VR(1,1,vz),ix,jx,2,lm_type)
     call limiter(V(1,1,pr),VL(1,1,pr),VR(1,1,pr),ix,jx,2,lm_type)
     call limiter(U(1,1,ro),VL(1,1,ro),VR(1,1,ro),ix,jx,2,lm_type)
     call limiter(U(1,1,bx),VL(1,1,bx),VR(1,1,bx),ix,jx,2,lm_type)
     call limiter(U(1,1,by),VL(1,1,by),VR(1,1,by),ix,jx,2,lm_type)
     call limiter(U(1,1,bz),VL(1,1,bz),VR(1,1,bz),ix,jx,2,lm_type)
     call limiter(U(1,1,ps),VL(1,1,ps),VR(1,1,ps),ix,jx,2,lm_type)
!    fix VL/VR for periodic bc (G)
!     write(6,*) 'fix VL/VR at MPI boundary'
     call parallel_exchange2(VL,VR,ix,jx,2)
     call mpibc_for_G(VL,VR,ix,jx)
!    Numerical flux in the Y direction (G)
!     write(6,*) 'VL, VR --> G'
     call flux_solver(G,VL,VR,ix,jx,2,flux_type)
     call flux_glm(G,VL,VR,ch,ix,jx,2)
     call flux_resistive(G,U,VL,VR,EG,dx,ix,jx,2)

     if( time_type == 0 ) then
!       write(6,*) 'U* = U + (dt/dx) (F-F)'
        call rk_tvd21(U,U1,F,G,dt,dx,ix,jx)
     elseif( time_type == 1 ) then
!       write(6,*) 'U*(n+1/2) = U + (0.5 dt/dx) (F-F)'
        call rk_std21(U,U1,F,G,dt,dx,ix,jx)
     endif
!    boundary conditions
     call parallel_exchange(U1,ix,jx,1)
     call parallel_exchange(U1,ix,jx,2)
     call mpibc_for_U(U1,ix,jx)
!     write(6,*) 'U* --> V'
     call u2v(U1,V,ix,jx)

!    Slope limiters on primitive variables
!     write(6,*) 'V --> VL, VR (F)'
     call limiter(V(1,1,vx),VL(1,1,vx),VR(1,1,vx),ix,jx,1,lm_type)
     call limiter(V(1,1,vy),VL(1,1,vy),VR(1,1,vy),ix,jx,1,lm_type)
     call limiter(V(1,1,vz),VL(1,1,vz),VR(1,1,vz),ix,jx,1,lm_type)
     call limiter(V(1,1,pr),VL(1,1,pr),VR(1,1,pr),ix,jx,1,lm_type)
     call limiter(U1(1,1,ro),VL(1,1,ro),VR(1,1,ro),ix,jx,1,lm_type)
     call limiter(U1(1,1,bx),VL(1,1,bx),VR(1,1,bx),ix,jx,1,lm_type)
     call limiter(U1(1,1,by),VL(1,1,by),VR(1,1,by),ix,jx,1,lm_type)
     call limiter(U1(1,1,bz),VL(1,1,bz),VR(1,1,bz),ix,jx,1,lm_type)
     call limiter(U1(1,1,ps),VL(1,1,ps),VR(1,1,ps),ix,jx,1,lm_type)
!    fix VL/VR for periodic bc (F)
!     write(6,*) 'fix VL/VR at MPI boundary'
     call parallel_exchange2(VL,VR,ix,jx,1)
!    Numerical flux in the X direction (F)
!     write(6,*) 'VL, VR --> F'
     call flux_solver(F,VL,VR,ix,jx,1,flux_type)
     call flux_glm(F,VL,VR,ch,ix,jx,1)
     call flux_resistive(F,U1,VL,VR,EF,dx,ix,jx,1)

!    Slope limiters on primitive variables
!     write(6,*) 'V --> VL, VR (G)'
     call limiter(V(1,1,vx),VL(1,1,vx),VR(1,1,vx),ix,jx,2,lm_type)
     call limiter(V(1,1,vy),VL(1,1,vy),VR(1,1,vy),ix,jx,2,lm_type)
     call limiter(V(1,1,vz),VL(1,1,vz),VR(1,1,vz),ix,jx,2,lm_type)
     call limiter(V(1,1,pr),VL(1,1,pr),VR(1,1,pr),ix,jx,2,lm_type)
     call limiter(U1(1,1,ro),VL(1,1,ro),VR(1,1,ro),ix,jx,2,lm_type)
     call limiter(U1(1,1,bx),VL(1,1,bx),VR(1,1,bx),ix,jx,2,lm_type)
     call limiter(U1(1,1,by),VL(1,1,by),VR(1,1,by),ix,jx,2,lm_type)
     call limiter(U1(1,1,bz),VL(1,1,bz),VR(1,1,bz),ix,jx,2,lm_type)
     call limiter(U1(1,1,ps),VL(1,1,ps),VR(1,1,ps),ix,jx,2,lm_type)
!    fix VL/VR for periodic bc (G)
!     write(6,*) 'fix VL/VR at MPI boundary'
     call parallel_exchange2(VL,VR,ix,jx,2)
     call mpibc_for_G(VL,VR,ix,jx)
!    Numerical flux in the Y direction (G)
!     write(6,*) 'VL, VR --> G'
     call flux_solver(G,VL,VR,ix,jx,2,flux_type)
     call flux_glm(G,VL,VR,ch,ix,jx,2)
     call flux_resistive(G,U1,VL,VR,EG,dx,ix,jx,2)

     if( time_type == 0 ) then
!       write(6,*) 'U_new = 0.5( U_old + U* + F dt )'
        call rk_tvd22(U,U1,F,G,dt,dx,ix,jx)
     elseif( time_type == 1 ) then
!       write(6,*) 'U_new = U + (dt/dx) (F-F) (n+1/2)'
        call rk_std22(U,F,G,dt,dx,ix,jx)
     endif

!    GLM solver for the second half timestep
     call glm_ss(U,ch,0.5d0*dt,ix,jx)

!    boundary conditions
     call parallel_exchange(U,ix,jx,1)
     call parallel_exchange(U,ix,jx,2)
     call mpibc_for_U(U,ix,jx)
     t=t+dt

  enddo
!-----------------------------------------------------------------------

  call mpi_finalize(merr)
  if( myrank == 0 )  write(6,*) '== end =='

980 format ('data/field-',i5.5,'.dat')
990 format ('data/field-rank',i4.4,'-',i5.5,'.dat')
!981 format ('data/field-',i5.5,'.dat.restart')
!991 format ('data/field-rank',i4.4,'-',i5.5,'.dat.restart')

end program main
!-----------------------------------------------------------------------
