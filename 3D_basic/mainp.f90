program main
!-----------------------------------------------------------------------
!     OpenMHD  Riemann solver (parallel version)
!-----------------------------------------------------------------------
!     2010/09/18  S. Zenitani  2D LLF/HLL/HLLC-G/HLLD solver
!     2015/04/09  S. Zenitani  MPI-IO
!     2018/05/02  S. Zenitani  parallel module (MPI-3 version)
!     2023/12/25  S. Zenitani  3-D version
!-----------------------------------------------------------------------
  use parallel
  implicit none
  include 'param.h'
  integer, parameter :: ix =  50 + 2
  integer, parameter :: jx =  50 + 2
  integer, parameter :: kx =  50 + 2
  integer, parameter :: mpi_nums(3)       = (/2, 2, 4/)  ! MPI numbers
  logical, parameter :: bc_periodicity(3) = (/.true., .true., .true./)
  integer, parameter :: loop_max = 200000
  real(8), parameter :: tend  = 1.01d0
  real(8), parameter :: dtout = 0.1d0 ! output interval
  real(8), parameter :: cfl   = 0.3d0 ! time step
! If non-zero, restart from a previous file. If negative, find a restart file backword in time.
  integer :: n_start = 0
! Slope limiter  (0: flat, 1: minmod, 2: MC, 3: van Leer, 4: Koren)
  integer, parameter :: lm_type   = 1
! Numerical flux (0: LLF, 1: HLL, 2: HLLC, 3: HLLD)
  integer, parameter :: flux_type = 3
! Time marching  (0: TVD RK2, 1: RK2)
  integer, parameter :: time_type = 0
! File I/O  (0: Standard, 1: MPI-IO)
  integer, parameter :: io_type   = 1
!-----------------------------------------------------------------------
! See also modelp.f90
!-----------------------------------------------------------------------
  integer :: k
  integer :: n_loop,n_output
  real(8) :: t, dt, t_output
  real(8) :: ch
  character*256 :: filename
  integer :: merr, myrank, mreq(2)   ! for MPI
!-----------------------------------------------------------------------
  real(8) :: x(ix), y(jx), z(kx), dx
  real(8) :: U(ix,jx,kx,var1)  ! conserved variables (U)
  real(8) :: U1(ix,jx,kx,var1) ! conserved variables: medium state (U*)
  real(8) :: V(ix,jx,kx,var2)  ! primitive variables (V)
  real(8) :: VL(ix,jx,kx,var1), VR(ix,jx,kx,var1) ! interpolated states
  real(8) :: F(ix,jx,kx,var1), G(ix,jx,kx,var1), H(ix,jx,kx,var1)   ! numerical flux (F,G,H)
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! for MPI
  call parallel_init(mpi_nums,bc_periodicity,ix,jx,kx)
  myrank = ranks%myrank
!-----------------------------------------------------------------------

  t    =  0.d0
  dt   =  0.d0
  call modelp(U,V,x,y,z,dx,ix,jx,kx)
! boundary conditions
  call parallel_exchange(U,ix,jx,kx,1) ! in the x direction
  call parallel_exchange(U,ix,jx,kx,2) ! in the y direction
  call parallel_exchange(U,ix,jx,kx,3) ! in the z direction
  call set_dt(U,V,ch,dt,dx,cfl,ix,jx,kx)

  call mpi_iallreduce(mpi_in_place,ch,1,mpi_real8,mpi_max,cart3d%comm,mreq(1),merr)
  call mpi_iallreduce(mpi_in_place,dt,1,mpi_real8,mpi_min,cart3d%comm,mreq(2),merr)
  call mpi_waitall(2,mreq,mpi_statuses_ignore,merr)
  call mpi_barrier(cart3d%comm,merr)

  if ( dt > dtout ) then
     write(6,*) 'error: ', dt, '>', dtout
     stop
  endif
  if( myrank == 0 ) then
     write(6,*) '[Params]'
     write(6,997) version, ranks%size, cart3d%sizes
     write(6,998) dt, dtout, cart3d%sizes(1)*(ix-2)+2, ix, cart3d%sizes(2)*(jx-2)+2, jx, cart3d%sizes(3)*(kx-2)+2, kx
     write(6,999) lm_type, flux_type, time_type
997  format ('Code version: ', i8, '  MPI node # : ', i5,' (',i4,' x ',i4,' x ',i4,')' )
998  format (' dt:', 1p, e10.3, ' dtout:', 1p, e10.3, ' grids:',i6,' (',i5,') x ',i6,' (',i5,') x ',i6,' (',i5,') ')
999  format (' limiter: ', i1, '  flux: ', i1, '  time-marching: ', i1 )
     write(6,*) '== start =='
  endif

  ! If n_start is negative, look for a latest restart file.
  if ( n_start < 0 ) then
     if ( myrank == 0 ) then
        do k = floor(tend/dtout),0,-1
           n_start = k
           if ( io_type == 0 )  write(filename,990) myrank, n_start
           if ( io_type /= 0 )  write(filename,980) n_start
           open(15,file=filename,form='unformatted',access='stream',status='old',err=100)
           close(15)
           exit
100        continue
        enddo
     endif
     call mpi_bcast(n_start,1,mpi_integer,0,cart3d%comm,merr)
  endif
  call mpi_barrier(cart3d%comm,merr)

  ! If n_start is non-zero, restart from a previous file.
  if ( n_start == 0 ) then
     n_output = 0
     t_output = -dt/3.d0
  else
     if ( io_type == 0 ) then
        write(6,*) 'reading data ...   rank = ', myrank
        write(filename,990) myrank, n_start
        call fileio_input(filename,ix,jx,kx,t,x,y,z,U)
     else
        if( myrank == 0 )  write(6,*) 'reading data ...'
        write(filename,980) n_start
        call mpiio_input(filename,ix,jx,kx,t,x,y,z,U)
     endif
     n_output = n_start + 1
     t_output = t + dtout
  endif

!-----------------------------------------------------------------------
  do n_loop=1,loop_max

     if( myrank == 0 ) then
        write(6,*) ' t = ', t
     endif
!    Recovering primitive variables
!     write(6,*) 'U --> V'
     call u2v(U,V,ix,jx,kx)
!   -----------------
!    [ output ]
     if ( t >= t_output ) then
        if (( n_loop > 1 ).or.( n_start == 0 )) then
           if ( io_type == 0 ) then
              write(6,*) 'writing data ...   t = ', t, ' rank = ', myrank
              write(filename,990) myrank, n_output
              call fileio_output(filename,ix,jx,kx,t,x,y,z,U,V)
           else
              if( myrank == 0 )  write(6,*) 'writing data ...   t = ', t
              write(filename,980) n_output
              call mpiio_output(filename,ix,jx,kx,t,x,y,z,U,V)
           endif
        endif
        n_output = n_output + 1
        t_output = t_output + dtout
        call mpi_barrier(cart3d%comm,merr)
     endif
!    [ end? ]
     if ( t >= tend )  exit
     if ( n_loop >= loop_max ) then
        if( myrank == 0 )  write(6,*) 'max loop'
        exit
     endif
!   -----------------
!    CFL condition
     call set_dt(U,V,ch,dt,dx,cfl,ix,jx,kx)
     call mpi_iallreduce(mpi_in_place,ch,1,mpi_real8,mpi_max,cart3d%comm,mreq(1),merr)
     call mpi_iallreduce(mpi_in_place,dt,1,mpi_real8,mpi_min,cart3d%comm,mreq(2),merr)
     call mpi_waitall(2,mreq,mpi_statuses_ignore,merr)

!    GLM solver for the first half timestep
!    This should be done after set_dt()
!     write(6,*) 'U --> SS'
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
!    Interpolated primitive variables at MPI boundaries
!     write(6,*) 'adjusting VL/VR (F)'
     call parallel_exchange2(VL,VR,ix,jx,kx,1) ! in the x direction
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
!    Interpolated primitive variables at MPI boundaries
!     write(6,*) 'adjusting VL/VR (G)'
     call parallel_exchange2(VL,VR,ix,jx,kx,2) ! in the y direction
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
!    Interpolated primitive variables at MPI boundaries
!     write(6,*) 'adjusting VL/VR (H)'
     call parallel_exchange2(VL,VR,ix,jx,kx,3) ! in the z direction
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
!    boundary conditions
     call parallel_exchange(U1,ix,jx,kx,1) ! in the x direction
     call parallel_exchange(U1,ix,jx,kx,2) ! in the y direction
     call parallel_exchange(U1,ix,jx,kx,3) ! in the z direction
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
!    Interpolated primitive variables at MPI boundaries
!     write(6,*) 'adjusting VL/VR (F)'
     call parallel_exchange2(VL,VR,ix,jx,kx,1) ! in the x direction
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
!    Interpolated primitive variables at MPI boundaries
!     write(6,*) 'adjusting VL/VR (G)'
     call parallel_exchange2(VL,VR,ix,jx,kx,2) ! in the y direction
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
!    Interpolated primitive variables at MPI boundaries
!     write(6,*) 'adjusting VL/VR (H)'
     call parallel_exchange2(VL,VR,ix,jx,kx,3) ! in the z direction
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

!    boundary conditions
     call parallel_exchange(U,ix,jx,kx,1) ! in the x direction
     call parallel_exchange(U,ix,jx,kx,2) ! in the y direction
     call parallel_exchange(U,ix,jx,kx,3) ! in the z direction

!    GLM solver for the second half timestep
     call glm_ss2(U,ch,dt,ix,jx,kx)

     t=t+dt
  enddo
!-----------------------------------------------------------------------

  call parallel_finalize()
  if( myrank == 0 )  write(6,*) '== end =='


980 format ('data/field-',i5.5,'.dat')
990 format ('data/field-rank',i5.5,'-',i5.5,'.dat')
!981 format ('data/field-',i5.5,'.dat.restart')
!991 format ('data/field-rank',i5.5,'-',i5.5,'.dat.restart')

end program main
!-----------------------------------------------------------------------
