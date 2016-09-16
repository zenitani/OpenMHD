subroutine model(U,V,x,y,dx,ix,jx)
  implicit none
  include 'param.h'
  real(8) :: U(ix,jx,var1)
  real(8) :: V(ix,jx,var2)
  real(8) :: x(ix), y(jx), dx
  integer :: ix, jx
  integer :: i, j, izero
  real(8) :: ro0, vx0, vx1, pr0, bn0, bt0
  real(8) :: B2, v2, f1, fluxm, b2m, bm, phix, ptot
! ---------------------------------------------------
  real(8), parameter :: pi  = acos(-1.d0)
! velocity transition layer (width and position)
  real(8), parameter :: xtr =  400.d0
  real(8), parameter :: xc  = 1500.d0
! Alfven fluctuation (width and position)
  real(8), parameter :: xar =  100.d0
  real(8), parameter :: xa  =  600.d0
! ---------------------------------------------------
! HSS parameters
  real(8), parameter :: am    = 10.d0  ! Alfven Mach number
  real(8), parameter :: alpha = 0.5d0  ! ratio of to LSS velocity
  real(8), parameter :: beta =  0.8d0  ! plasma beta
  real(8), parameter :: theta = 60.d0 * pi /180.d0  ! magnetic angle
! ---------------------------------------------------

!  dx=1.d0/dble(ix-2)
  dx=1.d0
!  izero=ix/2
  izero=1
  x(izero)=-dx/2
  do i=izero+1,ix
     x(i) = x(i-1)+dx
  enddo
  do i=izero-1,1,-1
     x(i) = x(i+1)-dx
  enddo
  y(1) = 0.d0

! HSS parameters
  vx0 = am * sqrt(2./gamma/beta)
  ro0 = 1.d0
  pr0 = 1.d0/gamma
! LSS parameters
  vx1 = vx0 * alpha
!  ro1 = ro0 / alpha

  fluxm = ro0*vx0
  ptot  = pr0 * (1.d0 + 1.d0/beta)

  bn0 = sqrt(2./gamma/beta) * cos(theta)
  bt0 = sqrt(2./gamma/beta) * sin(theta)

  j = 1
  do i = 1, ix
     V(i,j,vx) = 0.5*(vx1+vx0) + 0.5*(vx1-vx0)*tanh( (x(i)-xc)/xtr )
     U(i,j,ro) = fluxm / V(i,j,vx)  ! constant mass flux assumption
     U(i,j,bx) = bn0

     phix = 0.5*pi* (1. - tanh( (x(i)-xa)/xar ) )
! constant Mach number (pressure balance) model
     bm  = sqrt( fluxm * V(i,j,vx) )/am
     b2m = sqrt( bm**2 - bn0**2 )    ! b2m

! simple CIR model
!     U(i,j,by) = 0.d0
!     U(i,j,bz) = b2m
! B,p constant model
!     U(i,j,by) = bt0 * sin(phix)
!     U(i,j,bz) = bt0 * cos(phix)
! constant Mach number (pressure balance) model
     U(i,j,by) = b2m * sin(phix)
     U(i,j,bz) = b2m * cos(phix)

     B2 = dot_product( U(i,j,bx:bz), U(i,j,bx:bz) )

! pressure constant model (B is also constant)
!    V(i,j,pr) = pr0
! pressure balance model (MA is constant... B is variable)
     V(i,j,pr) = ptot - 0.5d0 * B2
     if( V(i,j,pr) .lt. 0 ) then
        write(6,*) 'init error'
        stop
     endif
     
! simple CIR model
!     V(i,j,vy) = 0.d0
!     V(i,j,vz) = 0.d0
! include Alfven term
     V(i,j,vy) = -U(i,j,by)/sqrt(U(i,j,ro))
     V(i,j,vz) = -U(i,j,bz)/sqrt(U(i,j,ro))

     f1 = 1.d0 / ( gamma - 1 )
     v2 = dot_product( V(i,j,vx:vz), V(i,j,vx:vz) )
     U(i,j,mx:mz) = U(i,j,ro) * V(i,j,vx:vz)
     U(i,j,en)    = 0.5d0 * ( U(i,j,ro)*v2 + B2 ) + f1*V(i,j,pr)
     U(i,j,ps)    = 0.d0

  enddo
     
  return
end subroutine model
