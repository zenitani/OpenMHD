subroutine limiter_f(wk,wL,wR,ix,jx,type)
!-----------------------------------------------------------------------
!     Slope limiters
!-----------------------------------------------------------------------
!     2010/01/22  S. Zenitani  2nd order limiters (minmod, MC)
!     2010/05/14  S. Zenitani  added van Leer limiter
!     2012/07/05  S. Zenitani  added Koren limiter (bug fixed)
!     2015/07/19  S. Zenitani  removed if-statements from minmod/MC limiters
!-----------------------------------------------------------------------
  implicit none
  integer, intent(in)  :: ix, jx
  real(8), intent(in)  :: wk(ix,jx)            ! work array (input)
  real(8), intent(out) :: wL(ix,jx), wR(ix,jx) ! interpolated values (output)
!-----------------------------------------------------------------------
! 0 (flat), 1 (minmod), 2 (MC), 3 (van Leer), 4 (Koren)
  integer, intent(in)  :: type
!-----------------------------------------------------------------------
  integer :: i, j
  real(8) :: gA, gB, gC      ! average, left, and right gradients
  real(8) :: grad
  real(8), parameter :: f1 = 1.d0 / 6.d0

! 1st order
  if( type .eq. 0 ) then

     do j=1,jx
        do i=1,ix-1
           wL(i,j) = wk(i,j)
           wR(i,j) = wk(i+1,j)
        enddo
     enddo

! 2nd order
! minmod limiter
  elseif( type .eq. 1 ) then

     do i=1,1
        wL(i,:)   = wk(i,:)
     enddo
     do j=1,jx
        do i=2,ix-1
           gA =     ( wk(i,j)  -wk(i-1,j))
           gB =     ( wk(i+1,j)-wk(i,j)  )
!           if( gA*gB .le. 0 ) then
!              grad = 0.d0
!           else
!              if( gA .gt. 0 ) then
!                 grad = min(gA,gB)
!              else
!                 grad = max(gA,gB)
!              endif
!           endif
!           wR(i-1,j) = wk(i,j) - 0.5d0 * grad
!           wL(i,j)   = wk(i,j) + 0.5d0 * grad
           grad = (sign(0.25d0,gA)+sign(0.25d0,gB))*min(abs(gA),abs(gB))
           wR(i-1,j) = wk(i,j) - grad
           wL(i,j)   = wk(i,j) + grad
        enddo
     enddo
     do i=ix,ix
        wR(i-1,:) = wk(i,:)
     enddo
         
! 2nd order
! Monotonized central (MC) limiter
  elseif( type .eq. 2 ) then
     
     do i=1,1
        wL(i,:)   = wk(i,:)
     enddo
     do j=1,jx
        do i=2,ix-1
           gA =       ( wk(i,j)  -wk(i-1,j))
           gB =       ( wk(i+1,j)-wk(i,j)  )
           gC = 0.25d0*( wk(i+1,j)-wk(i-1,j))
!           gC = 0.5d0*( wk(i+1,j)-wk(i-1,j))
!           if( gA*gB .le. 0 ) then
!              grad = 0.d0
!           else
!              if( gA .gt. 0 ) then
!                 grad = min(gA,gB,gC)
!!                 grad = min(2*gA,2*gB,gC)
!              else
!                 grad = max(gA,gB,gC)
!!                 grad = max(2*gA,2*gB,gC)
!              endif
!           endif
           grad = (sign(0.5d0,gA)+sign(0.5d0,gB))*min(abs(gA),abs(gB),abs(gC))
           wR(i-1,j) = wk(i,j) - grad
           wL(i,j)   = wk(i,j) + grad
!           wR(i-1,j) = wk(i,j) - 0.5d0 * grad
!           wL(i,j)   = wk(i,j) + 0.5d0 * grad
        enddo
     enddo
     do i=ix,ix
        wR(i-1,:) = wk(i,:)
     enddo
     
! 2nd order
! van Leer limiter
  elseif( type .eq. 3 ) then

     do i=1,1
        wL(i,:)   = wk(i,:)
     enddo
     do j=1,jx
        do i=2,ix-1
           gA =     ( wk(i,j)  -wk(i-1,j))
           gB =     ( wk(i+1,j)-wk(i,j)  )
           if( gA*gB .le. 0 ) then
              grad = 0.d0
           else
              grad = gA*gB / ( gA + gB )
!              grad = 2*gA*gB / ( gA + gB )
           endif
           wR(i-1,j) = wk(i,j) - grad
           wL(i,j)   = wk(i,j) + grad
!           wR(i-1,j) = wk(i,j) - 0.5d0 * grad
!           wL(i,j)   = wk(i,j) + 0.5d0 * grad
        enddo
     enddo
     do i=ix,ix
        wR(i-1,:) = wk(i,:)
     enddo
         
! 2nd (3rd?) order
! Koren limiter
  elseif( type .eq. 4 ) then

     do i=1,1
        wL(i,:)   = wk(i,:)
     enddo
     do j=1,jx
        do i=2,ix-1
           gA =       ( wk(i,j)  -wk(i-1,j))
           gB =       ( wk(i+1,j)-wk(i,j)  )

           if( gA*gB .le. 0 ) then
              wR(i-1,j) = wk(i,j)
              wL(i,j)   = wk(i,j)
           else

              gC = f1*( 2*gA+gB )
              if( gA .gt. 0 ) then
                 grad = min(gA,gB,gC)
              else
                 grad = max(gA,gB,gC)
              endif
              wR(i-1,j) = wk(i,j) - grad

              gC = f1*( gA+2*gB )
              if( gA .gt. 0 ) then
                 grad = min(gA,gB,gC)
              else
                 grad = max(gA,gB,gC)
              endif
              wL(i,j)   = wk(i,j) + grad

           endif
        enddo
     enddo
     do i=ix,ix
        wR(i-1,:) = wk(i,:)
     enddo
         
  else
     write(6,*) 'unknown limiter type'
     stop
     
  endif
  
  return
end subroutine limiter_f
