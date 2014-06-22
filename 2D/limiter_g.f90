subroutine limiter_g(wk,wD,wU,ix,jx,type)
!-----------------------------------------------------------------------
!     Slope limiters
!-----------------------------------------------------------------------
!     2010/01/22  S. Zenitani  2nd order limiters (minmod, MC)
!     2010/05/14  S. Zenitani  added van Leer limiter
!     2012/07/05  S. Zenitani  added Koren limiter (bug fixed)
!-----------------------------------------------------------------------
  implicit none
  integer, intent(in)  :: ix, jx
  real(8), intent(in)  :: wk(ix,jx)            ! work array (input)
  real(8), intent(out) :: wD(ix,jx), wU(ix,jx) ! interpolated values (output)
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
     do j=1,jx-1
        wD(:,j) = wk(:,j)
        wU(:,j) = wk(:,j+1)
     enddo

! 2nd order
! minmod limiter
  elseif( type .eq. 1 ) then

     do j=1,1
        wD(:,j)   = wk(:,j)
     enddo
     do j=2,jx-1
        do i=1,ix
           gA =     ( wk(i,j)  -wk(i,j-1))
           gB =     ( wk(i,j+1)-wk(i,j)  )
           if( gA*gB .le. 0 ) then
              grad = 0.d0
           else
              if( gA .gt. 0 ) then
                 grad = 0.5d0 * min(gA,gB)
!                 grad = min(gA,gB)
              else
                 grad = 0.5d0 * max(gA,gB)
!                 grad = max(gA,gB)
              endif
           endif
!           grad = sign(0.5d0,gA) * max( 0.d0, min(abs(gA),sign(1.d0,gA)*gB) )
           wU(i,j-1) = wk(i,j) - grad
           wD(i,j)   = wk(i,j) + grad
!           wU(i,j-1) = wk(i,j) - 0.5d0 * grad
!           wD(i,j)   = wk(i,j) + 0.5d0 * grad
        enddo
     enddo
     do j=jx,jx
        wU(:,j-1) = wk(:,j)
     enddo
         
! 2nd order
! Monotonized central (MC) limiter
  elseif( type .eq. 2 ) then
     
     do j=1,1
        wD(:,j)   = wk(:,j)
     enddo
     do j=2,jx-1
        do i=1,ix
           gA =       ( wk(i,j)  -wk(i,j-1))
           gB =       ( wk(i,j+1)-wk(i,j)  )
           gC = 0.25d0*( wk(i,j+1)-wk(i,j-1))
!           gC = 0.5d0*( wk(i,j+1)-wk(i,j-1))
           if( gA*gB .le. 0 ) then
              grad = 0.d0
           else
              if( gA .gt. 0 ) then
                 grad = min(gA,gB,gC)
!                 grad = min(2*gA,2*gB,gC)
              else
                 grad = max(gA,gB,gC)
!                 grad = max(2*gA,2*gB,gC)
              endif
           endif
           wU(i,j-1) = wk(i,j) - grad
           wD(i,j)   = wk(i,j) + grad
!           wU(i,j-1) = wk(i,j) - 0.5d0 * grad
!           wD(i,j)   = wk(i,j) + 0.5d0 * grad
        enddo
     enddo
     do j=jx,jx
        wU(:,j-1) = wk(:,j)
     enddo
     
! 2nd order
! van Leer limiter
  elseif( type .eq. 3 ) then

     do j=1,1
        wD(:,j)   = wk(:,j)
     enddo
     do j=2,jx-1
        do i=1,ix
           gA =     ( wk(i,j)  -wk(i,j-1))
           gB =     ( wk(i,j+1)-wk(i,j)  )
           if( gA*gB .le. 0 ) then
              grad = 0.d0
           else
              grad = gA*gB / ( gA + gB )
!              grad = 2*gA*gB / ( gA + gB )
           endif
           wU(i,j-1) = wk(i,j) - grad
           wD(i,j)   = wk(i,j) + grad
!           wU(i,j-1) = wk(i,j) - 0.5d0 * grad
!           wD(i,j)   = wk(i,j) + 0.5d0 * grad
        enddo
     enddo
     do j=jx,jx
        wU(:,j-1) = wk(:,j)
     enddo

! 2nd (3rd?) order
! Koren limiter
  elseif( type .eq. 4 ) then

     do j=1,1
        wD(:,j)   = wk(:,j)
     enddo
     do j=2,jx-1
        do i=1,ix
           gA =     ( wk(i,j)  -wk(i,j-1))
           gB =     ( wk(i,j+1)-wk(i,j)  )
           if( gA*gB .le. 0 ) then
              wU(i,j-1) = wk(i,j)
              wD(i,j)   = wk(i,j)
           else

              gC = f1*( 2*gA+gB )
              if( gA .gt. 0 ) then
                 grad = min(gA,gB,gC)
              else
                 grad = max(gA,gB,gC)
              endif
              wU(i,j-1) = wk(i,j) - grad

              gC = f1*( gA+2*gB )
              if( gA .gt. 0 ) then
                 grad = min(gA,gB,gC)
              else
                 grad = max(gA,gB,gC)
              endif
              wD(i,j)   = wk(i,j) + grad
           endif
        enddo
     enddo
     do j=jx,jx
        wU(:,j-1) = wk(:,j)
     enddo
         
  else
     write(6,*) 'unknown limiter type'
     stop
     
  endif
  
  return
end subroutine limiter_g
