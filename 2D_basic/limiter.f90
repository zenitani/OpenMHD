subroutine limiter(wk,wL,wR,ix,jx,dir,type)
!-----------------------------------------------------------------------
!     This routine interpolates variables with slope limiters
!-----------------------------------------------------------------------
!     2010/01/22  S. Zenitani  2nd order limiters (minmod, MC)
!     2010/05/14  S. Zenitani  added van Leer limiter
!     2012/07/05  S. Zenitani  added Koren limiter (bug fixed)
!     2016/09/07  S. Zenitani  X and Y directions
!-----------------------------------------------------------------------
  implicit none
!-----------------------------------------------------------------------
  integer, intent(in)  :: ix, jx               ! size of arrays [input]
  real(8), intent(in)  :: wk(ix,jx)            ! work array [input]
  real(8), intent(out) :: wL(ix,jx), wR(ix,jx) ! interpolated values [output]
! direction [input]: 1 (X), 2 (Y), 3 (Z)
  integer, intent(in)  :: dir
! slope limiters [input]: 0 (flat), 1 (minmod), 2 (MC), 3 (van Leer), 4 (Koren)
  integer, intent(in)  :: type
!-----------------------------------------------------------------------
  integer :: i, j, is=0, ie=0, js=0, je=0
  real(8) :: gA, gB, gC, grad                  ! gradients
  real(8), parameter :: f1 = 1.d0 / 6.d0
!-----------------------------------------------------------------------

!  wL = 0.d0
!  wR = 0.d0
!  wL = -999999999.d0
!  wR = -999999999.d0

  select case(dir)
  case(1)
     js = min(2,jx)
     je = max(1,jx-1)
!     if( ix < 3 ) then
!        write(6,*) 'Slope limiter: Too short in X: ', ix
!        stop
!     endif
  case(2)
     is = min(2,ix)
     ie = max(1,ix-1)
!     if( jx < 3 ) then
!        write(6,*) 'Slope limiter: Too short in Y: ', jx
!        stop
!     endif
  case(3)
  end select


  select case(type)
  !-----------------------------------------------------------------------
  !  1st order
  !  No OpenMP optimization - just for education
  !-----------------------------------------------------------------------
  case(0)

     select case(dir)
     case(1)
        do j=js,je
           do i=1,ix-1
              wL(i,j) = wk(i,j)
              wR(i,j) = wk(i+1,j)
           enddo
        enddo

     case(2)
        do j=1,jx-1
           wL(is:ie,j) = wk(is:ie,j)
           wR(is:ie,j) = wk(is:ie,j+1)
        enddo

     end select

  !-----------------------------------------------------------------------
  !  2nd order:  minmod limiter
  !-----------------------------------------------------------------------
  case(1)

     select case(dir)
     case(1)
!$omp parallel do private(i,j,gA,gB,grad)
        do j=js,je
           wL(1,j) = wk(1,j)
           do i=2,ix-1
              gA = ( wk(i,j)  -wk(i-1,j))
              gB = ( wk(i+1,j)-wk(i,j)  )
!             if( gA*gB <= 0 ) then
!                grad = 0.d0
!             else
!                if( gA > 0 ) then
!                   grad = min(gA,gB)
!                else
!                   grad = max(gA,gB)
!                endif
!             endif
!             wR(i-1,j) = wk(i,j) - 0.5d0 * grad
!             wL(i,j)   = wk(i,j) + 0.5d0 * grad
!             grad = sign(0.5d0,gA) * max( 0.d0,min(abs(gA),sign(1.d0,gA)*gB) )
              grad = (sign(0.25d0,gA)+sign(0.25d0,gB))*min(abs(gA),abs(gB))
!             grad = 0.5d0*((gA+gB)-min(0.d0,gA,gB)-max(0.d0,gA,gB))
              wR(i-1,j) = wk(i,j) - grad
              wL(i,j)   = wk(i,j) + grad
           enddo
           wR(ix-1,j) = wk(ix,j)
        enddo
!$omp end parallel do

     case(2)
!$omp parallel private(i,j,gA,gB,grad)
!$omp workshare
        wL(is:ie,1) = wk(is:ie,1)
!$omp end workshare
!$omp do
        do j=2,jx-1
           do i=is,ie
              gA = ( wk(i,j)  -wk(i,j-1))
              gB = ( wk(i,j+1)-wk(i,j)  )
              grad = (sign(0.25d0,gA)+sign(0.25d0,gB))*min(abs(gA),abs(gB))
              wR(i,j-1) = wk(i,j) - grad
              wL(i,j)   = wk(i,j) + grad
           enddo
        enddo
!$omp end do
!$omp workshare
        wR(is:ie,jx-1) = wk(is:ie,jx)
!$omp end workshare
!$omp end parallel

     end select

  !-----------------------------------------------------------------------
  !  2nd order:  Monotonized central (MC) limiter
  !-----------------------------------------------------------------------
  case(2)
     
     select case(dir)
     case(1)
!$omp parallel do private(i,j,gA,gB,gC,grad)
        do j=js,je
           wL(1,j) = wk(1,j)
           do i=2,ix-1
              gA =        ( wk(i,j)  -wk(i-1,j))
              gB =        ( wk(i+1,j)-wk(i,j)  )
              gC = 0.25d0*( wk(i+1,j)-wk(i-1,j))
!             gC = 0.5d0*( wk(i+1,j)-wk(i-1,j))
!             if( gA*gB <= 0 ) then
!                grad = 0.d0
!             else
!                if( gA > 0 ) then
!                   grad = min(gA,gB,gC)
!!                   grad = min(2*gA,2*gB,gC)
!                else
!                   grad = max(gA,gB,gC)
!!                   grad = max(2*gA,2*gB,gC)
!                endif
!             endif
              grad = (sign(0.5d0,gA)+sign(0.5d0,gB))*min(abs(gA),abs(gB),abs(gC))
              wR(i-1,j) = wk(i,j) - grad
              wL(i,j)   = wk(i,j) + grad
!             wR(i-1,j) = wk(i,j) - 0.5d0 * grad
!             wL(i,j)   = wk(i,j) + 0.5d0 * grad
           enddo
           wR(ix-1,j) = wk(ix,j)
        enddo
!$omp end parallel do

     case(2)
!$omp parallel private(i,j,gA,gB,gC,grad)
!$omp workshare
        wL(is:ie,1)   = wk(is:ie,1)
!$omp end workshare
!$omp do
        do j=2,jx-1
           do i=is,ie
              gA =       ( wk(i,j)  -wk(i,j-1))
              gB =       ( wk(i,j+1)-wk(i,j)  )
              gC = 0.25d0*( wk(i,j+1)-wk(i,j-1))
              grad = (sign(0.5d0,gA)+sign(0.5d0,gB))*min(abs(gA),abs(gB),abs(gC))
              wR(i,j-1) = wk(i,j) - grad
              wL(i,j)   = wk(i,j) + grad
           enddo
        enddo
!$omp end do
!$omp workshare
        wR(is:ie,jx-1) = wk(is:ie,jx)
!$omp end workshare
!$omp end parallel

     end select
     
  !-----------------------------------------------------------------------
  !  2nd order:  van Leer limiter
  !-----------------------------------------------------------------------
  case(3)

     select case(dir)
     case(1)
!$omp parallel do private(i,j,gA,gB,grad)
        do j=js,je
           wL(1,j) = wk(1,j)
           do i=2,ix-1
              gA = ( wk(i,j)  -wk(i-1,j))
              gB = ( wk(i+1,j)-wk(i,j)  )
              if( gA*gB <= 0 ) then
                 grad = 0.d0
              else
                 grad = gA*gB / ( gA + gB )
!                grad = 2*gA*gB / ( gA + gB )
              endif
              wR(i-1,j) = wk(i,j) - grad
              wL(i,j)   = wk(i,j) + grad
!             wR(i-1,j) = wk(i,j) - 0.5d0 * grad
!             wL(i,j)   = wk(i,j) + 0.5d0 * grad
           enddo
           wR(ix-1,j) = wk(ix,j)
        enddo
!$omp end parallel do

     case(2)
!$omp parallel private(i,j,gA,gB,grad)
!$omp workshare
        wL(is:ie,1)   = wk(is:ie,1)
!$omp end workshare
!$omp do
        do j=2,jx-1
           do i=is,ie
              gA = ( wk(i,j)  -wk(i,j-1))
              gB = ( wk(i,j+1)-wk(i,j)  )
              if( gA*gB <= 0 ) then
                 grad = 0.d0
              else
                 grad = gA*gB / ( gA + gB )
              endif
              wR(i,j-1) = wk(i,j) - grad
              wL(i,j)   = wk(i,j) + grad
           enddo
        enddo
!$omp end do
!$omp workshare
        wR(is:ie,jx-1) = wk(is:ie,jx)
!$omp end workshare
!$omp end parallel

     end select

  !-----------------------------------------------------------------------
  !  3rd order:  Koren limiter
  !-----------------------------------------------------------------------
  case(4)

     select case(dir)
     case(1)
!$omp parallel do private(i,j,gA,gB,gC,grad)
        do j=js,je
           wL(1,j) = wk(1,j)
           do i=2,ix-1
              gA = ( wk(i,j)  -wk(i-1,j))
              gB = ( wk(i+1,j)-wk(i,j)  )
!             if( gA*gB <= 0 ) then
!                wR(i-1,j) = wk(i,j)
!                wL(i,j)   = wk(i,j)
!             else
              gC = f1*( 2*gA+gB )
              wR(i-1,j) = wk(i,j) - (sign(0.5d0,gA)+sign(0.5d0,gB))*min(abs(gA),abs(gB),abs(gC))
!             if( gA > 0 ) then
!                grad = min(gA,gB,gC)
!             else
!                grad = max(gA,gB,gC)
!             endif
!             wR(i-1,j) = wk(i,j) - grad

              gC = f1*( gA+2*gB )
              wL(i,j) = wk(i,j) + (sign(0.5d0,gA)+sign(0.5d0,gB))*min(abs(gA),abs(gB),abs(gC))
!             if( gA > 0 ) then
!                grad = min(gA,gB,gC)
!             else
!                grad = max(gA,gB,gC)
!             endif
!             wL(i,j)   = wk(i,j) + grad
!             endif
           enddo
           wR(ix-1,j) = wk(ix,j)
        enddo
!$omp end parallel do
         
     case(2)
!$omp parallel private(i,j,gA,gB,gC,grad)
!$omp workshare
        wL(is:ie,1)   = wk(is:ie,1)
!$omp end workshare
!$omp do
        do j=2,jx-1
           do i=is,ie
              gA = ( wk(i,j)  -wk(i,j-1))
              gB = ( wk(i,j+1)-wk(i,j)  )
              gC = f1*( 2*gA+gB )
              wR(i,j-1) = wk(i,j) - (sign(0.5d0,gA)+sign(0.5d0,gB))*min(abs(gA),abs(gB),abs(gC))
              gC = f1*( gA+2*gB )
              wL(i,j) = wk(i,j) + (sign(0.5d0,gA)+sign(0.5d0,gB))*min(abs(gA),abs(gB),abs(gC))
           enddo
        enddo
!$omp end do
!$omp workshare
        wR(is:ie,jx-1) = wk(is:ie,jx)
!$omp end workshare
!$omp end parallel

     end select

  !-----------------------------------------------------------------------
  case default
     write(6,*) 'unknown limiter'
     stop
     
  end select
  !-----------------------------------------------------------------------
  
  return
end subroutine limiter
