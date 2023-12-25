subroutine limiter(wk,wL,wR,ix,jx,kx,dir,type)
!-----------------------------------------------------------------------
!     This routine interpolates variables with slope limiters
!-----------------------------------------------------------------------
!     2010/01/22  S. Zenitani  2nd order limiters (minmod, MC)
!     2010/05/14  S. Zenitani  added van Leer limiter
!     2012/07/05  S. Zenitani  added Koren limiter (bug fixed)
!     2016/09/07  S. Zenitani  X and Y directions
!     2023/12/24  S. Zenitani  3D: Z direction
!-----------------------------------------------------------------------
  implicit none
!-----------------------------------------------------------------------
  integer, intent(in)  :: ix, jx, kx                 ! size of arrays [input]
  real(8), intent(in)  :: wk(ix,jx,kx)               ! work array [input]
  real(8), intent(out) :: wL(ix,jx,kx), wR(ix,jx,kx) ! interpolated values [output]
! direction [input]: 1 (X), 2 (Y), 3 (Z)
  integer, intent(in)  :: dir
! slope limiters [input]: 0 (flat), 1 (minmod), 2 (MC), 3 (van Leer), 4 (Koren)
  integer, intent(in)  :: type
!-----------------------------------------------------------------------
  integer :: i, j, k, is=0, ie=0, js=0, je=0, ks=0, ke=0
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
     ks = min(2,kx)
     ke = max(1,kx-1)
!     if( ix < 3 ) then
!        write(6,*) 'Slope limiter: Too short in X: ', ix
!        stop
!     endif
  case(2)
     is = min(2,ix)
     ie = max(1,ix-1)
     ks = min(2,kx)
     ke = max(1,kx-1)
!     if( jx < 3 ) then
!        write(6,*) 'Slope limiter: Too short in Y: ', jx
!        stop
!     endif
  case(3)
     is = min(2,ix)
     ie = max(1,ix-1)
     js = min(2,jx)
     je = max(1,jx-1)
  end select


  select case(type)
  !-----------------------------------------------------------------------
  !  1st order
  !  No OpenMP optimization - just for education
  !-----------------------------------------------------------------------
  case(0)

     select case(dir)
     case(1)
        do k=ks,ke
           do j=js,je
              do i=1,ix-1
                 wL(i,j,k) = wk(i  ,j,k)
                 wR(i,j,k) = wk(i+1,j,k)
              enddo
           enddo
        enddo

     case(2)
        do k=ks,ke
           do j=1,jx-1
              wL(is:ie,j,k) = wk(is:ie,j  ,k)
              wR(is:ie,j,k) = wk(is:ie,j+1,k)
           enddo
        enddo

     case(3)
        do k=1,kx-1
           wL(is:ie,js:je,k) = wk(is:ie,js:je,k  )
           wR(is:ie,js:je,k) = wk(is:ie,js:je,k+1)
        enddo

     end select

  !-----------------------------------------------------------------------
  !  2nd order:  minmod limiter
  !-----------------------------------------------------------------------
  case(1)

     select case(dir)
     case(1)
!$omp parallel do private(i,j,k,gA,gB,grad)
        do k=ks,ke
           do j=js,je
              wL(1,j,k) = wk(1,j,k)
              do i=2,ix-1
                 gA = ( wk(i,j,k)  -wk(i-1,j,k) )
                 gB = ( wk(i+1,j,k)-wk(i  ,j,k) )
                 grad = (sign(0.25d0,gA)+sign(0.25d0,gB))*min(abs(gA),abs(gB))
                 wR(i-1,j,k) = wk(i,j,k) - grad
                 wL(i  ,j,k) = wk(i,j,k) + grad
              enddo
              wR(ix-1,j,k) = wk(ix,j,k)
           enddo
        enddo
!$omp end parallel do

     case(2)
!$omp parallel private(i,j,k,gA,gB,grad)
!$omp do
        do k=ks,ke
           wL(is:ie,1,k) = wk(is:ie,1,k)
        enddo
!$omp end do
!$omp do
        do k=ks,ke
           do j=2,jx-1
              do i=is,ie
                 gA = ( wk(i,j  ,k)-wk(i,j-1,k) )
                 gB = ( wk(i,j+1,k)-wk(i,j  ,k) )
                 grad = (sign(0.25d0,gA)+sign(0.25d0,gB))*min(abs(gA),abs(gB))
                 wR(i,j-1,k) = wk(i,j,k) - grad
                 wL(i,j  ,k) = wk(i,j,k) + grad
              enddo
           enddo
        enddo
!$omp end do
!$omp do
        do k=ks,ke
           wR(is:ie,jx-1,k) = wk(is:ie,jx,k)
        enddo
!$omp end do
!$omp end parallel
        
     case(3)
!$omp parallel private(i,j,k,gA,gB,grad)
!$omp workshare
        wL(is:ie,js:je,1) = wk(is:ie,js:je,1)
!$omp end workshare
!$omp do
        do k=2,kx-1
           do j=js,je
              do i=is,ie
                 gA = ( wk(i,j,k  )-wk(i,j,k-1) )
                 gB = ( wk(i,j,k+1)-wk(i,j,k  ) )
                 grad = (sign(0.25d0,gA)+sign(0.25d0,gB))*min(abs(gA),abs(gB))
                 wR(i,j,k-1) = wk(i,j,k) - grad
                 wL(i,j,k  ) = wk(i,j,k) + grad
              enddo
           enddo
        enddo
!$omp end do
!$omp workshare
        wR(is:ie,js:je,kx-1) = wk(is:ie,js:je,kx)
!$omp end workshare
!$omp end parallel

     end select

  !-----------------------------------------------------------------------
  !  2nd order:  Monotonized central (MC) limiter
  !-----------------------------------------------------------------------
  case(2)
     
     select case(dir)
     case(1)
!$omp parallel do private(i,j,k,gA,gB,gC,grad)
        do k=ks,ke
           do j=js,je
              wL(1,j,k) = wk(1,j,k)
              do i=2,ix-1
                 gA =        ( wk(i  ,j,k)-wk(i-1,j,k) )
                 gB =        ( wk(i+1,j,k)-wk(i  ,j,k) )
                 gC = 0.25d0*( wk(i+1,j,k)-wk(i-1,j,k) )
                 grad = (sign(0.5d0,gA)+sign(0.5d0,gB))*min(abs(gA),abs(gB),abs(gC))
                 wR(i-1,j,k) = wk(i,j,k) - grad
                 wL(i  ,j,k) = wk(i,j,k) + grad
              enddo
              wR(ix-1,j,k) = wk(ix,j,k)
           enddo
        enddo
!$omp end parallel do

     case(2)
!$omp parallel private(i,j,k,gA,gB,gC,grad)
!$omp do
        do k=ks,ke
           wL(is:ie,1,k) = wk(is:ie,1,k)
        enddo
!$omp end do
!$omp do
        do k=ks,ke
           do j=2,jx-1
              do i=is,ie
                 gA =        ( wk(i,j  ,k)-wk(i,j-1,k) )
                 gB =        ( wk(i,j+1,k)-wk(i,j  ,k) )
                 gC = 0.25d0*( wk(i,j+1,k)-wk(i,j-1,k) )
                 grad = (sign(0.5d0,gA)+sign(0.5d0,gB))*min(abs(gA),abs(gB),abs(gC))
                 wR(i,j-1,k) = wk(i,j,k) - grad
                 wL(i,j  ,k) = wk(i,j,k) + grad
              enddo
           enddo
        enddo
!$omp end do
!$omp do
        do k=ks,ke
           wR(is:ie,jx-1,k) = wk(is:ie,jx,k)
        enddo
!$omp end do
!$omp end parallel

     case(3)
!$omp parallel private(i,j,k,gA,gB,gC,grad)
!$omp workshare
        wL(is:ie,js:je,1) = wk(is:ie,js:je,1)
!$omp end workshare
!$omp do
        do k=2,kx-1
           do j=js,je
              do i=is,ie
                 gA =        ( wk(i,j,k  )-wk(i,j,k-1) )
                 gB =        ( wk(i,j,k+1)-wk(i,j,k  ) )
                 gC = 0.25d0*( wk(i,j,k+1)-wk(i,j,k-1) )
                 grad = (sign(0.5d0,gA)+sign(0.5d0,gB))*min(abs(gA),abs(gB),abs(gC))
                 wR(i,j,k-1) = wk(i,j,k) - grad
                 wL(i,j,k  ) = wk(i,j,k) + grad
              enddo
           enddo
        enddo
!$omp end do
!$omp workshare
        wR(is:ie,js:je,kx-1) = wk(is:ie,js:je,kx)
!$omp end workshare
!$omp end parallel

     end select
     
  !-----------------------------------------------------------------------
  !  2nd order:  van Leer limiter
  !-----------------------------------------------------------------------
  case(3)

     select case(dir)
     case(1)
!$omp parallel do private(i,j,k,gA,gB,grad)
        do k=ks,ke
           do j=js,je
              wL(1,j,k) = wk(1,j,k)
              do i=2,ix-1
                 gA = ( wk(i  ,j,k)-wk(i-1,j,k) )
                 gB = ( wk(i+1,j,k)-wk(i  ,j,k) )
                 if( gA*gB <= 0 ) then
                    grad = 0.d0
                 else
                    grad = gA*gB / ( gA + gB )
                 endif
                 wR(i-1,j,k) = wk(i,j,k) - grad
                 wL(i  ,j,k) = wk(i,j,k) + grad
              enddo
              wR(ix-1,j,k) = wk(ix,j,k)
           enddo
        enddo
!$omp end parallel do

     case(2)
!$omp parallel private(i,j,k,gA,gB,grad)
!$omp do
        do k=ks,ke
           wL(is:ie,1,k) = wk(is:ie,1,k)
        enddo
!$omp end do
!$omp do
        do k=ks,ke
           do j=2,jx-1
              do i=is,ie
                 gA = ( wk(i,j  ,k)-wk(i,j-1,k) )
                 gB = ( wk(i,j+1,k)-wk(i,j  ,k) )
                 if( gA*gB <= 0 ) then
                    grad = 0.d0
                 else
                    grad = gA*gB / ( gA + gB )
                 endif
                 wR(i,j-1,k) = wk(i,j,k) - grad
                 wL(i,j  ,k) = wk(i,j,k) + grad
              enddo
           enddo
        enddo
!$omp end do
!$omp do
        do k=ks,ke
           wR(is:ie,jx-1,k) = wk(is:ie,jx,k)
        enddo
!$omp end do
!$omp end parallel

     case(3)
!$omp parallel private(i,j,k,gA,gB,grad)
!$omp workshare
        wL(is:ie,js:je,1)   = wk(is:ie,js:je,1)
!$omp end workshare
!$omp do
        do k=2,kx-1
           do j=js,je
              do i=is,ie
                 gA = ( wk(i,j,k  )-wk(i,j,k-1) )
                 gB = ( wk(i,j,k+1)-wk(i,j,k  ) )
                 if( gA*gB <= 0 ) then
                    grad = 0.d0
                 else
                    grad = gA*gB / ( gA + gB )
                 endif
                 wR(i,j,k-1) = wk(i,j,k) - grad
                 wL(i,j,k  ) = wk(i,j,k) + grad
              enddo
           enddo
        enddo
!$omp end do
!$omp workshare
        wR(is:ie,js:je,kx-1) = wk(is:ie,js:je,kx)
!$omp end workshare
!$omp end parallel

     end select

  !-----------------------------------------------------------------------
  !  3rd order:  Koren limiter
  !-----------------------------------------------------------------------
  case(4)

     select case(dir)
     case(1)
!$omp parallel do private(i,j,k,gA,gB,gC,grad)
        do k=ks,ke
           do j=js,je
              wL(1,j,k) = wk(1,j,k)
              do i=2,ix-1
                 gA = ( wk(i  ,j,k)-wk(i-1,j,k) )
                 gB = ( wk(i+1,j,k)-wk(i  ,j,k) )
                 gC = f1*( 2*gA+gB )
                 wR(i-1,j,k) = wk(i,j,k) - (sign(0.5d0,gA)+sign(0.5d0,gB))*min(abs(gA),abs(gB),abs(gC))
                 gC = f1*( gA+2*gB )
                 wL(i  ,j,k) = wk(i,j,k) + (sign(0.5d0,gA)+sign(0.5d0,gB))*min(abs(gA),abs(gB),abs(gC))
              enddo
              wR(ix-1,j,k) = wk(ix,j,k)
           enddo
        enddo
!$omp end parallel do
         
     case(2)
!$omp parallel private(i,j,k,gA,gB,gC,grad)
!$omp do
        do k=ks,ke
           wL(is:ie,1,k) = wk(is:ie,1,k)
        enddo
!$omp end do
!$omp do
        do k=ks,ke
           do j=2,jx-1
              do i=is,ie
                 gA = ( wk(i,j  ,k)-wk(i,j-1,k) )
                 gB = ( wk(i,j+1,k)-wk(i,j  ,k) )
                 gC = f1*( 2*gA+gB )
                 wR(i,j-1,k) = wk(i,j,k) - (sign(0.5d0,gA)+sign(0.5d0,gB))*min(abs(gA),abs(gB),abs(gC))
                 gC = f1*( gA+2*gB )
                 wL(i,j  ,k) = wk(i,j,k) + (sign(0.5d0,gA)+sign(0.5d0,gB))*min(abs(gA),abs(gB),abs(gC))
              enddo
           enddo
        enddo
!$omp end do
!$omp do
        do k=ks,ke
           wR(is:ie,jx-1,k) = wk(is:ie,jx,k)
        enddo
!$omp end do
!$omp end parallel

     case(3)
!$omp parallel private(i,j,k,gA,gB,gC,grad)
!$omp workshare
        wL(is:ie,js:je,1) = wk(is:ie,js:je,1)
!$omp end workshare
!$omp do
        do k=2,kx-1
           do j=js,je
              do i=is,ie
                 gA = ( wk(i,j,k  )-wk(i,j,k-1) )
                 gB = ( wk(i,j,k+1)-wk(i,j,k  ) )
                 gC = f1*( 2*gA+gB )
                 wR(i,j,k-1) = wk(i,j,k) - (sign(0.5d0,gA)+sign(0.5d0,gB))*min(abs(gA),abs(gB),abs(gC))
                 gC = f1*( gA+2*gB )
                 wL(i,j,k  ) = wk(i,j,k) + (sign(0.5d0,gA)+sign(0.5d0,gB))*min(abs(gA),abs(gB),abs(gC))
              enddo
           enddo
        enddo
!$omp end do
!$omp workshare
        wR(is:ie,js:je,kx-1) = wk(is:ie,js:je,kx)
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
