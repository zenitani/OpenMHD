subroutine set_eta(E,EF,EG,x,y,dx,Rm1,Rm0,ix,jx)
!-----------------------------------------------------------------------
!     spatial profile of resistivity
!-----------------------------------------------------------------------
  implicit none
  real(8), intent(out) :: E(ix,jx), EF(ix,jx), EG(ix,jx)
  real(8), intent(in)  :: x(ix), y(jx), dx
  real(8), intent(in)  :: Rm1, Rm0
  integer, intent(in)  :: ix, jx
!-----------------------------------------------------------------------
  integer :: i, j
  real(8) :: dx2, eta0, eta01
  real(8) :: eta
!-----------------------------------------------------------------------

  dx2 = dx/2.d0
  eta0  = 1.d0 / Rm0
  eta01 = ( Rm0 - Rm1 ) / ( Rm0 * Rm1 )

  E(:,:) = 0.d0
  do j=1,jx
  do i=1,ix
     E(i,j) = eta( x(i), y(j), eta0, eta01 )
  enddo
  enddo

  EF(ix,:) = 0.d0
  do j=1,jx
  do i=1,ix-1
     EF(i,j) = eta( x(i)+dx2, y(j), eta0, eta01 )
  enddo
  enddo

  EG(:,jx) = 0.d0
  do j=1,jx-1
  do i=1,ix
     EG(i,j) = eta( x(i), y(j)+dx2, eta0, eta01 )
  enddo
  enddo

end subroutine set_eta

function eta(xx,yy,eta0,eta01)
  real(8) :: eta
  real(8), intent(in) :: xx, yy, eta0, eta01

!  eta = eta0 + eta01 * ( cosh( sqrt( xx**2+yy**2 )/2.d0 ) )**(-2)
  eta = eta0 + eta01 * ( cosh( sqrt( xx**2+yy**2 ) ) )**(-2)

end function eta
