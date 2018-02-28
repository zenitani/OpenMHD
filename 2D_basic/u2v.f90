subroutine u2v(U,V,ix,jx)
  implicit none
  include 'param.h'
!-----------------------------------------------------------------------
  integer, intent(in)  :: ix, jx
  real(8), intent(in)  :: U(ix,jx,var1)  ! conserved variables  [input]
  real(8), intent(out) :: V(ix,jx,var2)  ! primitive variables  [output]
!-----------------------------------------------------------------------
  integer :: i, j
  real(8) :: B2, rv2
  real(8), parameter :: f1 = gamma - 1
  real(8) :: prmin, rhomin
  integer :: mypos(2)

! V = 0.d0

!$omp parallel
!$omp do private(i,j,rv2,B2)
  do j=1,jx
  do i=1,ix

     V(i,j,vx:vz) = U(i,j,mx:mz) / U(i,j,ro)
     rv2 = dot_product( V(i,j,vx:vz), U(i,j,mx:mz) )  ! rho v**2
     B2  = dot_product( U(i,j,bx:bz), U(i,j,bx:bz) )  ! B**2
     V(i,j,pr)    = f1 * ( U(i,j,en) - 0.5d0*(rv2+B2) )

  enddo
  enddo
!$omp end do

!$omp workshare
  prmin = minval( V(:,:,pr) )
!$omp end workshare
!$omp workshare
  rhomin = minval( U(:,:,ro) )
!$omp end workshare
!$omp end parallel

  if( prmin <= 0 ) then
     mypos = minloc(V(:,:,pr))
     write(6,*) 'negative pressure at ',mypos,' P: ',V(mypos(1),mypos(2),pr)
     stop
  endif
  if( rhomin <= 0 ) then
     mypos = minloc(U(:,:,ro))
     write(6,*) 'negative density at ',mypos,' rho: ',U(mypos(1),mypos(2),ro)
     stop
  endif

  return
end subroutine u2v
