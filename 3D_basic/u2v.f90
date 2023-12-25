subroutine u2v(U,V,ix,jx,kx)
  implicit none
  include 'param.h'
!-----------------------------------------------------------------------
  integer, intent(in)  :: ix, jx, kx
  real(8), intent(in)  :: U(ix,jx,kx,var1)  ! conserved variables  [input]
  real(8), intent(out) :: V(ix,jx,kx,var2)  ! primitive variables  [output]
!-----------------------------------------------------------------------
  integer :: i, j, k
  real(8) :: B2, rv2
  real(8), parameter :: f1 = gamma - 1
  real(8) :: prmin, rhomin
  integer :: mypos(3)

! V = 0.d0
  prmin  = U(1,1,1,en)
  rhomin = U(1,1,1,ro)

!$omp parallel do private(i,j,k,rv2,B2) reduction(min: prmin,rhomin)
  do k=1,kx
  do j=1,jx
  do i=1,ix

     V(i,j,k,vx:vz) = U(i,j,k,mx:mz) / U(i,j,k,ro)
     rv2 = dot_product( V(i,j,k,vx:vz), U(i,j,k,mx:mz) )  ! rho v**2
     B2  = dot_product( U(i,j,k,bx:bz), U(i,j,k,bx:bz) )  ! B**2
     V(i,j,k,pr)    = f1 * ( U(i,j,k,en) - 0.5d0*(rv2+B2) )
     prmin  = min( V(i,j,k,pr), prmin )
     rhomin = min( U(i,j,k,ro), rhomin )

  enddo
  enddo
  enddo
!$omp end parallel do

  if( prmin <= 0 ) then
     mypos = minloc(V(:,:,:,pr))
     write(6,*) 'negative pressure at ',mypos,' P: ',V(mypos(1),mypos(2),mypos(3),pr)
     stop
  endif
  if( rhomin <= 0 ) then
     mypos = minloc(U(:,:,:,ro))
     write(6,*) 'negative density at ',mypos,' rho: ',U(mypos(1),mypos(2),mypos(3),ro)
     stop
  endif

  return
end subroutine u2v
