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
  integer :: pos1(2), pos2(2)

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

!$omp sections
!$omp section
  if( minval( V(:,:,pr) ) <= 0 ) then
     pos1 = minloc(V(:,:,pr))
     write(6,*) 'negative pressure at ',pos1,' P: ',V(pos1(1),pos1(2),pr)
     stop
  endif
!$omp section
  if( minval( U(:,:,ro) ) <= 0 ) then
     pos2 = minloc(U(:,:,ro))
     write(6,*) 'negative density at ',pos2,' rho: ',U(pos2(1),pos2(2),ro)
     stop
  endif
!$omp end sections
!$omp end parallel

  return
end subroutine u2v
