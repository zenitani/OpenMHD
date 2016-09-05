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
  integer :: pos(2)

  V = 0.d0

  do j=1,jx
  do i=1,ix

     V(i,j,vx:vz) = U(i,j,mx:mz) / U(i,j,ro)
     rv2 = dot_product( V(i,j,vx:vz), U(i,j,mx:mz) )  ! rho v**2
     B2  = dot_product( U(i,j,bx:bz), U(i,j,bx:bz) )  ! B**2
     V(i,j,pr)    = f1 * ( U(i,j,en) - 0.5d0*(rv2+B2) )

  enddo
  enddo

  if( minval( V(:,:,pr) ) .le. 0 ) then
     pos = minloc(V(:,:,pr))
     write(6,*) 'negative pressure at ',pos,' E: ',U(pos(1),pos(2),en)
     stop
  endif
  if( minval( U(:,:,ro) ) .le. 0 ) then
     pos = minloc(U(:,:,ro))
     write(6,*) 'negative density at ',pos,' ro: ',U(pos(1),pos(2),en)
     stop
  endif

  return
end subroutine u2v
