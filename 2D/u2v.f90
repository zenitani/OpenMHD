subroutine u2v(U,V,ix,jx)
  implicit none
  include 'param.h'
!-----------------------------------------------------------------------
  integer, intent(in)  :: ix, jx
  real(8), intent(in)  :: U(ix,jx,var1)  ! conserved variables  [input]
  real(8), intent(out) :: V(ix,jx,var2)  ! primitive variables  [output]
!-----------------------------------------------------------------------
  integer :: i, j
  real(8) :: B2, mv2
  real(8), parameter :: f1 = gamma - 1
  integer :: pos(2)

  V(:,:,:) = 0.d0
  do j=1,jx
  do i=1,ix

     V(i,j,vx:vz) = U(i,j,mx:mz) / U(i,j,ro)
     mv2 = dot_product( V(i,j,vx:vz), U(i,j,mx:mz) )  ! m v**2
     B2  = dot_product( U(i,j,bx:bz), U(i,j,bx:bz) )  ! B**2
     V(i,j,pr)    = f1 * ( U(i,j,en) - 0.5d0*(mv2+B2) )

  enddo
  enddo

  if( minval( V(:,:,pr) ) .le. 0 ) then
     pos = minloc(V(:,:,pr))
     i = pos(1)
     j = pos(2)
     write(6,*) 'negative pressure at ', i, j, ' E: ', U(i,j,en)
     stop
  endif
  if( minval( U(:,:,ro) ) .le. 0 ) then
     pos = minloc(U(:,:,ro))
     i = pos(1)
     j = pos(2)
     write(6,*) 'negative density at ', i, j, ' ro: ', U(i,j,ro)
     stop
  endif

  return
end subroutine u2v
