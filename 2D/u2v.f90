subroutine u2v(U,V,ix,jx)
  implicit none
  include 'param.h'
!-----------------------------------------------------------------------
  integer, intent(in)  :: ix, jx
  real(8), intent(in)  :: U(ix,jx,var1)  ! conserved variables  [input]
  real(8), intent(out) :: V(ix,jx,var2)  ! primitive variables  [output]
!-----------------------------------------------------------------------
  integer :: i, j
  real(8) :: B2, M2
  real(8), parameter :: f1 = gamma - 1

  do j=1,jx
  do i=1,ix

     B2 = dot_product( U(i,j,bx:bz), U(i,j,bx:bz) )  ! B**2
     M2 = dot_product( U(i,j,mx:mz), U(i,j,mx:mz) )  ! M**2
     V(i,j,vx:vz) = U(i,j,mx:mz) / U(i,j,ro)
     V(i,j,pr)    = f1 * ( U(i,j,en) - 0.5d0*(M2/U(i,j,ro)+B2) )

     if( U(i,j,ro) .le. 0 ) then
        write(6,*) 'negative density at ', i, j
        write(6,*) 'density: ',U(i,j,ro)
        stop
     endif
     if( V(i,j,pr) .le. 0 ) then
        write(6,*) 'negative pressure at ', i, j
        write(6,*) ' E: ', U(i,j,en), '  E_B: ', 0.5d0*B2, '  E_K: ', 0.5d0*(M2/U(i,j,ro))
        stop
     endif

  enddo
  enddo

  return
end subroutine u2v
