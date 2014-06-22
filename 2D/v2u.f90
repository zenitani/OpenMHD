subroutine v2u(V,U,ix,is,ie,jx,js,je)
!-----------------------------------------------------------------------
!     Convert primitive variables (VL|VR) to the conserved quantities (U)
!     This V vector contains some conserved variables (rho, B)
!-----------------------------------------------------------------------
  implicit none
  include 'param.h'
!-----------------------------------------------------------------------
  real(8), intent(in)    :: V(ix,jx,var1)  ! primitive variables (V) [input]
  real(8), intent(inout) :: U(ix,jx,var2)  ! subset of conserved variables (U) [output]
  integer, intent(in) :: ix, is, ie
  integer, intent(in) :: jx, js, je
!-----------------------------------------------------------------------
  integer :: i, j
  real(8) :: v2, B2, f1

  f1 = 1.d0 / ( gamma - 1 )

  do j=js,je
  do i=is,ie
     
     v2 = dot_product( V(i,j,vx:vz), V(i,j,vx:vz) )
     B2 = dot_product( V(i,j,bx:bz), V(i,j,bx:bz) )
     
     U(i,j,mx:mz) = V(i,j,ro) * V(i,j,vx:vz)
     U(i,j,en)    = 0.5d0 * ( V(i,j,ro)*v2 + B2 ) + f1*V(i,j,pr)

  enddo
  enddo

  return
end subroutine v2u
