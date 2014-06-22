subroutine v2g(V,G,ix,is,ie,jx,js,je)
!-----------------------------------------------------------------------
!     Convert primitive variables (VL|VR) to numerical flux (G)
!     This V vector contains some conserved variables (rho, B)
!-----------------------------------------------------------------------
  implicit none
  include 'param.h'
!-----------------------------------------------------------------------
  real(8), intent(in) :: V(ix,jx,var1)  ! primitive variables (V) [input]
  real(8) :: G(ix,jx,var1)  ! numerical flux (F) [output]
  integer, intent(in) :: ix, is, ie
  integer, intent(in) :: jx, js, je
!-----------------------------------------------------------------------
  integer :: i, j
  real(8) :: v2, vB, B2, f1

  f1 = gamma / ( gamma - 1 )

  do j=js,je
  do i=is,ie

     v2 = dot_product( V(i,j,vx:vz), V(i,j,vx:vz) )
     vB = dot_product( V(i,j,vx:vz), V(i,j,bx:bz) )
     B2 = dot_product( V(i,j,bx:bz), V(i,j,bx:bz) )
     
     G(i,j,ro) = V(i,j,ro) * V(i,j,vy)
     G(i,j,mx) = V(i,j,ro)*V(i,j,vy)*V(i,j,vx) - V(i,j,by)*V(i,j,bx)
     G(i,j,my) = V(i,j,ro)*V(i,j,vy)*V(i,j,vy) - V(i,j,by)*V(i,j,by) &
          + V(i,j,pr) + 0.5d0*B2
     G(i,j,mz) = V(i,j,ro)*V(i,j,vy)*V(i,j,vz) - V(i,j,by)*V(i,j,bz)
     G(i,j,en) = ( 0.5d0*V(i,j,ro)*v2 + f1*V(i,j,pr) + B2 ) * V(i,j,vy) - vB * V(i,j,by)
     G(i,j,bx) = V(i,j,vy)*V(i,j,bx)-V(i,j,by)*V(i,j,vx)
     G(i,j,by) = 0.d0
     G(i,j,bz) = V(i,j,vy)*V(i,j,bz)-V(i,j,by)*V(i,j,vz)
     G(i,j,ps) = 0.d0
     
  enddo
  enddo

  return
end subroutine v2g
