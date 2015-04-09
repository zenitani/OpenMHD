subroutine v2f(V,F,ix,is,ie,jx,js,je)
!-----------------------------------------------------------------------
!     Convert primitive variables (VL|VR) to numerical flux (F)
!     This V vector contains some conserved variables (rho, B)
!-----------------------------------------------------------------------
  implicit none
  include 'param.h'
!-----------------------------------------------------------------------
  real(8), intent(in) :: V(ix,jx,var1)  ! primitive variables (V) [input]
  real(8) :: F(ix,jx,var1)  ! numerical flux (F) [output]
  integer, intent(in) :: ix, is, ie
  integer, intent(in) :: jx, js, je
!-----------------------------------------------------------------------
  integer :: i, j
  real(8) :: v2, vB, B2
  real(8), parameter :: f1 = gamma / ( gamma - 1 )

  do j=js,je
  do i=is,ie

     v2 = dot_product( V(i,j,vx:vz), V(i,j,vx:vz) )
     vB = dot_product( V(i,j,vx:vz), V(i,j,bx:bz) )
     B2 = dot_product( V(i,j,bx:bz), V(i,j,bx:bz) )
     
     F(i,j,ro) = V(i,j,ro) * V(i,j,vx)
     F(i,j,mx) = V(i,j,ro)*V(i,j,vx)*V(i,j,vx) - V(i,j,bx)*V(i,j,bx) &
          + V(i,j,pr) + 0.5d0*B2
     F(i,j,my) = V(i,j,ro)*V(i,j,vx)*V(i,j,vy) - V(i,j,bx)*V(i,j,by)
     F(i,j,mz) = V(i,j,ro)*V(i,j,vx)*V(i,j,vz) - V(i,j,bx)*V(i,j,bz)
     F(i,j,en) = ( 0.5d0*V(i,j,ro)*v2 + f1*V(i,j,pr) + B2 ) * V(i,j,vx) - vB * V(i,j,bx)
     F(i,j,bx) = 0.d0
     F(i,j,by) = V(i,j,vx)*V(i,j,by)-V(i,j,bx)*V(i,j,vy)
     F(i,j,bz) = V(i,j,vx)*V(i,j,bz)-V(i,j,bx)*V(i,j,vz)
     F(i,j,ps) = 0.d0
     
  enddo
  enddo

  return
end subroutine v2f
