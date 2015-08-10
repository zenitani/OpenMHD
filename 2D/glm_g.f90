subroutine glm_g(G,VL,VR,ch,ix,jx)
!-----------------------------------------------------------------------
!     GLM flux terms in the Y direction
!-----------------------------------------------------------------------
!     This needs to be called right after LLF/HLL/HLLC/HLLD flux solver,
!     because these solvers may set G(:,:,by) and G(:,:,ps) to zero
!-----------------------------------------------------------------------
!     2010/09/25  S. Zenitani  GLM solver
!-----------------------------------------------------------------------
  implicit none
  include 'param.h'
  integer, intent(in) :: ix, jx
! numerical flux (G) [input/output]
  real(8), intent(inout) :: G(ix,jx,var1)
! left/right states (VL,VR) [input]
  real(8), intent(in) :: VL(ix,jx,var1), VR(ix,jx,var1)
! div cleaning wave speed
  real(8), intent(in) :: ch
  real(8) :: ch2
  integer :: i, j

  ch2 = ch**2

  do j=1,jx-1
  do i=1,ix

     G(i,j,by) = 0.5d0*( VL(i,j,ps) + VR(i,j,ps) &
          - ch*( VR(i,j,by)-VL(i,j,by) ) )
     G(i,j,ps) = 0.5d0*( ch2*( VL(i,j,by)+VR(i,j,by) ) &
          - ch*( VR(i,j,ps)-VL(i,j,ps) ) )

  enddo
  enddo

  return
end subroutine glm_g
