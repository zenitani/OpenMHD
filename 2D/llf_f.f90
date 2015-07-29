subroutine llf_f(F,VL,VR,ix,jx)
!-----------------------------------------------------------------------
!     Local Lax-Friedrich (Rusanov) solver in the X direction
!-----------------------------------------------------------------------
!     2010/05/13  S. Zenitani  LLF solver
!-----------------------------------------------------------------------
  implicit none
  include 'param.h'
  integer, intent(in)  :: ix, jx
! numerical flux (F) [output]
  real(8), intent(out) :: F(ix,jx,var1)
! left/right states (VL,VR) [input]
  real(8), intent(in)  :: VL(ix,jx,var1), VR(ix,jx,var1)
! left/right conserved variables (UL & UR; local)
  real(8) :: UL(ix,jx,var2), UR(ix,jx,var2)
! numerical flux (FL & FR local)
  real(8) :: FL(ix,jx,var1), FR(ix,jx,var1)
  integer :: i, j

  real(8) :: B2, f1, f2
  real(8) :: aLLF
  real(8) :: vfL, vfR

  F(:,:,:) = 0.d0

  call v2u(VL,UL,ix,1,ix-1,jx,1,jx)
  call v2f(VL,FL,ix,1,ix-1,jx,1,jx)
  call v2u(VR,UR,ix,1,ix-1,jx,1,jx)
  call v2f(VR,FR,ix,1,ix-1,jx,1,jx)

  do j=1,jx
  do i=1,ix-1

!    VL -> coefficients
     B2 = dot_product( VL(i,j,bx:bz), VL(i,j,bx:bz) )
!    f1: Gamma p
!    f2: 4 gamma p B_n^2
     f1 = gamma * VL(i,j,pr)
     f2 = 4 * f1 * VL(i,j,bx)**2
!    fast mode
!     vfL = sqrt( ( (f1+B2) + sqrt( (f1+B2)**2 - f2 )) / ( 2*VL(i,j,ro) ))
     vfL = sqrt( ( (f1+B2) + sqrt(max( (f1+B2)**2-f2, 0.d0 ))) / ( 2*VL(i,j,ro) ))

!    VR -> coefficients
     B2 = dot_product( VR(i,j,bx:bz), VR(i,j,bx:bz) )
!    f1: Gamma p
!    f2: 4 gamma p B_n^2
     f1 = gamma * VR(i,j,pr)
     f2 = 4 * f1 * VR(i,j,bx)**2
!    fast mode
!     vfR = sqrt( ( (f1+B2) + sqrt( (f1+B2)**2 - f2 )) / ( 2*VR(i,j,ro) ))
     vfR = sqrt( ( (f1+B2) + sqrt(max( (f1+B2)**2-f2, 0.d0 ))) / ( 2*VR(i,j,ro) ))

!    Local Lax Friedrich
     aLLF = max( abs(VL(i,j,vx))+vfL, abs(VR(i,j,vx))+vfR )

     F(i,j,mx:en) = 0.5d0 * ( FL(i,j,mx:en)+FR(i,j,mx:en) &
          - aLLF*( UR(i,j,mx:en)-UL(i,j,mx:en) ) )
     F(i,j,ro:bz) = 0.5d0 * ( FL(i,j,ro:bz)+FR(i,j,ro:bz) &
          - aLLF*( VR(i,j,ro:bz)-VL(i,j,ro:bz) ) )

  enddo
  enddo

  return
end subroutine llf_f
