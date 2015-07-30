subroutine hll_f(F,VL,VR,ix,jx)
!-----------------------------------------------------------------------
!     HLL solver in the X direction
!-----------------------------------------------------------------------
!     2010/09/18  S. Zenitani  HLL solver
!     2015/07/29  S. Zenitani  if-statements ==> max/min functions
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
  real(8) :: aL, aR
  real(8) :: vf, vfL2, vfR2

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
!    fast mode^2
!     vfL = sqrt( ( (f1+B2) + sqrt( (f1+B2)**2 - f2 )) / ( 2*VL(i,j,ro) ))
     vfL2 = ( (f1+B2) + sqrt(max( (f1+B2)**2-f2, 0.d0 ))) / ( 2*VL(i,j,ro) )

!    VR -> coefficients
     B2 = dot_product( VR(i,j,bx:bz), VR(i,j,bx:bz) )
!    f1: Gamma p
!    f2: 4 gamma p B_n^2
     f1 = gamma * VR(i,j,pr)
     f2 = 4 * f1 * VR(i,j,bx)**2
!    fast mode^2
!     vfR = sqrt( ( (f1+B2) + sqrt( (f1+B2)**2 - f2 )) / ( 2*VR(i,j,ro) ))
     vfR2 = ( (f1+B2) + sqrt(max( (f1+B2)**2-f2, 0.d0 ))) / ( 2*VR(i,j,ro) )

!    Riemann fan speed (MK05 eq. 67)
!     aL = min( VL(i,j,vx) - vfL, VR(i,j,vx) - vfR )
!     aR = max( VL(i,j,vx) + vfL, VR(i,j,vx) + vfR )
!     aL = min( VL(i,j,vx), VR(i,j,vx) ) - max( vfL, vfR )
!     aR = max( VL(i,j,vx), VR(i,j,vx) ) + max( vfL, vfR )
     vf = sqrt( max( vfL2, vfR2 ) )
     aL = min( min(VL(i,j,vx),VR(i,j,vx))-vf, 0.d0 )
     aR = max( max(VL(i,j,vx),VR(i,j,vx))+vf, 0.d0 )

!!    F = F(L)
!     if ( aL .ge. 0 ) then
!        F(i,j,:) = FL(i,j,:)
!!    F = F(R)
!     elseif ( aR .le. 0 ) then
!        F(i,j,:) = FR(i,j,:)
!!    F = F(*)
!     else

     f1 = 1.d0 / ( aR - aL )
     f2 = aL * aR
     F(i,j,mx:en) = f1*( aR*FL(i,j,mx:en) - aL*FR(i,j,mx:en) &
          + f2 *(UR(i,j,mx:en)-UL(i,j,mx:en)) )
     F(i,j,ro) = f1*( aR*FL(i,j,ro) - aL*FR(i,j,ro) + f2*(VR(i,j,ro)-VL(i,j,ro)) )
     F(i,j,by) = f1*( aR*FL(i,j,by) - aL*FR(i,j,by) + f2*(VR(i,j,by)-VL(i,j,by)) )
     F(i,j,bz) = f1*( aR*FL(i,j,bz) - aL*FR(i,j,bz) + f2*(VR(i,j,bz)-VL(i,j,bz)) )

!     endif

  enddo
  enddo

  return
end subroutine hll_f
