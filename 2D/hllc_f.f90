subroutine hllc_f(F,VL,VR,ix,jx)
!-----------------------------------------------------------------------
!     HLLC-G solver in the X direction
!       Ref: K. F. Gurski, SIAM J. Sci. Comput., 25, 2165 (2004)
!-----------------------------------------------------------------------
!     2010/05/11  S. Zenitani  HLLC-G solver
!     2015/08/15  S. Zenitani  opitimization
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
  real(8) :: aL, aR, aM, ay, az
  real(8) :: ptL, ptR, pt, roL, roR, enL, enR, vBL, vBR
  real(8) :: U_hll(var1)

  F(:,:,:) = 0.d0

  call v2u(VL,UL,ix,1,ix-1,jx,1,jx)
  call v2f(VL,FL,ix,1,ix-1,jx,1,jx)
  call v2u(VR,UR,ix,1,ix-1,jx,1,jx)
  call v2f(VR,FR,ix,1,ix-1,jx,1,jx)

!$omp parallel do &
!$omp private(i,j,B2,f1,f2,aL,aR,aM,ay,az) &
!$omp private(ptL,ptR,pt,roL,roR,enL,enR,vBL,vBR,U_hll)
  do j=1,jx
  do i=1,ix-1

!    VL -> coefficients
     B2 = dot_product( VL(i,j,bx:bz), VL(i,j,bx:bz) )
!    f1: Gamma p
!    f2: 4 gamma p B_n^2
     f1 = gamma * VL(i,j,pr)
     f2 = 4 * f1 * VL(i,j,bx)*VL(i,j,bx)
!    fast mode^2
!    vfL = sqrt( ( (f1+B2) + sqrt( (f1+B2)**2 - f2 )) / ( 2*VL(i,j,ro) ))
     aL  = ( (f1+B2) + sqrt(max( (f1+B2)**2-f2, 0.d0 ))) / ( 2*VL(i,j,ro) )
     ptL = VL(i,j,pr) + 0.5d0*B2

!    VR -> coefficients
     B2 = dot_product( VR(i,j,bx:bz), VR(i,j,bx:bz) )
!    f1: Gamma p
!    f2: 4 gamma p B_n^2
     f1 = gamma * VR(i,j,pr)
     f2 = 4 * f1 * VR(i,j,bx)*VR(i,j,bx)
!    fast mode^2
!    vfR = sqrt( ( (f1+B2) + sqrt( (f1+B2)**2 - f2 )) / ( 2*VR(i,j,ro) ))
     aR  = ( (f1+B2) + sqrt(max( (f1+B2)**2-f2, 0.d0 ))) / ( 2*VR(i,j,ro) )
     ptR = VR(i,j,pr) + 0.5d0*B2

!    Riemann fan speed (MK05 eq. 67)
!     aL = min( VL(i,j,vx) - vfL, VR(i,j,vx) - vfR )
!     aR = max( VL(i,j,vx) + vfL, VR(i,j,vx) + vfR )
!     aL = min( VL(i,j,vx), VR(i,j,vx) ) - max( vfL, vfR )
!     aR = max( VL(i,j,vx), VR(i,j,vx) ) + max( vfL, vfR )
     f1 = sqrt( max( aL, aR ) )  ! faster fast-wave (This aL [aR] stores vfL^2 [vfR^2])
     aL = min( min(VL(i,j,vx),VR(i,j,vx))-f1, 0.d0 ) ! *** if (aL > 0), then F = F(L) ***
     aR = max( max(VL(i,j,vx),VR(i,j,vx))+f1, 0.d0 ) ! *** if (aR < 0), then F = F(R) ***

!!    F = F(L)
!     if ( aL .ge. 0 ) then
!        F(i,j,:) = FL(i,j,:)
!!    F = F(R)
!     elseif ( aR .le. 0 ) then
!        F(i,j,:) = FR(i,j,:)
!!    F = F(HLLC)
!     else

     ! HLL
     f1 = 1.d0 / ( aR - aL )
     f2 = aL * aR
     U_hll(mx:mz) = f1*( aR*UR(i,j,mx:mz) - aL*UL(i,j,mx:mz) - FR(i,j,mx:mz) + FL(i,j,mx:mz) )
     U_hll(ro:bz) = f1*( aR*VR(i,j,ro:bz) - aL*VL(i,j,ro:bz) - FR(i,j,ro:bz) + FL(i,j,ro:bz) )
     F(i,j,by)    = f1*( aR*FL(i,j,by) - aL*FR(i,j,by) + f2*(VR(i,j,by)-VL(i,j,by)) )
     F(i,j,bz)    = f1*( aR*FL(i,j,bz) - aL*FR(i,j,bz) + f2*(VR(i,j,bz)-VL(i,j,bz)) )

!    entropy wave
     f2 = 1 / U_hll(ro)
     aM = U_hll(mx) * f2 ! vx_HLL
     ay = U_hll(my) * f2 ! vy_HLL
     az = U_hll(mz) * f2 ! vz_HLL
!    Total pressure (these two are identical)
     pt = ptL + VL(i,j,ro) * ( aL - VL(i,j,vx) ) * ( aM - VL(i,j,vx) )
!    pt = ptR + VR(i,j,ro) * ( aR - VR(i,j,vx) ) * ( aM - VR(i,j,vx) )

!!    F = F(L*) or F(L)
!     if ( aM .ge. 0 ) then

     roL = VL(i,j,ro) * ( aL - VL(i,j,vx) ) / ( aL - aM )
     vBL = dot_product( VL(i,j,vx:vz), VL(i,j,bx:bz) )
     enL = ( ( aL - VL(i,j,vx) )*UL(i,j,en) - ptL*VL(i,j,vx) + pt * aM + &
          U_hll(bx)*( vBL - aM*U_hll(bx) - ay*U_hll(by) - az*U_hll(bz) ) ) / ( aL - aM )
!     F(i,j,mx) = aL *( roL*aM - UL(i,j,mx) ) + FL(i,j,mx)
!     F(i,j,my) = aL *( roL*ay - UL(i,j,my) ) + FL(i,j,my)
!     F(i,j,mz) = aL *( roL*az - UL(i,j,mz) ) + FL(i,j,mz)
!     F(i,j,en) = aL *( enL    - UL(i,j,en) ) + FL(i,j,en)
!     F(i,j,ro) = aL *( roL    - VL(i,j,ro) ) + FL(i,j,ro)

!!    F = F(R*) or F(R)
!     else

     roR = VR(i,j,ro) * ( aR - VR(i,j,vx) ) / ( aR - aM )
     vBR = dot_product( VR(i,j,vx:vz), VR(i,j,bx:bz) )
     enR = ( ( aR - VR(i,j,vx) )*UR(i,j,en) - ptR*VR(i,j,vx) + pt * aM + &
          U_hll(bx)*( vBR - aM*U_hll(bx) - ay*U_hll(by) + az*U_hll(bz) ) ) / ( aR - aM )
!     F(i,j,mx) = aR *( roR*aM - UR(i,j,mx) ) + FR(i,j,mx)
!     F(i,j,my) = aR *( roR*ay - UR(i,j,my) ) + FR(i,j,my)
!     F(i,j,mz) = aR *( roR*az - UR(i,j,mz) ) + FR(i,j,mz)
!     F(i,j,en) = aR *( enR    - UR(i,j,en) ) + FR(i,j,en)
!     F(i,j,ro) = aR *( roR    - VR(i,j,ro) ) + FR(i,j,ro)

!    Weight factor, 0 or 1.  This looks tricky.
!    The code runs 1.0x times slower on Intel, but it runs 1.6 times faster on SPARC.
     f1 = 0.5d0 + sign(0.5d0,aM)  !!  F = F(L*) or F(L)
     f2 = 1.d0 - f1               !!  F = F(R*) or F(R)

     F(i,j,mx) = f1*( aL *( roL*aM - UL(i,j,mx) ) + FL(i,j,mx) ) &
               + f2*( aR *( roR*aM - UR(i,j,mx) ) + FR(i,j,mx) )
     F(i,j,my) = f1*( aL *( roL*ay - UL(i,j,my) ) + FL(i,j,my) ) &
               + f2*( aR *( roR*ay - UR(i,j,my) ) + FR(i,j,my) )
     F(i,j,mz) = f1*( aL *( roL*az - UL(i,j,mz) ) + FL(i,j,mz) ) &
               + f2*( aR *( roR*az - UR(i,j,mz) ) + FR(i,j,mz) )
     F(i,j,en) = f1*( aL *( enL    - UL(i,j,en) ) + FL(i,j,en) ) &
               + f2*( aR *( enR    - UR(i,j,en) ) + FR(i,j,en) )
     F(i,j,ro) = f1*( aL *( roL    - VL(i,j,ro) ) + FL(i,j,ro) ) &
               + f2*( aR *( roR    - VR(i,j,ro) ) + FR(i,j,ro) )

!     endif

!     endif

  enddo
  enddo
!$omp end parallel do

  return
end subroutine hllc_f
