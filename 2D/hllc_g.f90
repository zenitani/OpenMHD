subroutine hllc_g(G,VL,VR,ix,jx)
!-----------------------------------------------------------------------
!     HLLC-G solver in the Y direction
!       Ref: K. F. Gurski, SIAM J. Sci. Comput., 25, 2165 (2004)
!-----------------------------------------------------------------------
!     2010/05/11  S. Zenitani  HLLC-G solver
!     2015/08/15  S. Zenitani  opitimization
!-----------------------------------------------------------------------
  implicit none
  include 'param.h'
  integer, intent(in)  :: ix, jx
! numerical flux (F) [output]
  real(8), intent(out) :: G(ix,jx,var1)
! left/right states (VL,VR) [input]
  real(8), intent(in)  :: VL(ix,jx,var1), VR(ix,jx,var1)
! left/right conserved variables (UL & UR; local)
  real(8) :: UL(ix,jx,var2), UR(ix,jx,var2)
! numerical flux (GL & GR local)
  real(8) :: GL(ix,jx,var1), GR(ix,jx,var1)
  integer :: i, j

  real(8) :: B2, f1, f2
  real(8) :: aL, aR, aM, ax, az
  real(8) :: vfL2, vfR2
  real(8) :: ptL, ptR, pt, roL, roR, enL, enR, vBL, vBR
  real(8) :: U_hll(var1)

  G(:,:,:) = 0.d0

  call v2u(VL,UL,ix,1,ix,jx,1,jx-1)
  call v2g(VL,GL,ix,1,ix,jx,1,jx-1)
  call v2u(VR,UR,ix,1,ix,jx,1,jx-1)
  call v2g(VR,GR,ix,1,ix,jx,1,jx-1)

!$omp parallel do &
!$omp private(i,B2,f1,f2,aL,aR,aM,ax,az,vfL2,vfR2) &
!$omp private(ptL,ptR,pt,roL,roR,enL,enR,vBL,vBR,U_hll)
  do j=1,jx-1
  do i=1,ix

!    VL -> coefficients
     B2 = dot_product( VL(i,j,bx:bz), VL(i,j,bx:bz) )
!    f1: Gamma p
!    f2: 4 gamma p B_n^2
     f1 = gamma * VL(i,j,pr)
     f2 = 4 * f1 * VL(i,j,by)**2
!    fast mode^2
!     vfL = sqrt( ( (f1+B2) + sqrt( (f1+B2)**2 - f2 )) / ( 2*VL(i,j,ro) ))
     vfL2 = ( (f1+B2) + sqrt(max( (f1+B2)**2-f2, 0.d0 ))) / ( 2*VL(i,j,ro) )
     ptL = VL(i,j,pr) + 0.5d0*B2

!    VR -> coefficients
     B2 = dot_product( VR(i,j,bx:bz), VR(i,j,bx:bz) )
!    f1: Gamma p
!    f2: 4 gamma p B_n^2
     f1 = gamma * VR(i,j,pr)
     f2 = 4 * f1 * VR(i,j,by)**2
!    fast mode^2
!     vfR = sqrt( ( (f1+B2) + sqrt( (f1+B2)**2 - f2 )) / ( 2*VR(i,j,ro) ))
     vfR2 = ( (f1+B2) + sqrt(max( (f1+B2)**2-f2, 0.d0 ))) / ( 2*VR(i,j,ro) )
     ptR = VR(i,j,pr) + 0.5d0*B2

!    Riemann fan speed (MK05 eq. 67)
!     aL = min( VL(i,j,vy) - vfL, VR(i,j,vy) - vfR )
!     aR = max( VL(i,j,vy) + vfL, VR(i,j,vy) + vfR )
!     aL = min( VL(i,j,vy), VR(i,j,vy) ) - max( vfL, vfR )
!     aR = max( VL(i,j,vy), VR(i,j,vy) ) + max( vfL, vfR )
     f1 = sqrt( max( vfL2, vfR2 ) )
     aL = min( min(VL(i,j,vy),VR(i,j,vy))-f1, 0.d0 ) ! *** if (aL > 0), then G = G(L) ***
     aR = max( max(VL(i,j,vy),VR(i,j,vy))+f1, 0.d0 ) ! *** if (aR < 0), then G = G(R) ***

!!    G = G(L)
!     if ( aL .gt. 0 ) then
!        G(i,j,:) = GL(i,j,:)
!!    G = G(R)
!     elseif ( aR .lt. 0 ) then
!        G(i,j,:) = GR(i,j,:)
!!    G = G(HLLC)
!     else

     ! HLL
     f1 = 1.d0 / ( aR - aL )
     f2 = aL * aR
     U_hll(mx:mz) = f1*( aR*UR(i,j,mx:mz) - aL*UL(i,j,mx:mz) - GR(i,j,mx:mz) + GL(i,j,mx:mz) )
     U_hll(ro:bz) = f1*( aR*VR(i,j,ro:bz) - aL*VL(i,j,ro:bz) - GR(i,j,ro:bz) + GL(i,j,ro:bz) )
     G(i,j,bx)    = f1*( aR*GL(i,j,bx) - aL*GR(i,j,bx) + f2*(VR(i,j,bx)-VL(i,j,bx)) )
     G(i,j,bz)    = f1*( aR*GL(i,j,bz) - aL*GR(i,j,bz) + f2*(VR(i,j,bz)-VL(i,j,bz)) )

!    entropy wave
     f2 = 1 / U_hll(ro)
     ax = U_hll(mx) * f2 ! vx_HLL
     aM = U_hll(my) * f2 ! vy_HLL
     az = U_hll(mz) * f2 ! vz_HLL
!    Total pressure (these two are identical)
     pt = ptL + VL(i,j,ro) * ( aL - VL(i,j,vy) ) * ( aM - VL(i,j,vy) )
!    pt = ptR + VR(i,j,ro) * ( aR - VR(i,j,vy) ) * ( aM - VR(i,j,vy) )

!!    G = G(L*) or G(L)
!     if ( aM .ge. 0 ) then

     roL = VL(i,j,ro) * ( aL - VL(i,j,vy) ) / ( aL - aM )
     vBL = dot_product( VL(i,j,vx:vz), VL(i,j,bx:bz) )
     enL = ( ( aL - VL(i,j,vy) )*UL(i,j,en) - ptL*VL(i,j,vy) + pt * aM + &
             U_hll(by)*( vBL - ax*U_hll(bx) - aM*U_hll(by) - az*U_hll(bz) ) ) /  ( aL - aM )
!     G(i,j,mx) = aL *( ro_tmp*ax - UL(i,j,mx) ) + GL(i,j,mx)
!     G(i,j,my) = aL *( ro_tmp*aM - UL(i,j,my) ) + GL(i,j,my)
!     G(i,j,mz) = aL *( ro_tmp*az - UL(i,j,mz) ) + GL(i,j,mz)
!     G(i,j,en) = aL *( en_tmp    - UL(i,j,en) ) + GL(i,j,en)
!     G(i,j,ro) = aL *( ro_tmp    - VL(i,j,ro) ) + GL(i,j,ro)

!!    G = G(R*) or G(R)
!     else

     roR = VR(i,j,ro) * ( aR - VR(i,j,vy) ) / ( aR - aM )
     vBR = dot_product( VR(i,j,vx:vz), VR(i,j,bx:bz) )
     enR = ( ( aR - VR(i,j,vy) )*UR(i,j,en) - ptR*VR(i,j,vy) + pt * aM + &
          U_hll(by)*( vBR - ax*U_hll(bx) - aM*U_hll(by) - az*U_hll(bz) ) ) / ( aR - aM )
!     G(i,j,mx) = aR *( ro_tmp*ax - UR(i,j,mx) ) + GR(i,j,mx)
!     G(i,j,my) = aR *( ro_tmp*aM - UR(i,j,my) ) + GR(i,j,my)
!     G(i,j,mz) = aR *( ro_tmp*az - UR(i,j,mz) ) + GR(i,j,mz)
!     G(i,j,en) = aR *( en_tmp    - UR(i,j,en) ) + GR(i,j,en)
!     G(i,j,ro) = aR *( ro_tmp    - VR(i,j,ro) ) + GR(i,j,ro)

!    Weight factor, 0 or 1.  This looks tricky.
!    The code runs 1.0x times slower on Intel, but it runs 1.6 times faster on SPARC.
     f1 = 0.5d0 + sign(0.5d0,aM)  !!  G = G(L*) or G(L)
     f2 = 1.d0 - f1               !!  G = G(R*) or G(R)

     G(i,j,mx) = f1 * ( aL *( roL*ax - UL(i,j,mx) ) + GL(i,j,mx) ) &
               + f2 * ( aR *( roR*ax - UR(i,j,mx) ) + GR(i,j,mx) )
     G(i,j,my) = f1 * ( aL *( roL*aM - UL(i,j,my) ) + GL(i,j,my) ) &
               + f2 * ( aR *( roR*aM - UR(i,j,my) ) + GR(i,j,my) )
     G(i,j,mz) = f1 * ( aL *( roL*az - UL(i,j,mz) ) + GL(i,j,mz) ) &
               + f2 * ( aR *( roR*az - UR(i,j,mz) ) + GR(i,j,mz) )
     G(i,j,en) = f1 * ( aL *( enL    - UL(i,j,en) ) + GL(i,j,en) ) &
               + f2 * ( aR *( enR    - UR(i,j,en) ) + GR(i,j,en) )
     G(i,j,ro) = f1 * ( aL *( roL    - VL(i,j,ro) ) + GL(i,j,ro) ) &
               + f2 * ( aR *( roR    - VR(i,j,ro) ) + GR(i,j,ro) )

!     endif

!     endif

  enddo
  enddo
!$omp end parallel do

  return
end subroutine hllc_g
