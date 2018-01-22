subroutine flux_solver(F,VL,VR,ix,jx,dir,type)
!-----------------------------------------------------------------------
!     HLLC-G flux solver
!       Ref: K. F. Gurski, SIAM J. Sci. Comput., 25, 2165 (2004)
!     HLLD flux solver
!       Ref: T. Miyoshi, K. Kusano, J. Comput. Phys., 208, 315 (2005)
!-----------------------------------------------------------------------
!     2010/05/11  S. Zenitani  HLLC-G solver
!     2010/05/12  S. Zenitani  HLLD solver
!     2010/09/18  S. Zenitani  HLLD solver: fixed some bugs
!     2015/07/29  S. Zenitani  HLL solver: if-statements ==> max/min functions
!     2015/08/15  S. Zenitani  HLLC-G solver: optimization
!     2016/09/06  S. Zenitani  X/Y/Z directions
!-----------------------------------------------------------------------
  implicit none
  include 'param.h'
!-----------------------------------------------------------------------
! size of arrays [input]
  integer, intent(in)  :: ix, jx
! numerical flux (F) [output]
  real(8), intent(out) :: F(ix,jx,var1)
! left/right states (VL,VR) [input]
  real(8), intent(in)  :: VL(ix,jx,var1), VR(ix,jx,var1)
! direction [input]: 1 (X), 2 (Y), 3 (Z)
  integer, intent(in)  :: dir
! numerical flux [input]: 0 (LLF Rusanov), 1 (HLL), 2 (HLLC), 3 (HLLD)
  integer, intent(in)  :: type
!-----------------------------------------------------------------------
! left/right conserved variables (UL & UR; local)
  real(8) :: UL(var2), UR(var2)
! left/right numerical flux (FL & FR; local)
  real(8) :: FL(var1), FR(var1)
  integer :: i, j, is=0, ie=0, js=0, je=0
! directions
  integer :: vn, vt1, vt2, mn, mt1, mt2, bn, bt1, bt2
! constants
  real(8), parameter :: r1 =  1.d0 / ( gamma - 1.d0 )  ! v2u
  real(8), parameter :: r2 = gamma / ( gamma - 1.d0 )  ! v2f

! HLL solver
  real(8) :: f1, f2, B2, v2, vBL, vBR
  real(8) :: aL, aR

! HLLC solver
  real(8) :: aM, at1, at2
  real(8) :: ptL, ptR, pt, roL, roR, enL, enR
  real(8) :: U_hll(var1)

! HLLD solver
  real(8) :: aL1, aR1, aVL, aVR
  real(8) :: UL1(var1), UR1(var1), U2(var1)
  real(8) :: vt1L, vt1R, vt2L, vt2R, roLs, roRs
  real(8), parameter :: hllg_factor = 1.001d0   ! threshold to switch HLLC-G and HLLD
!-----------------------------------------------------------------------

  F = 0.d0
  FL = 0.d0; FR = 0.d0

! directions
  vn  = vx + mod(dir-1,3)
  vt1 = vx + mod(dir  ,3)
  vt2 = vx + mod(dir+1,3)
  bn  = bx + mod(dir-1,3)
  bt1 = bx + mod(dir  ,3)
  bt2 = bx + mod(dir+1,3)
  mn  = mx + mod(dir-1,3)
  mt1 = mx + mod(dir  ,3)
  mt2 = mx + mod(dir+1,3)

! boundaries
  select case(dir)
  case(1)
     is = 1; ie = ix-1
     js = min(2,jx); je = max(1,jx-1)
  case(2)
     is = min(2,ix); ie = max(1,ix-1)
     js = 1; je = jx-1
  case(3)
  endselect

! main switch
  select case(type)
  !-----------------------------------------------------------------------
  !  Local Lax Friedrich (Rusanov)
  !  No OpenMP optimization - just for education
  !-----------------------------------------------------------------------
  case(0)

     do j=js,je
     do i=is,ie

!    VL --> UL
     v2 = dot_product( VL(i,j,vx:vz), VL(i,j,vx:vz) )
     B2 = dot_product( VL(i,j,bx:bz), VL(i,j,bx:bz) )
     vBL= dot_product( VL(i,j,vx:vz), VL(i,j,bx:bz) )
     UL(mx:mz) = VL(i,j,ro) * VL(i,j,vx:vz)
     UL(en)    = 0.5d0 * ( VL(i,j,ro)*v2 + B2 ) + r1*VL(i,j,pr)
!    VL --> FL
     FL(ro)  = VL(i,j,ro) * VL(i,j,vn)
     FL(mn)  = VL(i,j,ro)*VL(i,j,vn)*VL(i,j,vn)  - VL(i,j,bn)*VL(i,j,bn) + VL(i,j,pr) + 0.5d0*B2
     FL(mt1) = VL(i,j,ro)*VL(i,j,vn)*VL(i,j,vt1) - VL(i,j,bn)*VL(i,j,bt1)
     FL(mt2) = VL(i,j,ro)*VL(i,j,vn)*VL(i,j,vt2) - VL(i,j,bn)*VL(i,j,bt2)
     FL(en)  = ( 0.5d0*VL(i,j,ro)*v2 + r2*VL(i,j,pr) + B2 ) * VL(i,j,vn) - vBL * VL(i,j,bn)
!    FL(bn)  = 0.d0
     FL(bt1) = VL(i,j,vn)*VL(i,j,bt1)-VL(i,j,bn)*VL(i,j,vt1)
     FL(bt2) = VL(i,j,vn)*VL(i,j,bt2)-VL(i,j,bn)*VL(i,j,vt2)
!    FL(ps)  = 0.d0

!    VL --> fast mode
     f1 = gamma * VL(i,j,pr)      ! f1: gamma p
     f2 = 4 * f1 * VL(i,j,bn)**2  ! f2: 4 gamma p B_n^2
     aL = sqrt( ( (f1+B2) + sqrt( (f1+B2)**2 - f2 )) / ( 2*VL(i,j,ro) ))

!    VR --> UR
     v2 = dot_product( VR(i,j,vx:vz), VR(i,j,vx:vz) )
     B2 = dot_product( VR(i,j,bx:bz), VR(i,j,bx:bz) )
     vBR= dot_product( VR(i,j,vx:vz), VR(i,j,bx:bz) )
     UR(mx:mz) = VR(i,j,ro) * VR(i,j,vx:vz)
     UR(en)    = 0.5d0 * ( VR(i,j,ro)*v2 + B2 ) + r1*VR(i,j,pr)
!    VR --> FR
     FR(ro)  = VR(i,j,ro) * VR(i,j,vn)
     FR(mn)  = VR(i,j,ro)*VR(i,j,vn)*VR(i,j,vn)  - VR(i,j,bn)*VR(i,j,bn) + VR(i,j,pr) + 0.5d0*B2
     FR(mt1) = VR(i,j,ro)*VR(i,j,vn)*VR(i,j,vt1) - VR(i,j,bn)*VR(i,j,bt1)
     FR(mt2) = VR(i,j,ro)*VR(i,j,vn)*VR(i,j,vt2) - VR(i,j,bn)*VR(i,j,bt2)
     FR(en)  = ( 0.5d0*VR(i,j,ro)*v2 + r2*VR(i,j,pr) + B2 ) * VR(i,j,vn) - vBR * VR(i,j,bn)
!    FR(bn)  = 0.d0
     FR(bt1) = VR(i,j,vn)*VR(i,j,bt1)-VR(i,j,bn)*VR(i,j,vt1)
     FR(bt2) = VR(i,j,vn)*VR(i,j,bt2)-VR(i,j,bn)*VR(i,j,vt2)
!    FR(ps)  = 0.d0

!    VR --> fast mode
     f1 = gamma * VR(i,j,pr)      ! f1: gamma p
     f2 = 4 * f1 * VR(i,j,bn)**2  ! f2: 4 gamma p B_n^2
     aR = sqrt( ( (f1+B2) + sqrt( (f1+B2)**2 - f2 )) / ( 2*VR(i,j,ro) ))

     f1 = max( abs(VL(i,j,vn))+aL, abs(VR(i,j,vn))+aR )
     F(i,j,mx:en) = 0.5d0 * ( FL(mx:en)+FR(mx:en) - f1*( UR(mx:en)-UL(mx:en) ) )
     F(i,j,ro)    = 0.5d0 * ( FL(ro)   +FR(ro)    - f1*( VR(i,j,ro) -VL(i,j,ro)  ) )
     F(i,j,bt1)   = 0.5d0 * ( FL(bt1)  +FR(ro)    - f1*( VR(i,j,bt1)-VL(i,j,bt1) ) )
     F(i,j,bt2)   = 0.5d0 * ( FL(bt2)  +FR(ro)    - f1*( VR(i,j,bt2)-VL(i,j,bt2) ) )

     enddo
     enddo

  !-----------------------------------------------------------------------
  !  HLL solver
  !-----------------------------------------------------------------------
  case(1)
!$omp parallel do private(i,j,v2,B2,vBL,vBR,UL,UR,FL,FR,f1,f2,aL,aR)
     do j=js,je
     do i=is,ie

!    VL --> UL
     v2 = dot_product( VL(i,j,vx:vz), VL(i,j,vx:vz) )
     B2 = dot_product( VL(i,j,bx:bz), VL(i,j,bx:bz) )
     vBL= dot_product( VL(i,j,vx:vz), VL(i,j,bx:bz) )
     UL(mx:mz) = VL(i,j,ro) * VL(i,j,vx:vz)
     UL(en)    = 0.5d0 * ( VL(i,j,ro)*v2 + B2 ) + r1*VL(i,j,pr)
!    VL --> FL
     FL(ro)  = VL(i,j,ro) * VL(i,j,vn)
     FL(mn)  = VL(i,j,ro)*VL(i,j,vn)*VL(i,j,vn)  - VL(i,j,bn)*VL(i,j,bn) + VL(i,j,pr) + 0.5d0*B2
     FL(mt1) = VL(i,j,ro)*VL(i,j,vn)*VL(i,j,vt1) - VL(i,j,bn)*VL(i,j,bt1)
     FL(mt2) = VL(i,j,ro)*VL(i,j,vn)*VL(i,j,vt2) - VL(i,j,bn)*VL(i,j,bt2)
     FL(en)  = ( 0.5d0*VL(i,j,ro)*v2 + r2*VL(i,j,pr) + B2 ) * VL(i,j,vn) - vBL * VL(i,j,bn)
     FL(bt1) = VL(i,j,vn)*VL(i,j,bt1)-VL(i,j,bn)*VL(i,j,vt1)
     FL(bt2) = VL(i,j,vn)*VL(i,j,bt2)-VL(i,j,bn)*VL(i,j,vt2)

!    VL --> fast mode^2
     f1 = gamma * VL(i,j,pr)      ! f1: gamma p
     f2 = 4 * f1 * VL(i,j,bn)**2  ! f2: 4 gamma p B_n^2
!    vfL= sqrt( ( (f1+B2) + sqrt( (f1+B2)**2 - f2 )) / ( 2*VL(i,j,ro) ))
     aL = ( (f1+B2) + sqrt(max( (f1+B2)**2-f2, 0.d0 ))) / ( 2*VL(i,j,ro) )

!    VR --> UR
     v2 = dot_product( VR(i,j,vx:vz), VR(i,j,vx:vz) )
     B2 = dot_product( VR(i,j,bx:bz), VR(i,j,bx:bz) )
     vBR= dot_product( VR(i,j,vx:vz), VR(i,j,bx:bz) )
     UR(mx:mz) = VR(i,j,ro) * VR(i,j,vx:vz)
     UR(en)    = 0.5d0 * ( VR(i,j,ro)*v2 + B2 ) + r1*VR(i,j,pr)
!    VR --> FR
     FR(ro)  = VR(i,j,ro) * VR(i,j,vn)
     FR(mn)  = VR(i,j,ro)*VR(i,j,vn)*VR(i,j,vn)  - VR(i,j,bn)*VR(i,j,bn) + VR(i,j,pr) + 0.5d0*B2
     FR(mt1) = VR(i,j,ro)*VR(i,j,vn)*VR(i,j,vt1) - VR(i,j,bn)*VR(i,j,bt1)
     FR(mt2) = VR(i,j,ro)*VR(i,j,vn)*VR(i,j,vt2) - VR(i,j,bn)*VR(i,j,bt2)
     FR(en)  = ( 0.5d0*VR(i,j,ro)*v2 + r2*VR(i,j,pr) + B2 ) * VR(i,j,vn) - vBR * VR(i,j,bn)
     FR(bt1) = VR(i,j,vn)*VR(i,j,bt1)-VR(i,j,bn)*VR(i,j,vt1)
     FR(bt2) = VR(i,j,vn)*VR(i,j,bt2)-VR(i,j,bn)*VR(i,j,vt2)

!    VR --> fast mode^2
     f1 = gamma * VR(i,j,pr)      ! f1: gamma p
     f2 = 4 * f1 * VR(i,j,bn)**2  ! f2: 4 gamma p B_n^2
!    vfR= sqrt( ( (f1+B2) + sqrt( (f1+B2)**2 - f2 )) / ( 2*VR(i,j,ro) ))
     aR = ( (f1+B2) + sqrt(max( (f1+B2)**2-f2, 0.d0 ))) / ( 2*VR(i,j,ro) )

!    Riemann fan speed (MK05 eq. 67)
!     aL = min( VL(i,j,vn) - vfL, VR(i,j,vn) - vfR )
!     aR = max( VL(i,j,vn) + vfL, VR(i,j,vn) + vfR )
!     aL = min( VL(i,j,vn), VR(i,j,vn) ) - max( vfL, vfR )
!     aR = max( VL(i,j,vn), VR(i,j,vn) ) + max( vfL, vfR )
     f1 = sqrt( max( aL, aR ) )  ! faster fast-wave (This aL [aR] stores vfL^2 [vfR^2])
     aL = min( min(VL(i,j,vn),VR(i,j,vn))-f1, 0.d0 ) ! *** if (aL > 0), then F = F(L) ***
     aR = max( max(VL(i,j,vn),VR(i,j,vn))+f1, 0.d0 ) ! *** if (aR < 0), then F = F(R) ***

     f1 = 1.d0 / ( aR - aL )
     f2 = aL * aR
     F(i,j,mx:en) = f1*( aR*FL(mx:en) - aL*FR(mx:en) + f2 *(UR(mx:en)-UL(mx:en)) )
     F(i,j,ro)    = f1*( aR*FL(ro)  - aL*FR(ro)  + f2*(VR(i,j,ro) -VL(i,j,ro) ) )
     F(i,j,bt1)   = f1*( aR*FL(bt1) - aL*FR(bt1) + f2*(VR(i,j,bt1)-VL(i,j,bt1)) )
     F(i,j,bt2)   = f1*( aR*FL(bt2) - aL*FR(bt2) + f2*(VR(i,j,bt2)-VL(i,j,bt2)) )

     enddo
     enddo
!$omp end parallel do

  !-----------------------------------------------------------------------
  !  HLLC solver
  !-----------------------------------------------------------------------
  case(2)
!$omp parallel do &
!$omp private(i,j,v2,B2,vBL,vBR,UL,UR,FL,FR,f1,f2,aL,aR,aM,at1,at2) &
!$omp private(ptL,ptR,pt,roL,roR,enL,enR,U_hll)
     do j=js,je
     do i=is,ie

!    VL --> UL
     v2 = dot_product( VL(i,j,vx:vz), VL(i,j,vx:vz) )
     B2 = dot_product( VL(i,j,bx:bz), VL(i,j,bx:bz) )
     vBL= dot_product( VL(i,j,vx:vz), VL(i,j,bx:bz) )
     UL(mx:mz) = VL(i,j,ro) * VL(i,j,vx:vz)
     UL(en)    = 0.5d0 * ( VL(i,j,ro)*v2 + B2 ) + r1*VL(i,j,pr)
!    VL --> FL
     FL(ro)  = VL(i,j,ro) * VL(i,j,vn)
     FL(mn)  = VL(i,j,ro)*VL(i,j,vn)*VL(i,j,vn)  - VL(i,j,bn)*VL(i,j,bn) + VL(i,j,pr) + 0.5d0*B2
     FL(mt1) = VL(i,j,ro)*VL(i,j,vn)*VL(i,j,vt1) - VL(i,j,bn)*VL(i,j,bt1)
     FL(mt2) = VL(i,j,ro)*VL(i,j,vn)*VL(i,j,vt2) - VL(i,j,bn)*VL(i,j,bt2)
     FL(en)  = ( 0.5d0*VL(i,j,ro)*v2 + r2*VL(i,j,pr) + B2 ) * VL(i,j,vn) - vBL * VL(i,j,bn)
     FL(bt1) = VL(i,j,vn)*VL(i,j,bt1)-VL(i,j,bn)*VL(i,j,vt1)
     FL(bt2) = VL(i,j,vn)*VL(i,j,bt2)-VL(i,j,bn)*VL(i,j,vt2)

!    VL --> fast mode^2, total pressure
     f1 = gamma * VL(i,j,pr)      ! f1: gamma p
     f2 = 4 * f1 * VL(i,j,bn)**2  ! f2: 4 gamma p B_n^2
!    vfL = sqrt( ( (f1+B2) + sqrt( (f1+B2)**2 - f2 )) / ( 2*VL(i,j,ro) ))
     aL  = ( (f1+B2) + sqrt(max( (f1+B2)**2-f2, 0.d0 ))) / ( 2*VL(i,j,ro) )
     ptL = VL(i,j,pr) + 0.5d0*B2

!    VR --> UR
     v2 = dot_product( VR(i,j,vx:vz), VR(i,j,vx:vz) )
     B2 = dot_product( VR(i,j,bx:bz), VR(i,j,bx:bz) )
     vBR= dot_product( VR(i,j,vx:vz), VR(i,j,bx:bz) )
     UR(mx:mz) = VR(i,j,ro) * VR(i,j,vx:vz)
     UR(en)    = 0.5d0 * ( VR(i,j,ro)*v2 + B2 ) + r1*VR(i,j,pr)
!    VR --> FR
     FR(ro)  = VR(i,j,ro) * VR(i,j,vn)
     FR(mn)  = VR(i,j,ro)*VR(i,j,vn)*VR(i,j,vn)  - VR(i,j,bn)*VR(i,j,bn) + VR(i,j,pr) + 0.5d0*B2
     FR(mt1) = VR(i,j,ro)*VR(i,j,vn)*VR(i,j,vt1) - VR(i,j,bn)*VR(i,j,bt1)
     FR(mt2) = VR(i,j,ro)*VR(i,j,vn)*VR(i,j,vt2) - VR(i,j,bn)*VR(i,j,bt2)
     FR(en)  = ( 0.5d0*VR(i,j,ro)*v2 + r2*VR(i,j,pr) + B2 ) * VR(i,j,vn) - vBR * VR(i,j,bn)
     FR(bt1) = VR(i,j,vn)*VR(i,j,bt1)-VR(i,j,bn)*VR(i,j,vt1)
     FR(bt2) = VR(i,j,vn)*VR(i,j,bt2)-VR(i,j,bn)*VR(i,j,vt2)

!    VR --> fast mode^2, total pressure
     f1 = gamma * VR(i,j,pr)      ! f1: gamma p
     f2 = 4 * f1 * VR(i,j,bn)**2  ! f2: 4 gamma p B_n^2
!    vfR = sqrt( ( (f1+B2) + sqrt( (f1+B2)**2 - f2 )) / ( 2*VR(i,j,ro) ))
     aR  = ( (f1+B2) + sqrt(max( (f1+B2)**2-f2, 0.d0 ))) / ( 2*VR(i,j,ro) )
     ptR = VR(i,j,pr) + 0.5d0*B2

!    Riemann fan speed (MK05 eq. 67)
!     aL = min( VL(i,j,vn) - vfL, VR(i,j,vn) - vfR )
!     aR = max( VL(i,j,vn) + vfL, VR(i,j,vn) + vfR )
!     aL = min( VL(i,j,vn), VR(i,j,vn) ) - max( vfL, vfR )
!     aR = max( VL(i,j,vn), VR(i,j,vn) ) + max( vfL, vfR )
     f1 = sqrt( max( aL, aR ) )  ! faster fast-wave (This aL [aR] stores vfL^2 [vfR^2])
     aL = min( min(VL(i,j,vn),VR(i,j,vn))-f1, 0.d0 ) ! *** if (aL > 0), then F = F(L) ***
     aR = max( max(VL(i,j,vn),VR(i,j,vn))+f1, 0.d0 ) ! *** if (aR < 0), then F = F(R) ***

!!    F = F(L)
!     if ( aL >= 0 ) then
!        F(i,j,:) = FL(:)
!!    F = F(R)
!     elseif ( aR <= 0 ) then
!        F(i,j,:) = FR(:)
!!    F = F(HLLC)
!     else

     ! HLL
     f1 = 1.d0 / ( aR - aL )
     f2 = aL * aR
     U_hll(mx:mz) = f1*( aR*UR(mx:mz) - aL*UL(mx:mz) - FR(mx:mz) + FL(mx:mz) )
     U_hll(ro:bz) = f1*( aR*VR(i,j,ro:bz) - aL*VL(i,j,ro:bz) - FR(ro:bz) + FL(ro:bz) )
     F(i,j,bt1)   = f1*( aR*FL(bt1) - aL*FR(bt1) + f2*(VR(i,j,bt1)-VL(i,j,bt1)) )
     F(i,j,bt2)   = f1*( aR*FL(bt2) - aL*FR(bt2) + f2*(VR(i,j,bt2)-VL(i,j,bt2)) )

!    entropy wave
     f2  = 1 / U_hll(ro)
     aM  = U_hll(mn ) * f2 ! vx_HLL
     at1 = U_hll(mt1) * f2 ! vy_HLL
     at2 = U_hll(mt2) * f2 ! vz_HLL
!    Total pressure (these two are identical)
     pt = ptL + VL(i,j,ro) * ( aL - VL(i,j,vn) ) * ( aM - VL(i,j,vn) )
!    pt = ptR + VR(i,j,ro) * ( aR - VR(i,j,vn) ) * ( aM - VR(i,j,vn) )

!!    F = F(L*) or F(L)
!     if ( aM >= 0 ) then

     roL = VL(i,j,ro) * ( aL - VL(i,j,vn) ) / ( aL - aM )
     enL = ( ( aL - VL(i,j,vn) )*UL(en) - ptL*VL(i,j,vn) + pt * aM + &
          U_hll(bn)*( vBL - aM*U_hll(bn) - at1*U_hll(bt1) - at2*U_hll(bt2) ) ) / ( aL - aM )
!     F(i,j,mn) = aL *( roL*aM - UL(mn) ) + FL(mn)
!     F(i,j,mt1)= aL *( roL*at1- UL(mt1)) + FL(mt1)
!     F(i,j,mt2)= aL *( roL*at2- UL(mt2)) + FL(mt2)
!     F(i,j,en) = aL *( enL    - UL(en) ) + FL(en)
!     F(i,j,ro) = aL *( roL    - VL(i,j,ro) ) + FL(ro)

!!    F = F(R*) or F(R)
!     else

     roR = VR(i,j,ro) * ( aR - VR(i,j,vn) ) / ( aR - aM )
     enR = ( ( aR - VR(i,j,vn) )*UR(en) - ptR*VR(i,j,vn) + pt * aM + &
          U_hll(bn)*( vBR - aM*U_hll(bn) - at1*U_hll(bt1) - at2*U_hll(bt2) ) ) / ( aR - aM )
!     F(i,j,mn) = aR *( roR*aM - UR(mn) ) + FR(mn)
!     F(i,j,mt1)= aR *( roR*at1- UR(mt1)) + FR(mt1)
!     F(i,j,mt2)= aR *( roR*at2- UR(mt2)) + FR(mt2)
!     F(i,j,en) = aR *( enR    - UR(en) ) + FR(en)
!     F(i,j,ro) = aR *( roR    - VR(i,j,ro) ) + FR(ro)

!    Weight factor, 0 or 1.  This looks tricky.
!    The code runs 1.0x times slower on Intel, but it runs 1.6 times faster on SPARC.
     f1 = 0.5d0 + sign(0.5d0,aM)  !!  F = F(L*) or F(L)
     f2 = 1.d0 - f1               !!  F = F(R*) or F(R)

     F(i,j,mn)  = f1*( aL *( roL*aM - UL(mn) ) + FL(mn) ) &
                + f2*( aR *( roR*aM - UR(mn) ) + FR(mn) )
     F(i,j,mt1) = f1*( aL *( roL*at1- UL(mt1)) + FL(mt1)) &
                + f2*( aR *( roR*at1- UR(mt1)) + FR(mt1))
     F(i,j,mt2) = f1*( aL *( roL*at2- UL(mt2)) + FL(mt2)) &
                + f2*( aR *( roR*at2- UR(mt2)) + FR(mt2))
     F(i,j,en)  = f1*( aL *( enL    - UL(en) ) + FL(en) ) &
                + f2*( aR *( enR    - UR(en) ) + FR(en) )
     F(i,j,ro)  = f1*( aL *( roL    - VL(i,j,ro) ) + FL(ro) ) &
                + f2*( aR *( roR    - VR(i,j,ro) ) + FR(ro) )

!     endif

  enddo
  enddo
!$omp end parallel do

  !-----------------------------------------------------------------------
  !   HLLD solver
  !-----------------------------------------------------------------------
  case(3)
!$omp parallel do &
!$omp private(i,j,v2,B2,vBL,vBR,UL,UR,FL,FR,f1,f2,aL,aL1,aM,aR1,aR,aVL,aVR) &
!$omp private(ptL,ptR,UL1,UR1,U2,U_hll,pt) &
!$omp private(roL,roR,vt1L,vt1R,vt2L,vt2R,at1,at2,roLs,roRs)
  do j=js,je
  do i=is,ie

!    VL --> UL
     v2 = dot_product( VL(i,j,vx:vz), VL(i,j,vx:vz) )
     B2 = dot_product( VL(i,j,bx:bz), VL(i,j,bx:bz) )
     vBL= dot_product( VL(i,j,vx:vz), VL(i,j,bx:bz) )
     UL(mx:mz) = VL(i,j,ro) * VL(i,j,vx:vz)
     UL(en)    = 0.5d0 * ( VL(i,j,ro)*v2 + B2 ) + r1*VL(i,j,pr)
!    VL --> FL
     FL(ro)  = VL(i,j,ro) * VL(i,j,vn)
     FL(mn)  = VL(i,j,ro)*VL(i,j,vn)*VL(i,j,vn)  - VL(i,j,bn)*VL(i,j,bn) + VL(i,j,pr) + 0.5d0*B2
     FL(mt1) = VL(i,j,ro)*VL(i,j,vn)*VL(i,j,vt1) - VL(i,j,bn)*VL(i,j,bt1)
     FL(mt2) = VL(i,j,ro)*VL(i,j,vn)*VL(i,j,vt2) - VL(i,j,bn)*VL(i,j,bt2)
     FL(en)  = ( 0.5d0*VL(i,j,ro)*v2 + r2*VL(i,j,pr) + B2 ) * VL(i,j,vn) - vBL * VL(i,j,bn)
     FL(bt1) = VL(i,j,vn)*VL(i,j,bt1)-VL(i,j,bn)*VL(i,j,vt1)
     FL(bt2) = VL(i,j,vn)*VL(i,j,bt2)-VL(i,j,bn)*VL(i,j,vt2)

!    VL --> fast mode^2, total pressure
     f1 = gamma * VL(i,j,pr)      ! f1: gamma p
     f2 = 4 * f1 * VL(i,j,bn)**2  ! f2: 4 gamma p B_n^2
!    vfL = sqrt( ( (f1+B2) + sqrt( (f1+B2)**2 - f2 )) / ( 2*VL(i,j,ro) ))
     aL  = ( (f1+B2) + sqrt(max( (f1+B2)**2-f2, 0.d0 ))) / ( 2*VL(i,j,ro) )
     ptL = VL(i,j,pr) + 0.5d0*B2

!    VR --> UR
     v2 = dot_product( VR(i,j,vx:vz), VR(i,j,vx:vz) )
     B2 = dot_product( VR(i,j,bx:bz), VR(i,j,bx:bz) )
     vBR= dot_product( VR(i,j,vx:vz), VR(i,j,bx:bz) )
     UR(mx:mz) = VR(i,j,ro) * VR(i,j,vx:vz)
     UR(en)    = 0.5d0 * ( VR(i,j,ro)*v2 + B2 ) + r1*VR(i,j,pr)
!    VR --> FR
     FR(ro)  = VR(i,j,ro) * VR(i,j,vn)
     FR(mn)  = VR(i,j,ro)*VR(i,j,vn)*VR(i,j,vn)  - VR(i,j,bn)*VR(i,j,bn) + VR(i,j,pr) + 0.5d0*B2
     FR(mt1) = VR(i,j,ro)*VR(i,j,vn)*VR(i,j,vt1) - VR(i,j,bn)*VR(i,j,bt1)
     FR(mt2) = VR(i,j,ro)*VR(i,j,vn)*VR(i,j,vt2) - VR(i,j,bn)*VR(i,j,bt2)
     FR(en)  = ( 0.5d0*VR(i,j,ro)*v2 + r2*VR(i,j,pr) + B2 ) * VR(i,j,vn) - vBR * VR(i,j,bn)
     FR(bt1) = VR(i,j,vn)*VR(i,j,bt1)-VR(i,j,bn)*VR(i,j,vt1)
     FR(bt2) = VR(i,j,vn)*VR(i,j,bt2)-VR(i,j,bn)*VR(i,j,vt2)

!    VR --> fast mode^2, total pressure
     f1 = gamma * VR(i,j,pr)      ! f1: gamma p
     f2 = 4 * f1 * VR(i,j,bn)**2  ! f2: 4 gamma p B_n^2
!    vfR = sqrt( ( (f1+B2) + sqrt( (f1+B2)**2 - f2 )) / ( 2*VR(i,j,ro) ))
     aR  = ( (f1+B2) + sqrt(max( (f1+B2)**2-f2, 0.d0 ))) / ( 2*VR(i,j,ro) )
     ptR = VR(i,j,pr) + 0.5d0*B2

!    Riemann fan speed (MK05 eq. 67)
!     aL = min( VL(i,j,vn) - vfL, VR(i,j,vn) - vfR )
!     aR = max( VL(i,j,vn) + vfL, VR(i,j,vn) + vfR )
!     aL = min( VL(i,j,vn), VR(i,j,vn) ) - max( vfL, vfR )
!     aR = max( VL(i,j,vn), VR(i,j,vn) ) + max( vfL, vfR )
     f1 = sqrt( max( aL, aR ) )  ! faster fast-wave (This aL [aR] stores vfL^2 [vfR^2])
     aL = min( VL(i,j,vn), VR(i,j,vn) ) - f1
     aR = max( VL(i,j,vn), VR(i,j,vn) ) + f1

!    F = F(L)
     if ( aL >= 0 ) then
        F(i,j,:) = FL(:)
!    F = F(R)
     elseif ( aR <= 0 ) then
        F(i,j,:) = FR(:)

!    HLLC/HLLD flux
     else

!       HLL state
        f1 = 1.d0 / ( aR - aL )
        U_hll(mx:en) = f1*( aR*UR(mx:en) - aL*UL(mx:en) - FR(mx:en) + FL(mx:en) )
        U_hll(ro:bz) = f1*( aR*VR(i,j,ro:bz) - aL*VL(i,j,ro:bz) - FR(ro:bz) + FL(ro:bz) )

!       entropy wave
        aM = U_hll(mn) / U_hll(ro)
        
!       Total pressure
        pt  = ptL + VL(i,j,ro) * ( aL - VL(i,j,vn) ) * ( aM - VL(i,j,vn) )
        roL = VL(i,j,ro) * ( aL - VL(i,j,vn) ) / ( aL - aM ) ! ro(L*)
        roR = VR(i,j,ro) * ( aR - VR(i,j,vn) ) / ( aR - aM ) ! ro(R*)

!       For logical consistency, we employ Bn_hll as B_x in the intermediate states
        aVL = abs( U_hll(bn) )/sqrt( roL )  ! Alfven wave (L)
        aVR = abs( U_hll(bn) )/sqrt( roR )  ! Alfven wave (R)

! ========== revert to HLLC-G ==========
!       When SR*(SL*) --> SR(SL), we switch to a HLLC solver.
!       In these limits, the HLLD denominator  rho_L (SL-uL)(SL-SM) - B_x^2,
!       which can be written as  rho_L* (SL+SL*-2SM) (SL-SL*), becomes zero.
!       Since the two HLLC states L* and R* are relevant to L** and R** in HLLD,
!       we should employ HLLC-G method for logical consistency:
!       i.e.:  vy_L* = vy_R* (HLLC-G)  <==>  vy_L** = vy_R** (HLLD).
!       if ( .true. ) then
        if ( ( aL >= ( aM - hllg_factor*aVL ) ) .or. &
             ( ( aM + hllg_factor*aVR ) >= aR ) ) then

           at1 = U_hll(mt1) / U_hll(ro) ! vt1_HLL
           at2 = U_hll(mt2) / U_hll(ro) ! vt2_HLL

!          F = F(L*)
           if ( aM >= 0 ) then

              U2(en) = ( ( aL-VL(i,j,vn) )*UL(en) - ptL*VL(i,j,vn) + pt*aM + &
                   U_hll(bn)*( vBL - aM*U_hll(bn) - at1*U_hll(bt1) - at2*U_hll(bt2) ) ) / ( aL - aM )
              U2(mn)  = roL * aM
              U2(mt1) = roL * at1
              U2(mt2) = roL * at2

              F(i,j,mx:en) = aL *( U2(mx:en)  - UL(mx:en)   ) + FL(mx:en)
              F(i,j,ro)    = aL *( roL        - VL(i,j,ro)  ) + FL(ro)
              F(i,j,bt1)   = aL *( U_hll(bt1) - VL(i,j,bt1) ) + FL(bt1)
              F(i,j,bt2)   = aL *( U_hll(bt2) - VL(i,j,bt2) ) + FL(bt2)

!          F = F(R*)
           else

              U2(en) = ( ( aR-VR(i,j,vn) )*UR(en) - ptR*VR(i,j,vn) + pt*aM + &
                   U_hll(bn)*( vBR - aM*U_hll(bn) - at1*U_hll(bt1) - at2*U_hll(bt2) ) ) / ( aR - aM )
              U2(mn)  = roR * aM
              U2(mt1) = roR * at1
              U2(mt2) = roR * at2

              F(i,j,mx:en) = aR *( U2(mx:en)  - UR(mx:en)   ) + FR(mx:en)
              F(i,j,ro)    = aR *( roR        - VR(i,j,ro)  ) + FR(ro)
              F(i,j,bt1)   = aR *( U_hll(bt1) - VR(i,j,bt1) ) + FR(bt1)
              F(i,j,bt2)   = aR *( U_hll(bt2) - VR(i,j,bt2) ) + FR(bt2)

           endif

! ========== HLLD flux ==========
        else

!          Intermediate state (L*)
           f1 = 1.d0 / ( VL(i,j,ro)*(aL-VL(i,j,vn))*(aL-aM) - U_hll(bn)**2 ) ! HLLD denominator
           UL1(bt1) = VL(i,j,bt1) * f1 * ( VL(i,j,ro)*(aL-VL(i,j,vn))**2 - U_hll(bn)**2 )
           UL1(bt2) = VL(i,j,bt2) * f1 * ( VL(i,j,ro)*(aL-VL(i,j,vn))**2 - U_hll(bn)**2 )
           vt1L     = VL(i,j,vt1) - f1 * U_hll(bn)*VL(i,j,bt1)*(aM-VL(i,j,vn))
           vt2L     = VL(i,j,vt2) - f1 * U_hll(bn)*VL(i,j,bt2)*(aM-VL(i,j,vn))
           
           UL1(ro)  = roL
           UL1(en)  = ( ( aL-VL(i,j,vn) )*UL(en) - ptL*VL(i,j,vn) + pt*aM + &
                U_hll(bn)*( vBL - aM*U_hll(bn) - vt1L*UL1(bt1) - vt2L*UL1(bt2)) ) / ( aL - aM )
           UL1(mn)  = roL * aM
           UL1(mt1) = roL * vt1L
           UL1(mt2) = roL * vt2L

!          Intermediate state (R*)
           f1 = 1.d0 / ( VR(i,j,ro)*(aR-VR(i,j,vn))*(aR-aM) - U_hll(bn)**2 ) ! HLLD denominator
           UR1(bt1) = VR(i,j,bt1) * f1 * ( VR(i,j,ro)*(aR-VR(i,j,vn))**2 - U_hll(bn)**2 )
           UR1(bt2) = VR(i,j,bt2) * f1 * ( VR(i,j,ro)*(aR-VR(i,j,vn))**2 - U_hll(bn)**2 )
           vt1R     = VR(i,j,vt1) - f1 * U_hll(bn)*VR(i,j,bt1)*(aM-VR(i,j,vn))
           vt2R     = VR(i,j,vt2) - f1 * U_hll(bn)*VR(i,j,bt2)*(aM-VR(i,j,vn))
           
           UR1(ro)  = roR
           UR1(en)  = ( ( aR-VR(i,j,vn) )*UR(en) - ptR*VR(i,j,vn) + pt*aM + &
                U_hll(bn)*( vBR - aM*U_hll(bn) - vt1R*UR1(bt1) - vt2R*UR1(bt2)) ) / ( aR - aM )
           UR1(mn)  = roR * aM
           UR1(mt1) = roR * vt1R
           UR1(mt2) = roR * vt2R
           
!          rotational waves
           aL1 = aM - aVL
           aR1 = aM + aVR

!          F = F(L*)
           if( aL1 >= 0 ) then

              F(i,j,mx:en) = aL *( UL1(mx:en) - UL(mx:en)   ) + FL(mx:en)
              F(i,j,ro)    = aL *( UL1(ro)    - VL(i,j,ro)  ) + FL(ro)
              F(i,j,bt1)   = aL *( UL1(bt1)   - VL(i,j,bt1) ) + FL(bt1)
              F(i,j,bt2)   = aL *( UL1(bt2)   - VL(i,j,bt2) ) + FL(bt2)

!          F = F(R*)
           elseif( aR1 <= 0 ) then

              F(i,j,mx:en) = aR *( UR1(mx:en) - UR(mx:en)   ) + FR(mx:en)
              F(i,j,ro)    = aR *( UR1(ro)    - VR(i,j,ro)  ) + FR(ro)
              F(i,j,bt1)   = aR *( UR1(bt1)   - VR(i,j,bt1) ) + FR(bt1)
              F(i,j,bt2)   = aR *( UR1(bt2)   - VR(i,j,bt2) ) + FR(bt2)

!          Central states are tricky
!          Question: can we really rule out U_hll(bn) == 0 ?
           else

              roLs = sqrt( roL )  ! sqrt(ro(L*))
              roRs = sqrt( roR )  ! sqrt(ro(R*))
              f1   = 1.d0 / ( roLs + roRs )
              f2   = sign( 1.d0, U_hll(bn) )
              at1  = f1 * ( roLs*vt1L + roRs*vt1R + ( UR1(bt1)-UL1(bt1) )*f2 )
              at2  = f1 * ( roLs*vt2L + roRs*vt2R + ( UR1(bt2)-UL1(bt2) )*f2 )
              U2(bt1) = f1 * ( roLs*UR1(bt1) + roRs*UL1(bt1) + roLs*roRs*( vt1R-vt1L )*f2 )
              U2(bt2) = f1 * ( roLs*UR1(bt2) + roRs*UL1(bt2) + roLs*roRs*( vt2R-vt2L )*f2 )

!             F = F(L**)
              if( aM >= 0 ) then

                 U2(ro) = roL
                 U2(en) = UL1(en) - roLs * ( vt1L*UL1(bt1) + vt2L*UL1(bt2) &
                      - at1*U2(bt1) - at2*U2(bt2) ) * f2
                 U2(mn)  = roL * aM
                 U2(mt1) = roL * at1
                 U2(mt2) = roL * at2

                 F(i,j,mx:en) = aL1*U2(mx:en) - (aL1-aL)*UL1(mx:en) - aL*UL(mx:en) + FL(mx:en)
                 F(i,j,ro)    = aL1*U2(ro)  - (aL1-aL)*UL1(ro)  - aL*VL(i,j,ro)  + FL(ro)
                 F(i,j,bt1)   = aL1*U2(bt1) - (aL1-aL)*UL1(bt1) - aL*VL(i,j,bt1) + FL(bt1)
                 F(i,j,bt2)   = aL1*U2(bt2) - (aL1-aL)*UL1(bt2) - aL*VL(i,j,bt2) + FL(bt2)

!             F = F(R**)
              else

                 U2(ro) = roR
                 U2(en) = UR1(en) + roRs * ( vt1R*UR1(bt1) + vt2R*UR1(bt2) &
                      - at1*U2(bt1) - at2*U2(bt2) ) * f2
                 U2(mn)  = roR * aM
                 U2(mt1) = roR * at1
                 U2(mt2) = roR * at2

                 F(i,j,mx:en) = aR1*U2(mx:en) - (aR1-aR)*UR1(mx:en) - aR*UR(mx:en) + FR(mx:en)
                 F(i,j,ro)    = aR1*U2(ro)  - (aR1-aR)*UR1(ro)  - aR*VR(i,j,ro)  + FR(ro)
                 F(i,j,bt1)   = aR1*U2(bt1) - (aR1-aR)*UR1(bt1) - aR*VR(i,j,bt1) + FR(bt1)
                 F(i,j,bt2)   = aR1*U2(bt2) - (aR1-aR)*UR1(bt2) - aR*VR(i,j,bt2) + FR(bt2)

              endif

           endif

        endif

     endif

  enddo
  enddo
!$omp end parallel do

  endselect
  !-----------------------------------------------------------------------

  return
end subroutine flux_solver
