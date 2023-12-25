subroutine flux_solver(F,VL,VR,ch,ix,jx,kx,dir,type)
!-----------------------------------------------------------------------
!     HLLC-G flux solver
!       Ref: K. F. Gurski, SIAM J. Sci. Comput., 25, 2165 (2004)
!     HLLD flux solver
!       Ref: T. Miyoshi, K. Kusano, J. Comput. Phys., 208, 315 (2005)
!-----------------------------------------------------------------------
!     2010/05/11  S. Zenitani  HLLC-G solver
!     2010/05/12  S. Zenitani  HLLD solver
!     2016/09/06  S. Zenitani  X/Y/Z directions
!     2023/12/25  S. Zenitani  3-D version
!-----------------------------------------------------------------------
  implicit none
  include 'param.h'
!-----------------------------------------------------------------------
! size of arrays [input]
  integer, intent(in)  :: ix, jx, kx
! numerical flux (F) [output]
  real(8), intent(out) :: F(ix,jx,kx,var1)
! left/right states (VL,VR) [input]
  real(8), intent(in)  :: VL(ix,jx,kx,var1), VR(ix,jx,kx,var1)
! div cleaning wave speed
  real(8), intent(in) :: ch
  real(8) :: ch2
! direction [input]: 1 (X), 2 (Y), 3 (Z)
  integer, intent(in)  :: dir
! numerical flux [input]:
! 0 (LLF Rusanov), 1 (HLL), 2 (HLLC), 3 (HLLD), 4 (LHLLD experimental)
  integer, intent(in)  :: type
!-----------------------------------------------------------------------
! left/right conserved variables (UL & UR; local)
  real(8) :: UL(var2), UR(var2)
! left/right numerical flux (FL & FR; local)
  real(8) :: FL(var1), FR(var1)
  integer :: i, j, k, is=0, ie=0, js=0, je=0, ks=0, ke=0
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
  ch2 = ch**2

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
     ks = min(2,kx); ke = max(1,kx-1)
  case(2)
     is = min(2,ix); ie = max(1,ix-1)
     js = 1; je = jx-1
     ks = min(2,kx); ke = max(1,kx-1)
  case(3)
     is = min(2,ix); ie = max(1,ix-1)
     js = min(2,jx); je = max(1,jx-1)
     ks = 1; ke = kx-1
  end select

! main switch
  select case(type)
  !-----------------------------------------------------------------------
  !  Local Lax Friedrich (Rusanov)
  !  No OpenMP optimization - just for education
  !-----------------------------------------------------------------------
  case(0)

     do k=ks,ke
     do j=js,je
     do i=is,ie

!    VL --> UL
     v2 = dot_product( VL(i,j,k,vx:vz), VL(i,j,k,vx:vz) )
     B2 = dot_product( VL(i,j,k,bx:bz), VL(i,j,k,bx:bz) )
     vBL= dot_product( VL(i,j,k,vx:vz), VL(i,j,k,bx:bz) )
     UL(mx:mz) = VL(i,j,k,ro) * VL(i,j,k,vx:vz)
     UL(en)    = 0.5d0 * ( VL(i,j,k,ro)*v2 + B2 ) + r1*VL(i,j,k,pr)
!    VL --> FL
     FL(ro)  = VL(i,j,k,ro) * VL(i,j,k,vn)
     FL(mn)  = VL(i,j,k,ro)*VL(i,j,k,vn)*VL(i,j,k,vn)  - VL(i,j,k,bn)*VL(i,j,k,bn) + VL(i,j,k,pr) + 0.5d0*B2
     FL(mt1) = VL(i,j,k,ro)*VL(i,j,k,vn)*VL(i,j,k,vt1) - VL(i,j,k,bn)*VL(i,j,k,bt1)
     FL(mt2) = VL(i,j,k,ro)*VL(i,j,k,vn)*VL(i,j,k,vt2) - VL(i,j,k,bn)*VL(i,j,k,bt2)
     FL(en)  = ( 0.5d0*VL(i,j,k,ro)*v2 + r2*VL(i,j,k,pr) + B2 ) * VL(i,j,k,vn) - vBL * VL(i,j,k,bn)
!    FL(bn)  = 0.d0
     FL(bt1) = VL(i,j,k,vn)*VL(i,j,k,bt1)-VL(i,j,k,bn)*VL(i,j,k,vt1)
     FL(bt2) = VL(i,j,k,vn)*VL(i,j,k,bt2)-VL(i,j,k,bn)*VL(i,j,k,vt2)
!    FL(ps)  = 0.d0

!    VL --> fast mode
     f1 = gamma * VL(i,j,k,pr) ! f1: gamma p
     f2 = 4 * VL(i,j,k,bn)**2  ! f2: 4 B_n^2
     aL = sqrt( ( (f1+B2) + sqrt( (f1+B2)**2 - f1*f2 )) / ( 2*VL(i,j,k,ro) ))

!    VR --> UR
     v2 = dot_product( VR(i,j,k,vx:vz), VR(i,j,k,vx:vz) )
     B2 = dot_product( VR(i,j,k,bx:bz), VR(i,j,k,bx:bz) )
     vBR= dot_product( VR(i,j,k,vx:vz), VR(i,j,k,bx:bz) )
     UR(mx:mz) = VR(i,j,k,ro) * VR(i,j,k,vx:vz)
     UR(en)    = 0.5d0 * ( VR(i,j,k,ro)*v2 + B2 ) + r1*VR(i,j,k,pr)
!    VR --> FR
     FR(ro)  = VR(i,j,k,ro) * VR(i,j,k,vn)
     FR(mn)  = VR(i,j,k,ro)*VR(i,j,k,vn)*VR(i,j,k,vn)  - VR(i,j,k,bn)*VR(i,j,k,bn) + VR(i,j,k,pr) + 0.5d0*B2
     FR(mt1) = VR(i,j,k,ro)*VR(i,j,k,vn)*VR(i,j,k,vt1) - VR(i,j,k,bn)*VR(i,j,k,bt1)
     FR(mt2) = VR(i,j,k,ro)*VR(i,j,k,vn)*VR(i,j,k,vt2) - VR(i,j,k,bn)*VR(i,j,k,bt2)
     FR(en)  = ( 0.5d0*VR(i,j,k,ro)*v2 + r2*VR(i,j,k,pr) + B2 ) * VR(i,j,k,vn) - vBR * VR(i,j,k,bn)
!    FR(bn)  = 0.d0
     FR(bt1) = VR(i,j,k,vn)*VR(i,j,k,bt1)-VR(i,j,k,bn)*VR(i,j,k,vt1)
     FR(bt2) = VR(i,j,k,vn)*VR(i,j,k,bt2)-VR(i,j,k,bn)*VR(i,j,k,vt2)
!    FR(ps)  = 0.d0

!    VR --> fast mode
     f1 = gamma * VR(i,j,k,pr) ! f1: gamma p
     f2 = 4 * VR(i,j,k,bn)**2  ! f2: 4 B_n^2
     aR = sqrt( ( (f1+B2) + sqrt( (f1+B2)**2 - f1*f2 )) / ( 2*VR(i,j,k,ro) ))

     f1 = max( abs(VL(i,j,k,vn))+aL, abs(VR(i,j,k,vn))+aR )
     F(i,j,k,mx:en) = 0.5d0 * ( FL(mx:en)+FR(mx:en) - f1*( UR(mx:en)-UL(mx:en) ) )
     F(i,j,k,ro)    = 0.5d0 * ( FL(ro)   +FR(ro)    - f1*( VR(i,j,k,ro) -VL(i,j,k,ro)  ) )
     F(i,j,k,bt1)   = 0.5d0 * ( FL(bt1)  +FR(ro)    - f1*( VR(i,j,k,bt1)-VL(i,j,k,bt1) ) )
     F(i,j,k,bt2)   = 0.5d0 * ( FL(bt2)  +FR(ro)    - f1*( VR(i,j,k,bt2)-VL(i,j,k,bt2) ) )

     F(i,j,k,bn) = 0.5d0*(     ( VL(i,j,k,ps)+VR(i,j,k,ps) ) - ch*( VR(i,j,k,bn)-VL(i,j,k,bn) ) )
     F(i,j,k,ps) = 0.5d0*( ch2*( VL(i,j,k,bn)+VR(i,j,k,bn) ) - ch*( VR(i,j,k,ps)-VL(i,j,k,ps) ) )

     enddo
     enddo
     enddo

  !-----------------------------------------------------------------------
  !  HLL solver
  !-----------------------------------------------------------------------
  case(1)
!$omp parallel do private(i,j,k,v2,B2,vBL,vBR,UL,UR,FL,FR,f1,f2,aL,aR)
     do k=ks,ke
     do j=js,je
     do i=is,ie

!    VL --> UL
     v2 = dot_product( VL(i,j,k,vx:vz), VL(i,j,k,vx:vz) )
     B2 = dot_product( VL(i,j,k,bx:bz), VL(i,j,k,bx:bz) )
     vBL= dot_product( VL(i,j,k,vx:vz), VL(i,j,k,bx:bz) )
     UL(mx:mz) = VL(i,j,k,ro) * VL(i,j,k,vx:vz)
     UL(en)    = 0.5d0 * ( VL(i,j,k,ro)*v2 + B2 ) + r1*VL(i,j,k,pr)
!    VL --> FL
     FL(ro)  = VL(i,j,k,ro) * VL(i,j,k,vn)
     FL(mn)  = VL(i,j,k,ro)*VL(i,j,k,vn)*VL(i,j,k,vn)  - VL(i,j,k,bn)*VL(i,j,k,bn) + VL(i,j,k,pr) + 0.5d0*B2
     FL(mt1) = VL(i,j,k,ro)*VL(i,j,k,vn)*VL(i,j,k,vt1) - VL(i,j,k,bn)*VL(i,j,k,bt1)
     FL(mt2) = VL(i,j,k,ro)*VL(i,j,k,vn)*VL(i,j,k,vt2) - VL(i,j,k,bn)*VL(i,j,k,bt2)
     FL(en)  = ( 0.5d0*VL(i,j,k,ro)*v2 + r2*VL(i,j,k,pr) + B2 ) * VL(i,j,k,vn) - vBL * VL(i,j,k,bn)
     FL(bt1) = VL(i,j,k,vn)*VL(i,j,k,bt1)-VL(i,j,k,bn)*VL(i,j,k,vt1)
     FL(bt2) = VL(i,j,k,vn)*VL(i,j,k,bt2)-VL(i,j,k,bn)*VL(i,j,k,vt2)

!    VL --> fast mode^2
     f1 = gamma * VL(i,j,k,pr) ! f1: gamma p
     f2 = 4 * VL(i,j,k,bn)**2  ! f2: 4 B_n^2
!    vfL= sqrt( ( (f1+B2) + sqrt( (f1+B2)**2 - f1*f2 )) / ( 2*VL(i,j,k,ro) ))
     aL = ( (f1+B2) + sqrt(max( (f1+B2)**2 - f1*f2, 0.d0 ))) / ( 2*VL(i,j,k,ro) )

!    VR --> UR
     v2 = dot_product( VR(i,j,k,vx:vz), VR(i,j,k,vx:vz) )
     B2 = dot_product( VR(i,j,k,bx:bz), VR(i,j,k,bx:bz) )
     vBR= dot_product( VR(i,j,k,vx:vz), VR(i,j,k,bx:bz) )
     UR(mx:mz) = VR(i,j,k,ro) * VR(i,j,k,vx:vz)
     UR(en)    = 0.5d0 * ( VR(i,j,k,ro)*v2 + B2 ) + r1*VR(i,j,k,pr)
!    VR --> FR
     FR(ro)  = VR(i,j,k,ro) * VR(i,j,k,vn)
     FR(mn)  = VR(i,j,k,ro)*VR(i,j,k,vn)*VR(i,j,k,vn)  - VR(i,j,k,bn)*VR(i,j,k,bn) + VR(i,j,k,pr) + 0.5d0*B2
     FR(mt1) = VR(i,j,k,ro)*VR(i,j,k,vn)*VR(i,j,k,vt1) - VR(i,j,k,bn)*VR(i,j,k,bt1)
     FR(mt2) = VR(i,j,k,ro)*VR(i,j,k,vn)*VR(i,j,k,vt2) - VR(i,j,k,bn)*VR(i,j,k,bt2)
     FR(en)  = ( 0.5d0*VR(i,j,k,ro)*v2 + r2*VR(i,j,k,pr) + B2 ) * VR(i,j,k,vn) - vBR * VR(i,j,k,bn)
     FR(bt1) = VR(i,j,k,vn)*VR(i,j,k,bt1)-VR(i,j,k,bn)*VR(i,j,k,vt1)
     FR(bt2) = VR(i,j,k,vn)*VR(i,j,k,bt2)-VR(i,j,k,bn)*VR(i,j,k,vt2)

!    VR --> fast mode^2
     f1 = gamma * VR(i,j,k,pr) ! f1: gamma p
     f2 = 4 * VR(i,j,k,bn)**2  ! f2: 4 B_n^2
!    vfR= sqrt( ( (f1+B2) + sqrt( (f1+B2)**2 - f1*f2 )) / ( 2*VR(i,j,k,ro) ))
     aR = ( (f1+B2) + sqrt(max( (f1+B2)**2 - f1*f2, 0.d0 ))) / ( 2*VR(i,j,k,ro) )

!    Riemann fan speed (MK05 eq. 67)
!     aL = min( VL(i,j,k,vn) - vfL, VR(i,j,k,vn) - vfR )
!     aR = max( VL(i,j,k,vn) + vfL, VR(i,j,k,vn) + vfR )
!     aL = min( VL(i,j,k,vn), VR(i,j,k,vn) ) - max( vfL, vfR )
!     aR = max( VL(i,j,k,vn), VR(i,j,k,vn) ) + max( vfL, vfR )
     f1 = sqrt( max( aL, aR ) )  ! faster fast-wave (This aL [aR] stores vfL^2 [vfR^2])
     aL = min( min(VL(i,j,k,vn),VR(i,j,k,vn))-f1, 0.d0 ) ! *** if (aL > 0), then F = F(L) ***
     aR = max( max(VL(i,j,k,vn),VR(i,j,k,vn))+f1, 0.d0 ) ! *** if (aR < 0), then F = F(R) ***

     f1 = 1.d0 / ( aR - aL )
     f2 = aL * aR
     F(i,j,k,mx:en) = f1*( aR*FL(mx:en) - aL*FR(mx:en) + f2 *(UR(mx:en)-UL(mx:en)) )
     F(i,j,k,ro)    = f1*( aR*FL(ro)  - aL*FR(ro)  + f2*(VR(i,j,k,ro) -VL(i,j,k,ro) ) )
     F(i,j,k,bt1)   = f1*( aR*FL(bt1) - aL*FR(bt1) + f2*(VR(i,j,k,bt1)-VL(i,j,k,bt1)) )
     F(i,j,k,bt2)   = f1*( aR*FL(bt2) - aL*FR(bt2) + f2*(VR(i,j,k,bt2)-VL(i,j,k,bt2)) )

     F(i,j,k,bn) = 0.5d0*(     ( VL(i,j,k,ps)+VR(i,j,k,ps) ) - ch*( VR(i,j,k,bn)-VL(i,j,k,bn) ) )
     F(i,j,k,ps) = 0.5d0*( ch2*( VL(i,j,k,bn)+VR(i,j,k,bn) ) - ch*( VR(i,j,k,ps)-VL(i,j,k,ps) ) )

     enddo
     enddo
     enddo
!$omp end parallel do

  !-----------------------------------------------------------------------
  !  HLLC solver
  !-----------------------------------------------------------------------
  case(2)
!$omp parallel do private(i,j,k,v2,B2,vBL,vBR,UL,UR,FL,FR,f1,f2) &
!$omp private(aL,aR,aM,at1,at2,ptL,ptR,pt,roL,roR,enL,enR,U_hll)
     do k=ks,ke
     do j=js,je
     do i=is,ie

!    VL --> UL
     v2 = dot_product( VL(i,j,k,vx:vz), VL(i,j,k,vx:vz) )
     B2 = dot_product( VL(i,j,k,bx:bz), VL(i,j,k,bx:bz) )
     vBL= dot_product( VL(i,j,k,vx:vz), VL(i,j,k,bx:bz) )
     UL(mx:mz) = VL(i,j,k,ro) * VL(i,j,k,vx:vz)
     UL(en)    = 0.5d0 * ( VL(i,j,k,ro)*v2 + B2 ) + r1*VL(i,j,k,pr)
!    VL --> FL
     FL(ro)  = VL(i,j,k,ro) * VL(i,j,k,vn)
     FL(mn)  = VL(i,j,k,ro)*VL(i,j,k,vn)*VL(i,j,k,vn)  - VL(i,j,k,bn)*VL(i,j,k,bn) + VL(i,j,k,pr) + 0.5d0*B2
     FL(mt1) = VL(i,j,k,ro)*VL(i,j,k,vn)*VL(i,j,k,vt1) - VL(i,j,k,bn)*VL(i,j,k,bt1)
     FL(mt2) = VL(i,j,k,ro)*VL(i,j,k,vn)*VL(i,j,k,vt2) - VL(i,j,k,bn)*VL(i,j,k,bt2)
     FL(en)  = ( 0.5d0*VL(i,j,k,ro)*v2 + r2*VL(i,j,k,pr) + B2 ) * VL(i,j,k,vn) - vBL * VL(i,j,k,bn)
     FL(bt1) = VL(i,j,k,vn)*VL(i,j,k,bt1)-VL(i,j,k,bn)*VL(i,j,k,vt1)
     FL(bt2) = VL(i,j,k,vn)*VL(i,j,k,bt2)-VL(i,j,k,bn)*VL(i,j,k,vt2)

!    VL --> fast mode^2, total pressure
     f1 = gamma * VL(i,j,k,pr) ! f1: gamma p
     f2 = 4 * VL(i,j,k,bn)**2  ! f2: 4 B_n^2
!    vfL = sqrt( ( (f1+B2) + sqrt( (f1+B2)**2 - f1*f2 )) / ( 2*VL(i,j,k,ro) ))
     aL  = ( (f1+B2) + sqrt(max( (f1+B2)**2 - f1*f2, 0.d0 ))) / ( 2*VL(i,j,k,ro) )
     ptL = VL(i,j,k,pr) + 0.5d0*B2

!    VR --> UR
     v2 = dot_product( VR(i,j,k,vx:vz), VR(i,j,k,vx:vz) )
     B2 = dot_product( VR(i,j,k,bx:bz), VR(i,j,k,bx:bz) )
     vBR= dot_product( VR(i,j,k,vx:vz), VR(i,j,k,bx:bz) )
     UR(mx:mz) = VR(i,j,k,ro) * VR(i,j,k,vx:vz)
     UR(en)    = 0.5d0 * ( VR(i,j,k,ro)*v2 + B2 ) + r1*VR(i,j,k,pr)
!    VR --> FR
     FR(ro)  = VR(i,j,k,ro) * VR(i,j,k,vn)
     FR(mn)  = VR(i,j,k,ro)*VR(i,j,k,vn)*VR(i,j,k,vn)  - VR(i,j,k,bn)*VR(i,j,k,bn) + VR(i,j,k,pr) + 0.5d0*B2
     FR(mt1) = VR(i,j,k,ro)*VR(i,j,k,vn)*VR(i,j,k,vt1) - VR(i,j,k,bn)*VR(i,j,k,bt1)
     FR(mt2) = VR(i,j,k,ro)*VR(i,j,k,vn)*VR(i,j,k,vt2) - VR(i,j,k,bn)*VR(i,j,k,bt2)
     FR(en)  = ( 0.5d0*VR(i,j,k,ro)*v2 + r2*VR(i,j,k,pr) + B2 ) * VR(i,j,k,vn) - vBR * VR(i,j,k,bn)
     FR(bt1) = VR(i,j,k,vn)*VR(i,j,k,bt1)-VR(i,j,k,bn)*VR(i,j,k,vt1)
     FR(bt2) = VR(i,j,k,vn)*VR(i,j,k,bt2)-VR(i,j,k,bn)*VR(i,j,k,vt2)

!    VR --> fast mode^2, total pressure
     f1 = gamma * VR(i,j,k,pr) ! f1: gamma p
     f2 = 4 * VR(i,j,k,bn)**2  ! f2: 4 B_n^2
!    vfR = sqrt( ( (f1+B2) + sqrt( (f1+B2)**2 - f1*f2 )) / ( 2*VR(i,j,k,ro) ))
     aR  = ( (f1+B2) + sqrt(max( (f1+B2)**2 - f1*f2, 0.d0 ))) / ( 2*VR(i,j,k,ro) )
     ptR = VR(i,j,k,pr) + 0.5d0*B2

!    Riemann fan speed (MK05 eq. 67)
!     aL = min( VL(i,j,k,vn) - vfL, VR(i,j,k,vn) - vfR )
!     aR = max( VL(i,j,k,vn) + vfL, VR(i,j,k,vn) + vfR )
!     aL = min( VL(i,j,k,vn), VR(i,j,k,vn) ) - max( vfL, vfR )
!     aR = max( VL(i,j,k,vn), VR(i,j,k,vn) ) + max( vfL, vfR )
     f1 = sqrt( max( aL, aR ) )  ! faster fast-wave (This aL [aR] stores vfL^2 [vfR^2])
     aL = min( min(VL(i,j,k,vn),VR(i,j,k,vn))-f1, 0.d0 ) ! *** if (aL > 0), then F = F(L) ***
     aR = max( max(VL(i,j,k,vn),VR(i,j,k,vn))+f1, 0.d0 ) ! *** if (aR < 0), then F = F(R) ***

     ! HLL
     f1 = 1.d0 / ( aR - aL )
     f2 = aL * aR
     U_hll(mx:mz) = f1*( aR*UR(mx:mz) - aL*UL(mx:mz) - FR(mx:mz) + FL(mx:mz) )
     U_hll(ro:bz) = f1*( aR*VR(i,j,k,ro:bz) - aL*VL(i,j,k,ro:bz) - FR(ro:bz) + FL(ro:bz) )
     F(i,j,k,bt1)   = f1*( aR*FL(bt1) - aL*FR(bt1) + f2*(VR(i,j,k,bt1)-VL(i,j,k,bt1)) )
     F(i,j,k,bt2)   = f1*( aR*FL(bt2) - aL*FR(bt2) + f2*(VR(i,j,k,bt2)-VL(i,j,k,bt2)) )

!    entropy wave
     f2  = 1 / U_hll(ro)
     aM  = U_hll(mn ) * f2 ! vx_HLL
     at1 = U_hll(mt1) * f2 ! vy_HLL
     at2 = U_hll(mt2) * f2 ! vz_HLL
!    Total pressure (these two are identical)
     pt = ptL + VL(i,j,k,ro) * ( aL - VL(i,j,k,vn) ) * ( aM - VL(i,j,k,vn) )
!    pt = ptR + VR(i,j,k,ro) * ( aR - VR(i,j,k,vn) ) * ( aM - VR(i,j,k,vn) )

!!    F = F(L*) or F(L)
!     if ( aM >= 0 ) then

     roL = VL(i,j,k,ro) * ( aL - VL(i,j,k,vn) ) / ( aL - aM )
     enL = ( ( aL - VL(i,j,k,vn) )*UL(en) - ptL*VL(i,j,k,vn) + pt * aM + &
          U_hll(bn)*( vBL - aM*U_hll(bn) - at1*U_hll(bt1) - at2*U_hll(bt2) ) ) / ( aL - aM )
!     F(i,j,k,mn) = aL *( roL*aM - UL(mn) ) + FL(mn)
!     F(i,j,k,mt1)= aL *( roL*at1- UL(mt1)) + FL(mt1)
!     F(i,j,k,mt2)= aL *( roL*at2- UL(mt2)) + FL(mt2)
!     F(i,j,k,en) = aL *( enL    - UL(en) ) + FL(en)
!     F(i,j,k,ro) = aL *( roL    - VL(i,j,k,ro) ) + FL(ro)

!!    F = F(R*) or F(R)
!     else

     roR = VR(i,j,k,ro) * ( aR - VR(i,j,k,vn) ) / ( aR - aM )
     enR = ( ( aR - VR(i,j,k,vn) )*UR(en) - ptR*VR(i,j,k,vn) + pt * aM + &
          U_hll(bn)*( vBR - aM*U_hll(bn) - at1*U_hll(bt1) - at2*U_hll(bt2) ) ) / ( aR - aM )
!     F(i,j,k,mn) = aR *( roR*aM - UR(mn) ) + FR(mn)
!     F(i,j,k,mt1)= aR *( roR*at1- UR(mt1)) + FR(mt1)
!     F(i,j,k,mt2)= aR *( roR*at2- UR(mt2)) + FR(mt2)
!     F(i,j,k,en) = aR *( enR    - UR(en) ) + FR(en)
!     F(i,j,k,ro) = aR *( roR    - VR(i,j,k,ro) ) + FR(ro)

!    Weight factor, 0 or 1.  This looks tricky.
!    The code runs 1.0x times slower on Intel, but it runs 1.6 times faster on SPARC.
     f1 = 0.5d0 + sign(0.5d0,aM)  !!  F = F(L*) or F(L)
     f2 = 1.d0 - f1               !!  F = F(R*) or F(R)

     F(i,j,k,mn)  = f1*( aL *( roL*aM - UL(mn) ) + FL(mn) ) &
                  + f2*( aR *( roR*aM - UR(mn) ) + FR(mn) )
     F(i,j,k,mt1) = f1*( aL *( roL*at1- UL(mt1)) + FL(mt1)) &
                  + f2*( aR *( roR*at1- UR(mt1)) + FR(mt1))
     F(i,j,k,mt2) = f1*( aL *( roL*at2- UL(mt2)) + FL(mt2)) &
                  + f2*( aR *( roR*at2- UR(mt2)) + FR(mt2))
     F(i,j,k,en)  = f1*( aL *( enL    - UL(en) ) + FL(en) ) &
                  + f2*( aR *( enR    - UR(en) ) + FR(en) )
     F(i,j,k,ro)  = f1*( aL *( roL    - VL(i,j,k,ro) ) + FL(ro) ) &
                  + f2*( aR *( roR    - VR(i,j,k,ro) ) + FR(ro) )

     F(i,j,k,bn) = 0.5d0*(     ( VL(i,j,k,ps)+VR(i,j,k,ps) ) - ch*( VR(i,j,k,bn)-VL(i,j,k,bn) ) )
     F(i,j,k,ps) = 0.5d0*( ch2*( VL(i,j,k,bn)+VR(i,j,k,bn) ) - ch*( VR(i,j,k,ps)-VL(i,j,k,ps) ) )

  enddo
  enddo
  enddo
!$omp end parallel do

  !-----------------------------------------------------------------------
  !   HLLD solver
  !-----------------------------------------------------------------------
  case(3)
!$omp parallel do private(i,j,k,v2,B2,vBL,vBR,UL,UR,FL,FR,f1,f2) &
!$omp private(aL,aR,aM,at1,at2,ptL,ptR,pt,roL,roR,enL,enR,U_hll) &
!$omp private(aL1,aR1,aVL,aVR,UL1,UR1,U2,vt1L,vt1R,vt2L,vt2R,roLs,roRs)
  do k=ks,ke
  do j=js,je
  do i=is,ie

!    VL --> UL
     v2 = dot_product( VL(i,j,k,vx:vz), VL(i,j,k,vx:vz) )
     B2 = dot_product( VL(i,j,k,bx:bz), VL(i,j,k,bx:bz) )
     vBL= dot_product( VL(i,j,k,vx:vz), VL(i,j,k,bx:bz) )
     UL(mx:mz) = VL(i,j,k,ro) * VL(i,j,k,vx:vz)
     UL(en)    = 0.5d0 * ( VL(i,j,k,ro)*v2 + B2 ) + r1*VL(i,j,k,pr)
!    VL --> FL
     FL(ro)  = VL(i,j,k,ro) * VL(i,j,k,vn)
     FL(mn)  = VL(i,j,k,ro)*VL(i,j,k,vn)*VL(i,j,k,vn)  - VL(i,j,k,bn)*VL(i,j,k,bn) + VL(i,j,k,pr) + 0.5d0*B2
     FL(mt1) = VL(i,j,k,ro)*VL(i,j,k,vn)*VL(i,j,k,vt1) - VL(i,j,k,bn)*VL(i,j,k,bt1)
     FL(mt2) = VL(i,j,k,ro)*VL(i,j,k,vn)*VL(i,j,k,vt2) - VL(i,j,k,bn)*VL(i,j,k,bt2)
     FL(en)  = ( 0.5d0*VL(i,j,k,ro)*v2 + r2*VL(i,j,k,pr) + B2 ) * VL(i,j,k,vn) - vBL * VL(i,j,k,bn)
     FL(bt1) = VL(i,j,k,vn)*VL(i,j,k,bt1)-VL(i,j,k,bn)*VL(i,j,k,vt1)
     FL(bt2) = VL(i,j,k,vn)*VL(i,j,k,bt2)-VL(i,j,k,bn)*VL(i,j,k,vt2)

!    VL --> fast mode^2, total pressure
     f1 = gamma * VL(i,j,k,pr) ! f1: gamma p
     f2 = 4 * VL(i,j,k,bn)**2  ! f2: 4 B_n^2
!    vfL = sqrt( ( (f1+B2) + sqrt( (f1+B2)**2 - f1*f2 )) / ( 2*VL(i,j,k,ro) ))
     aL  = ( (f1+B2) + sqrt(max( (f1+B2)**2 - f1*f2, 0.d0 ))) / ( 2*VL(i,j,k,ro) )
     ptL = VL(i,j,k,pr) + 0.5d0*B2

!    VR --> UR
     v2 = dot_product( VR(i,j,k,vx:vz), VR(i,j,k,vx:vz) )
     B2 = dot_product( VR(i,j,k,bx:bz), VR(i,j,k,bx:bz) )
     vBR= dot_product( VR(i,j,k,vx:vz), VR(i,j,k,bx:bz) )
     UR(mx:mz) = VR(i,j,k,ro) * VR(i,j,k,vx:vz)
     UR(en)    = 0.5d0 * ( VR(i,j,k,ro)*v2 + B2 ) + r1*VR(i,j,k,pr)
!    VR --> FR
     FR(ro)  = VR(i,j,k,ro) * VR(i,j,k,vn)
     FR(mn)  = VR(i,j,k,ro)*VR(i,j,k,vn)*VR(i,j,k,vn)  - VR(i,j,k,bn)*VR(i,j,k,bn) + VR(i,j,k,pr) + 0.5d0*B2
     FR(mt1) = VR(i,j,k,ro)*VR(i,j,k,vn)*VR(i,j,k,vt1) - VR(i,j,k,bn)*VR(i,j,k,bt1)
     FR(mt2) = VR(i,j,k,ro)*VR(i,j,k,vn)*VR(i,j,k,vt2) - VR(i,j,k,bn)*VR(i,j,k,bt2)
     FR(en)  = ( 0.5d0*VR(i,j,k,ro)*v2 + r2*VR(i,j,k,pr) + B2 ) * VR(i,j,k,vn) - vBR * VR(i,j,k,bn)
     FR(bt1) = VR(i,j,k,vn)*VR(i,j,k,bt1)-VR(i,j,k,bn)*VR(i,j,k,vt1)
     FR(bt2) = VR(i,j,k,vn)*VR(i,j,k,bt2)-VR(i,j,k,bn)*VR(i,j,k,vt2)

!    VR --> fast mode^2, total pressure
     f1 = gamma * VR(i,j,k,pr) ! f1: gamma p
     f2 = 4 * VR(i,j,k,bn)**2  ! f2: 4 B_n^2
!    vfR = sqrt( ( (f1+B2) + sqrt( (f1+B2)**2 - f1*f2 )) / ( 2*VR(i,j,k,ro) ))
     aR  = ( (f1+B2) + sqrt(max( (f1+B2)**2 - f1*f2, 0.d0 ))) / ( 2*VR(i,j,k,ro) )
     ptR = VR(i,j,k,pr) + 0.5d0*B2

!    Riemann fan speed (MK05 eq. 67)
!     aL = min( VL(i,j,k,vn) - vfL, VR(i,j,k,vn) - vfR )
!     aR = max( VL(i,j,k,vn) + vfL, VR(i,j,k,vn) + vfR )
!     aL = min( VL(i,j,k,vn), VR(i,j,k,vn) ) - max( vfL, vfR )
!     aR = max( VL(i,j,k,vn), VR(i,j,k,vn) ) + max( vfL, vfR )
     f1 = sqrt( max( aL, aR ) )  ! faster fast-wave (This aL [aR] stores vfL^2 [vfR^2])
     aL = min( VL(i,j,k,vn), VR(i,j,k,vn) ) - f1
     aR = max( VL(i,j,k,vn), VR(i,j,k,vn) ) + f1

!    F = F(L)
     if ( aL >= 0 ) then
        F(i,j,k,:) = FL(:)
!    F = F(R)
     elseif ( aR <= 0 ) then
        F(i,j,k,:) = FR(:)

!    HLLC/HLLD flux
     else

!       HLL state
        f1 = 1.d0 / ( aR - aL )
        U_hll(mx:mz) = f1*( aR*UR(mx:mz) - aL*UL(mx:mz) - FR(mx:mz) + FL(mx:mz) )
        U_hll(ro:bz) = f1*( aR*VR(i,j,k,ro:bz) - aL*VL(i,j,k,ro:bz) - FR(ro:bz) + FL(ro:bz) )

!       entropy wave
        aM = U_hll(mn) / U_hll(ro)

!       Total pressure
        pt  = ptL + VL(i,j,k,ro) * ( aL - VL(i,j,k,vn) ) * ( aM - VL(i,j,k,vn) )
        roL = VL(i,j,k,ro) * ( aL - VL(i,j,k,vn) ) / ( aL - aM ) ! ro(L*)
        roR = VR(i,j,k,ro) * ( aR - VR(i,j,k,vn) ) / ( aR - aM ) ! ro(R*)

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

           !f1 = 1.d0 / ( aR - aL )
           f2 = aL * aR
           F(i,j,k,bt1)   = f1*( aR*FL(bt1) - aL*FR(bt1) + f2*(VR(i,j,k,bt1)-VL(i,j,k,bt1)) )
           F(i,j,k,bt2)   = f1*( aR*FL(bt2) - aL*FR(bt2) + f2*(VR(i,j,k,bt2)-VL(i,j,k,bt2)) )

           at1 = U_hll(mt1) / U_hll(ro) ! vt1_HLL
           at2 = U_hll(mt2) / U_hll(ro) ! vt2_HLL

           enL = ( ( aL-VL(i,j,k,vn) )*UL(en) - ptL*VL(i,j,k,vn) + pt*aM + &
                U_hll(bn)*( vBL - aM*U_hll(bn) - at1*U_hll(bt1) - at2*U_hll(bt2) ) ) / ( aL - aM )
           enR = ( ( aR-VR(i,j,k,vn) )*UR(en) - ptR*VR(i,j,k,vn) + pt*aM + &
                U_hll(bn)*( vBR - aM*U_hll(bn) - at1*U_hll(bt1) - at2*U_hll(bt2) ) ) / ( aR - aM )

!          Weight factor, 0 or 1.  This looks tricky.
           f1 = 0.5d0 + sign(0.5d0,aM)  !!  F = F(L*)
           f2 = 1.d0 - f1               !!  F = F(R*)

           F(i,j,k,mn)  = f1*( aL *( roL*aM - UL(mn) ) + FL(mn) ) &
                        + f2*( aR *( roR*aM - UR(mn) ) + FR(mn) )
           F(i,j,k,mt1) = f1*( aL *( roL*at1- UL(mt1)) + FL(mt1)) &
                        + f2*( aR *( roR*at1- UR(mt1)) + FR(mt1))
           F(i,j,k,mt2) = f1*( aL *( roL*at2- UL(mt2)) + FL(mt2)) &
                        + f2*( aR *( roR*at2- UR(mt2)) + FR(mt2))
           F(i,j,k,en)  = f1*( aL *( enL    - UL(en) ) + FL(en) ) &
                        + f2*( aR *( enR    - UR(en) ) + FR(en) )
           F(i,j,k,ro)  = f1*( aL *( roL    - VL(i,j,k,ro) ) + FL(ro) ) &
                        + f2*( aR *( roR    - VR(i,j,k,ro) ) + FR(ro) )

! ========== HLLD flux ==========
        else

!          Intermediate state (L*)
           f1 = 1.d0 / ( VL(i,j,k,ro)*(aL-VL(i,j,k,vn))*(aL-aM) - U_hll(bn)**2 ) ! HLLD denominator
           UL1(bt1) = VL(i,j,k,bt1) * f1 * ( VL(i,j,k,ro)*(aL-VL(i,j,k,vn))**2 - U_hll(bn)**2 )
           UL1(bt2) = VL(i,j,k,bt2) * f1 * ( VL(i,j,k,ro)*(aL-VL(i,j,k,vn))**2 - U_hll(bn)**2 )
           vt1L     = VL(i,j,k,vt1) - f1 * U_hll(bn)*VL(i,j,k,bt1)*(aM-VL(i,j,k,vn))
           vt2L     = VL(i,j,k,vt2) - f1 * U_hll(bn)*VL(i,j,k,bt2)*(aM-VL(i,j,k,vn))

           UL1(ro)  = roL
           UL1(en)  = ( ( aL-VL(i,j,k,vn) )*UL(en) - ptL*VL(i,j,k,vn) + pt*aM + &
                U_hll(bn)*( vBL - aM*U_hll(bn) - vt1L*UL1(bt1) - vt2L*UL1(bt2)) ) / ( aL - aM )
           UL1(mn)  = roL * aM
           UL1(mt1) = roL * vt1L
           UL1(mt2) = roL * vt2L

!          Intermediate state (R*)
           f1 = 1.d0 / ( VR(i,j,k,ro)*(aR-VR(i,j,k,vn))*(aR-aM) - U_hll(bn)**2 ) ! HLLD denominator
           UR1(bt1) = VR(i,j,k,bt1) * f1 * ( VR(i,j,k,ro)*(aR-VR(i,j,k,vn))**2 - U_hll(bn)**2 )
           UR1(bt2) = VR(i,j,k,bt2) * f1 * ( VR(i,j,k,ro)*(aR-VR(i,j,k,vn))**2 - U_hll(bn)**2 )
           vt1R     = VR(i,j,k,vt1) - f1 * U_hll(bn)*VR(i,j,k,bt1)*(aM-VR(i,j,k,vn))
           vt2R     = VR(i,j,k,vt2) - f1 * U_hll(bn)*VR(i,j,k,bt2)*(aM-VR(i,j,k,vn))

           UR1(ro)  = roR
           UR1(en)  = ( ( aR-VR(i,j,k,vn) )*UR(en) - ptR*VR(i,j,k,vn) + pt*aM + &
                U_hll(bn)*( vBR - aM*U_hll(bn) - vt1R*UR1(bt1) - vt2R*UR1(bt2)) ) / ( aR - aM )
           UR1(mn)  = roR * aM
           UR1(mt1) = roR * vt1R
           UR1(mt2) = roR * vt2R

!          rotational waves
           aL1 = aM - aVL
           aR1 = aM + aVR

!          F = F(L*)
           if( aL1 >= 0 ) then

              F(i,j,k,mx:en) = aL *( UL1(mx:en) - UL(mx:en)   ) + FL(mx:en)
              F(i,j,k,ro)    = aL *( UL1(ro)    - VL(i,j,k,ro)  ) + FL(ro)
              F(i,j,k,bt1)   = aL *( UL1(bt1)   - VL(i,j,k,bt1) ) + FL(bt1)
              F(i,j,k,bt2)   = aL *( UL1(bt2)   - VL(i,j,k,bt2) ) + FL(bt2)

!          F = F(R*)
           elseif( aR1 <= 0 ) then

              F(i,j,k,mx:en) = aR *( UR1(mx:en) - UR(mx:en)   ) + FR(mx:en)
              F(i,j,k,ro)    = aR *( UR1(ro)    - VR(i,j,k,ro)  ) + FR(ro)
              F(i,j,k,bt1)   = aR *( UR1(bt1)   - VR(i,j,k,bt1) ) + FR(bt1)
              F(i,j,k,bt2)   = aR *( UR1(bt2)   - VR(i,j,k,bt2) ) + FR(bt2)

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

                 F(i,j,k,mx:en) = aL1*U2(mx:en) - (aL1-aL)*UL1(mx:en) - aL*UL(mx:en) + FL(mx:en)
                 F(i,j,k,ro)    = aL1*U2(ro)  - (aL1-aL)*UL1(ro)  - aL*VL(i,j,k,ro)  + FL(ro)
                 F(i,j,k,bt1)   = aL1*U2(bt1) - (aL1-aL)*UL1(bt1) - aL*VL(i,j,k,bt1) + FL(bt1)
                 F(i,j,k,bt2)   = aL1*U2(bt2) - (aL1-aL)*UL1(bt2) - aL*VL(i,j,k,bt2) + FL(bt2)

!             F = F(R**)
              else

                 U2(ro) = roR
                 U2(en) = UR1(en) + roRs * ( vt1R*UR1(bt1) + vt2R*UR1(bt2) &
                      - at1*U2(bt1) - at2*U2(bt2) ) * f2
                 U2(mn)  = roR * aM
                 U2(mt1) = roR * at1
                 U2(mt2) = roR * at2

                 F(i,j,k,mx:en) = aR1*U2(mx:en) - (aR1-aR)*UR1(mx:en) - aR*UR(mx:en) + FR(mx:en)
                 F(i,j,k,ro)    = aR1*U2(ro)  - (aR1-aR)*UR1(ro)  - aR*VR(i,j,k,ro)  + FR(ro)
                 F(i,j,k,bt1)   = aR1*U2(bt1) - (aR1-aR)*UR1(bt1) - aR*VR(i,j,k,bt1) + FR(bt1)
                 F(i,j,k,bt2)   = aR1*U2(bt2) - (aR1-aR)*UR1(bt2) - aR*VR(i,j,k,bt2) + FR(bt2)

              endif

           endif

        endif

     endif

     F(i,j,k,bn) = 0.5d0*(     ( VL(i,j,k,ps)+VR(i,j,k,ps) ) - ch*( VR(i,j,k,bn)-VL(i,j,k,bn) ) )
     F(i,j,k,ps) = 0.5d0*( ch2*( VL(i,j,k,bn)+VR(i,j,k,bn) ) - ch*( VR(i,j,k,ps)-VL(i,j,k,ps) ) )

  enddo
  enddo
  enddo
!$omp end parallel do

  !-----------------------------------------------------------------------
  !   LHLLD solver (experimental)
  !-----------------------------------------------------------------------
  case(4)
!$omp parallel do private(i,j,k,v2,B2,vBL,vBR,UL,UR,FL,FR,f1,f2) &
!$omp private(aL,aR,aM,at1,at2,ptL,ptR,pt,roL,roR,enL,enR,U_hll) &
!$omp private(aL1,aR1,aVL,aVR,UL1,UR1,U2,vt1L,vt1R,vt2L,vt2R,roLs,roRs)
  do k=ks,ke
  do j=js,je
  do i=is,ie

!    VL --> UL
     v2 = dot_product( VL(i,j,k,vx:vz), VL(i,j,k,vx:vz) )
     B2 = dot_product( VL(i,j,k,bx:bz), VL(i,j,k,bx:bz) )
     vBL= dot_product( VL(i,j,k,vx:vz), VL(i,j,k,bx:bz) )
     UL(mx:mz) = VL(i,j,k,ro) * VL(i,j,k,vx:vz)
     UL(en)    = 0.5d0 * ( VL(i,j,k,ro)*v2 + B2 ) + r1*VL(i,j,k,pr)
!    VL --> FL
     FL(ro)  = VL(i,j,k,ro) * VL(i,j,k,vn)
     FL(mn)  = VL(i,j,k,ro)*VL(i,j,k,vn)*VL(i,j,k,vn)  - VL(i,j,k,bn)*VL(i,j,k,bn) + VL(i,j,k,pr) + 0.5d0*B2
     FL(mt1) = VL(i,j,k,ro)*VL(i,j,k,vn)*VL(i,j,k,vt1) - VL(i,j,k,bn)*VL(i,j,k,bt1)
     FL(mt2) = VL(i,j,k,ro)*VL(i,j,k,vn)*VL(i,j,k,vt2) - VL(i,j,k,bn)*VL(i,j,k,bt2)
     FL(en)  = ( 0.5d0*VL(i,j,k,ro)*v2 + r2*VL(i,j,k,pr) + B2 ) * VL(i,j,k,vn) - vBL * VL(i,j,k,bn)
     FL(bt1) = VL(i,j,k,vn)*VL(i,j,k,bt1)-VL(i,j,k,bn)*VL(i,j,k,vt1)
     FL(bt2) = VL(i,j,k,vn)*VL(i,j,k,bt2)-VL(i,j,k,bn)*VL(i,j,k,vt2)

!    VL --> fast mode^2, total pressure
     f1 = gamma * VL(i,j,k,pr) ! f1: gamma p
     f2 = 4 * VL(i,j,k,bn)**2  ! f2: 4 B_n^2
!    vfL = sqrt( ( (f1+B2) + sqrt( (f1+B2)**2 - f1*f2 )) / ( 2*VL(i,j,k,ro) ))
     aL  = ( (f1+B2) + sqrt(max( (f1+B2)**2 - f1*f2, 0.d0 ))) / ( 2*VL(i,j,k,ro) )
     ptL = VL(i,j,k,pr) + 0.5d0*B2
!    LHLLD
     f1  = VL(i,j,k,ro) * v2
     aL1 = ( (f1+B2) + sqrt(max( (f1+B2)**2 - f1*f2, 0.d0 ))) / ( 2*VL(i,j,k,ro) )

!    VR --> UR
     v2 = dot_product( VR(i,j,k,vx:vz), VR(i,j,k,vx:vz) )
     B2 = dot_product( VR(i,j,k,bx:bz), VR(i,j,k,bx:bz) )
     vBR= dot_product( VR(i,j,k,vx:vz), VR(i,j,k,bx:bz) )
     UR(mx:mz) = VR(i,j,k,ro) * VR(i,j,k,vx:vz)
     UR(en)    = 0.5d0 * ( VR(i,j,k,ro)*v2 + B2 ) + r1*VR(i,j,k,pr)
!    VR --> FR
     FR(ro)  = VR(i,j,k,ro) * VR(i,j,k,vn)
     FR(mn)  = VR(i,j,k,ro)*VR(i,j,k,vn)*VR(i,j,k,vn)  - VR(i,j,k,bn)*VR(i,j,k,bn) + VR(i,j,k,pr) + 0.5d0*B2
     FR(mt1) = VR(i,j,k,ro)*VR(i,j,k,vn)*VR(i,j,k,vt1) - VR(i,j,k,bn)*VR(i,j,k,bt1)
     FR(mt2) = VR(i,j,k,ro)*VR(i,j,k,vn)*VR(i,j,k,vt2) - VR(i,j,k,bn)*VR(i,j,k,bt2)
     FR(en)  = ( 0.5d0*VR(i,j,k,ro)*v2 + r2*VR(i,j,k,pr) + B2 ) * VR(i,j,k,vn) - vBR * VR(i,j,k,bn)
     FR(bt1) = VR(i,j,k,vn)*VR(i,j,k,bt1)-VR(i,j,k,bn)*VR(i,j,k,vt1)
     FR(bt2) = VR(i,j,k,vn)*VR(i,j,k,bt2)-VR(i,j,k,bn)*VR(i,j,k,vt2)

!    VR --> fast mode^2, total pressure
     f1 = gamma * VR(i,j,k,pr) ! f1: gamma p
     f2 = 4 * VR(i,j,k,bn)**2  ! f2: 4 B_n^2
!    vfR = sqrt( ( (f1+B2) + sqrt( (f1+B2)**2 - f1*f2 )) / ( 2*VR(i,j,k,ro) ))
     aR  = ( (f1+B2) + sqrt(max( (f1+B2)**2 - f1*f2, 0.d0 ))) / ( 2*VR(i,j,k,ro) )
     ptR = VR(i,j,k,pr) + 0.5d0*B2
!    LHLLD
     f1  = VR(i,j,k,ro) * v2
     aR1 = ( (f1+B2) + sqrt(max( (f1+B2)**2 - f1*f2, 0.d0 ))) / ( 2*VR(i,j,k,ro) )

!    pressure correction factor psi (f2)
     f1 = min( 1.d0, max( aL1, aR1 ) / max( aL, aR ) )  ! chi**2
     f2 = 2 * sqrt(f1) - f1  ! psi

!    Riemann fan speed (MK05 eq. 67)
     f1 = sqrt( max( aL, aR ) )  ! faster fast-wave (This aL [aR] stores vfL^2 [vfR^2])
     aL = min( VL(i,j,k,vn), VR(i,j,k,vn) ) - f1
     aR = max( VL(i,j,k,vn), VR(i,j,k,vn) ) + f1

!    F = F(L)
     if ( aL >= 0 ) then
        F(i,j,k,:) = FL(:)
!    F = F(R)
     elseif ( aR <= 0 ) then
        F(i,j,k,:) = FR(:)

!    HLLC/HLLD flux
     else

!       HLL state
        f1 = 1.d0 / ( aR - aL )
        U_hll(mx:mz) = f1*( aR*UR(mx:mz) - aL*UL(mx:mz) - FR(mx:mz) + FL(mx:mz) )
        U_hll(ro:bz) = f1*( aR*VR(i,j,k,ro:bz) - aL*VL(i,j,k,ro:bz) - FR(ro:bz) + FL(ro:bz) )

!       entropy wave
        aM = U_hll(mn) / U_hll(ro)

!       Total pressure
!       pt  = ptL + VL(i,j,k,ro) * ( aL - VL(i,j,k,vn) ) * ( aM - VL(i,j,k,vn) )
!       f2 (psi) is defined 25 lines earlier
        pt  = ( VR(i,j,k,ro)*ptL*(aR-VR(i,j,k,vn)) - VL(i,j,k,ro)*ptR*(aL-VL(i,j,k,vn)) &
             + f2*VL(i,j,k,ro)*VR(i,j,k,ro)*(aR-VR(i,j,k,vn))*(aL-VL(i,j,k,vn))*(VR(i,j,k,vn)-VL(i,j,k,vn)) ) &
             / ( VR(i,j,k,ro)*(aR-VR(i,j,k,vn)) - VL(i,j,k,ro)*(aL-VL(i,j,k,vn)) )
        roL = VL(i,j,k,ro) * ( aL - VL(i,j,k,vn) ) / ( aL - aM ) ! ro(L*)
        roR = VR(i,j,k,ro) * ( aR - VR(i,j,k,vn) ) / ( aR - aM ) ! ro(R*)

!       For logical consistency, we employ Bn_hll as B_x in the intermediate states
        aVL = abs( U_hll(bn) )/sqrt( roL )  ! Alfven wave (L)
        aVR = abs( U_hll(bn) )/sqrt( roR )  ! Alfven wave (R)

!       Intermediate state (L*)
        f1 = 1.d0 / ( VL(i,j,k,ro)*(aL-VL(i,j,k,vn))*(aL-aM) - U_hll(bn)**2 ) ! HLLD denominator
        UL1(bt1) = VL(i,j,k,bt1) * f1 * ( VL(i,j,k,ro)*(aL-VL(i,j,k,vn))**2 - U_hll(bn)**2 )
        UL1(bt2) = VL(i,j,k,bt2) * f1 * ( VL(i,j,k,ro)*(aL-VL(i,j,k,vn))**2 - U_hll(bn)**2 )
        vt1L     = VL(i,j,k,vt1) - f1 * U_hll(bn)*VL(i,j,k,bt1)*(aM-VL(i,j,k,vn))
        vt2L     = VL(i,j,k,vt2) - f1 * U_hll(bn)*VL(i,j,k,bt2)*(aM-VL(i,j,k,vn))

        UL1(ro)  = roL
        UL1(en)  = ( ( aL-VL(i,j,k,vn) )*UL(en) - ptL*VL(i,j,k,vn) + pt*aM + &
             U_hll(bn)*( vBL - aM*U_hll(bn) - vt1L*UL1(bt1) - vt2L*UL1(bt2)) ) / ( aL - aM )
        UL1(mn)  = roL * aM
        UL1(mt1) = roL * vt1L
        UL1(mt2) = roL * vt2L

!       Intermediate state (R*)
        f1 = 1.d0 / ( VR(i,j,k,ro)*(aR-VR(i,j,k,vn))*(aR-aM) - U_hll(bn)**2 ) ! HLLD denominator
        UR1(bt1) = VR(i,j,k,bt1) * f1 * ( VR(i,j,k,ro)*(aR-VR(i,j,k,vn))**2 - U_hll(bn)**2 )
        UR1(bt2) = VR(i,j,k,bt2) * f1 * ( VR(i,j,k,ro)*(aR-VR(i,j,k,vn))**2 - U_hll(bn)**2 )
        vt1R     = VR(i,j,k,vt1) - f1 * U_hll(bn)*VR(i,j,k,bt1)*(aM-VR(i,j,k,vn))
        vt2R     = VR(i,j,k,vt2) - f1 * U_hll(bn)*VR(i,j,k,bt2)*(aM-VR(i,j,k,vn))

        UR1(ro)  = roR
        UR1(en)  = ( ( aR-VR(i,j,k,vn) )*UR(en) - ptR*VR(i,j,k,vn) + pt*aM + &
             U_hll(bn)*( vBR - aM*U_hll(bn) - vt1R*UR1(bt1) - vt2R*UR1(bt2)) ) / ( aR - aM )
        UR1(mn)  = roR * aM
        UR1(mt1) = roR * vt1R
        UR1(mt2) = roR * vt2R

!       rotational waves
        aL1 = aM - aVL
        aR1 = aM + aVR

!       F = F(L*)
        if( aL1 >= 0 ) then

           F(i,j,k,mx:en) = aL *( UL1(mx:en) - UL(mx:en)   ) + FL(mx:en)
           F(i,j,k,ro)    = aL *( UL1(ro)    - VL(i,j,k,ro)  ) + FL(ro)
           F(i,j,k,bt1)   = aL *( UL1(bt1)   - VL(i,j,k,bt1) ) + FL(bt1)
           F(i,j,k,bt2)   = aL *( UL1(bt2)   - VL(i,j,k,bt2) ) + FL(bt2)

!       F = F(R*)
        elseif( aR1 <= 0 ) then

           F(i,j,k,mx:en) = aR *( UR1(mx:en) - UR(mx:en)   ) + FR(mx:en)
           F(i,j,k,ro)    = aR *( UR1(ro)    - VR(i,j,k,ro)  ) + FR(ro)
           F(i,j,k,bt1)   = aR *( UR1(bt1)   - VR(i,j,k,bt1) ) + FR(bt1)
           F(i,j,k,bt2)   = aR *( UR1(bt2)   - VR(i,j,k,bt2) ) + FR(bt2)

!       Central states are tricky
!       Question: can we really rule out U_hll(bn) == 0 ?
        else

           roLs = sqrt( roL )  ! sqrt(ro(L*))
           roRs = sqrt( roR )  ! sqrt(ro(R*))
           f1   = 1.d0 / ( roLs + roRs )
           f2   = sign( 1.d0, U_hll(bn) )
           at1  = f1 * ( roLs*vt1L + roRs*vt1R + ( UR1(bt1)-UL1(bt1) )*f2 )
           at2  = f1 * ( roLs*vt2L + roRs*vt2R + ( UR1(bt2)-UL1(bt2) )*f2 )
           U2(bt1) = f1 * ( roLs*UR1(bt1) + roRs*UL1(bt1) + roLs*roRs*( vt1R-vt1L )*f2 )
           U2(bt2) = f1 * ( roLs*UR1(bt2) + roRs*UL1(bt2) + roLs*roRs*( vt2R-vt2L )*f2 )

!          F = F(L**)
           if( aM >= 0 ) then

              U2(ro) = roL
              U2(en) = UL1(en) - roLs * ( vt1L*UL1(bt1) + vt2L*UL1(bt2) &
                   - at1*U2(bt1) - at2*U2(bt2) ) * f2
              U2(mn)  = roL * aM
              U2(mt1) = roL * at1
              U2(mt2) = roL * at2

              F(i,j,k,mx:en) = aL1*U2(mx:en) - (aL1-aL)*UL1(mx:en) - aL*UL(mx:en) + FL(mx:en)
              F(i,j,k,ro)    = aL1*U2(ro)  - (aL1-aL)*UL1(ro)  - aL*VL(i,j,k,ro)  + FL(ro)
              F(i,j,k,bt1)   = aL1*U2(bt1) - (aL1-aL)*UL1(bt1) - aL*VL(i,j,k,bt1) + FL(bt1)
              F(i,j,k,bt2)   = aL1*U2(bt2) - (aL1-aL)*UL1(bt2) - aL*VL(i,j,k,bt2) + FL(bt2)

!          F = F(R**)
           else

              U2(ro) = roR
              U2(en) = UR1(en) + roRs * ( vt1R*UR1(bt1) + vt2R*UR1(bt2) &
                   - at1*U2(bt1) - at2*U2(bt2) ) * f2
              U2(mn)  = roR * aM
              U2(mt1) = roR * at1
              U2(mt2) = roR * at2

              F(i,j,k,mx:en) = aR1*U2(mx:en) - (aR1-aR)*UR1(mx:en) - aR*UR(mx:en) + FR(mx:en)
              F(i,j,k,ro)    = aR1*U2(ro)  - (aR1-aR)*UR1(ro)  - aR*VR(i,j,k,ro)  + FR(ro)
              F(i,j,k,bt1)   = aR1*U2(bt1) - (aR1-aR)*UR1(bt1) - aR*VR(i,j,k,bt1) + FR(bt1)
              F(i,j,k,bt2)   = aR1*U2(bt2) - (aR1-aR)*UR1(bt2) - aR*VR(i,j,k,bt2) + FR(bt2)

           endif

        endif

     endif

     F(i,j,k,bn) = 0.5d0*(     ( VL(i,j,k,ps)+VR(i,j,k,ps) ) - ch*( VR(i,j,k,bn)-VL(i,j,k,bn) ) )
     F(i,j,k,ps) = 0.5d0*( ch2*( VL(i,j,k,bn)+VR(i,j,k,bn) ) - ch*( VR(i,j,k,ps)-VL(i,j,k,ps) ) )

  enddo
  enddo
  enddo
!$omp end parallel do

  end select
  !-----------------------------------------------------------------------

  return
end subroutine flux_solver
