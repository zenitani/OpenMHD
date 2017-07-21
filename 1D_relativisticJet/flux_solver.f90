subroutine flux_solver(F,VL,VR,ix,jx,dir,hlld)
!-----------------------------------------------------------------------
!     relativistic HLL/HLLD solver
!     Ref. Mignone & Bodo (2006)
!     Ref. Mignone et al. (2009) section 3.4.1
!-----------------------------------------------------------------------
!     2009/08/01  S. Zenitani  standard HLL solver
!     2009/08/05  S. Zenitani   + simplified HLLD solver
!     2010/04/26  S. Zenitani  a bug fixed (sound speed in HLL flux)
!     2017/07/20  S. Zenitani  u-based solver
!-----------------------------------------------------------------------
  implicit none
  include 'param_rela.h'
!-----------------------------------------------------------------------
! size of arrays [input]
  integer, intent(in)  :: ix, jx
! numerical flux (F) [output]
  real(8), intent(out) :: F(ix,jx,var1)
! left/right states (VL,VR) [input]
  real(8), intent(in)  :: VL(ix,jx,var1), VR(ix,jx,var1)
! direction [input]: 1 (X), 2 (Y), 3 (Z)
  integer, intent(in)  :: dir
! left/right conserved variables (UL & UR; local)
  real(8) :: UL(var1), UR(var1)
! left/right numerical flux (FL & FR; local)
  real(8) :: FL(var1), FR(var1)
! .true. ==> HLLD, .false. ==> HLL
  logical :: hlld

  real(8), parameter :: threshold = 1.0d-5
! numerical flux (F) [output] -- to be removed
  real(8) :: Fd(ix), Fm(ix,3), Fe(ix), Fb(ix,3)
  integer :: i, j

  real(8) :: v2, u0, u02, u2, W, B2, vB, uB
  real(8) :: aL, aR, aLR, a
  real(8) :: f1L, f1R, f3L, f3R
  real(8) :: a0L, a0R, a1L, a1R, a2L, a2R
! for hlld states
  real(8) :: pt, ehll, mxhll, a1, a0
  integer :: un, ut1, ut2, bn, bt1, bt2, mn, mt1, mt2
  real(8) :: aM, at1, at2, mybn, mybt1, mybt2, f1, f2
  real(8) :: Rd, Re, Rmn, Rmt1, Rmt2, Rbt1, Rbt2


  F = 0.d0
  FL = 0.d0; FR = 0.d0

! directions
  un  = ux + mod(dir-1,3)
  ut1 = ux + mod(dir  ,3)
  ut2 = ux + mod(dir+1,3)
  bn  = bx + mod(dir-1,3)
  bt1 = bx + mod(dir  ,3)
  bt2 = bx + mod(dir+1,3)
  mn  = mx + mod(dir-1,3)
  mt1 = mx + mod(dir  ,3)
  mt2 = mx + mod(dir+1,3)

  j=1
  do i=1,ix-1

!! ========== This requires B=0 ========== !!

!    VL --> UL
     B2 = dot_product( VL(i,j,bx:bz), VL(i,j,bx:bz) )
     uB = dot_product( VL(i,j,bx:bz), VL(i,j,ux:uz) )
     u2 = dot_product( VL(i,j,ux:uz), VL(i,j,ux:uz) )
     u02= 1.d0 + u2
     u0 = sqrt( u02 )
     vB = uB / u0
     v2 = u2 / u02
     UL(mx:mz) = u0 * ( VL(i,j,ro) + 4.d0 * VL(i,j,pr) ) * VL(i,j,ux:uz) + (B2*VL(i,j,ux:uz) - uB*VL(i,j,bx:bz))/u0
     UL(en) = u02*(VL(i,j,ro)+4.d0*VL(i,j,pr))-VL(i,j,pr) + B2 - 0.5d0*(B2+uB**2)/(u02)
     UL(de) = u0 * VL(i,j,ro)
     UL(bx:bz) = VL(i,j,bx:bz)
     UL(ps) = VL(i,j,ps)

!    VL --> FL
     W = ( VL(i,j,ro) + 4.d0 * VL(i,j,pr) )
     FL(de) = VL(i,j,ro) * VL(i,j,un)
     FL(mn)  = (W+B2/u02) * VL(i,j,un)*VL(i,j,un) + VL(i,j,pr) &
          - 2*vB*( VL(i,j,bn)*VL(i,j,un)/u0 ) &
          + 0.5d0*( B2*(1.d0-v2) + vB**2 )
     FL(mt1) = (W+B2/u02) * VL(i,j,un)*VL(i,j,ut1) - (1-v2)*VL(i,j,bn)*VL(i,j,bt1) &
          - vB*( VL(i,j,bn)*VL(i,j,ut1)+VL(i,j,bt1)*VL(i,j,un) )/u0
     FL(mt2) = (W+B2/u02) * VL(i,j,un)*VL(i,j,ut2) - (1-v2)*VL(i,j,bn)*VL(i,j,bt2) &
          - vB*( VL(i,j,bn)*VL(i,j,ut2)+VL(i,j,bt2)*VL(i,j,un) )/u0
     FL(en)  = u0*(W+B2/u02) * VL(i,j,un) - vB * VL(i,j,bn)
!     FL(bn) = 0.d0
     FL(bt1) = ( VL(i,j,un)*VL(i,j,bt1)-VL(i,j,bn)*VL(i,j,ut1) )/u0
     FL(bt2) = ( VL(i,j,un)*VL(i,j,bt2)-VL(i,j,bn)*VL(i,j,ut2) )/u0

!    VL -> coefficients
!    f1: cs^2 = (Gamma p / w) = (4p)/(3ro+12p)
!    f3: Q/w
     f1L = 4*VL(i,j,pr) / ( 3*VL(i,j,ro)+12*VL(i,j,pr) )
     f3L = ( B2-v2*B2+vB**2-f1L*(vB**2) )/( VL(i,j,ro)+4*VL(i,j,pr) )
     a2L = f1L + u02*(1-f1L) + f3L
     a1L = -2*u0*VL(i,j,un)*(1-f1L)
!     a0L = -f1L + u02*(VLv(i,1)**2)*(1-f1L) - f3L
     a0L = -f1L + (VL(i,j,un)**2)*(1-f1L) - f3L


!    VR --> UR
     B2 = dot_product( VR(i,j,bx:bz), VR(i,j,bx:bz) )
     uB = dot_product( VR(i,j,bx:bz), VR(i,j,ux:uz) )
     u2 = dot_product( VR(i,j,ux:uz), VR(i,j,ux:uz) )
     u02= 1.d0 + u2
     u0 = sqrt( u02 )
     vB = uB / u0
     v2 = u2 / u02
     UR(mx:mz) = u0 * ( VR(i,j,ro) + 4.d0 * VR(i,j,pr) ) * VR(i,j,ux:uz) + (B2*VR(i,j,ux:uz) - uB*VR(i,j,bx:bz))/u0
     UR(en) = u02*(VR(i,j,ro)+4.d0*VR(i,j,pr))-VR(i,j,pr) + B2 - 0.5d0*(B2+uB**2)/(u02)
     UR(de) = u0 * VR(i,j,ro)
     UR(bx:bz) = VL(i,j,bx:bz)
     UR(ps) = VL(i,j,ps)

!    VR --> FR
     W = ( VR(i,j,ro) + 4.d0 * VR(i,j,pr) )
     FR(de) = VR(i,j,ro) * VR(i,j,un)
     FR(mn)  = (W+B2/u02) * VR(i,j,un)*VR(i,j,un) + VR(i,j,pr) &
          - 2*vB*( VR(i,j,bn)*VR(i,j,un)/u0 ) &
          + 0.5d0*( B2*(1.d0-v2) + vB**2 )
     FR(mt1) = (W+B2/u02) * VR(i,j,un)*VR(i,j,ut1) - (1-v2)*VR(i,j,bn)*VR(i,j,bt1) &
          - vB*( VR(i,j,bn)*VR(i,j,ut1)+VR(i,j,bt1)*VR(i,j,un) )/u0
     FR(mt2) = (W+B2/u02) * VR(i,j,un)*VR(i,j,ut2) - (1-v2)*VR(i,j,bn)*VR(i,j,bt2) &
          - vB*( VR(i,j,bn)*VR(i,j,ut2)+VR(i,j,bt2)*VR(i,j,un) )/u0
     FR(en)  = u0*(W+B2/u02) * VR(i,j,un) - vB * VR(i,j,bn)
!     FR(bn) = 0.d0
     FR(bt1) = ( VR(i,j,un)*VR(i,j,bt1)-VR(i,j,bn)*VR(i,j,ut1) )/u0
     FR(bt2) = ( VR(i,j,un)*VR(i,j,bt2)-VR(i,j,bn)*VR(i,j,ut2) )/u0

!    VR -> coefficients
     f1R = 4*VR(i,j,pr) / ( 3*VR(i,j,ro)+12*VR(i,j,pr) )
     f3R = ( B2-v2*B2+vB**2-f1R*(vB**2) )/( VR(i,j,ro)+4*VR(i,j,pr) )
     a2R = f1R + u02*(1-f1R) + f3R
     a1R = -2*u0*VR(i,j,un)*(1-f1R)
!     a0R = -f1R + u02*(VRv(i,1)**2)*(1-f1R) - f3R
     a0R = -f1R + (VR(i,j,un)**2)*(1-f1R) - f3R
!! ========== This requires B=0 ========== !!

!    Riemann fan speed
     aL = min( &
          (-a1L-sqrt(a1L**2-4*a2L*a0L))/(2*a2L), &
          (-a1R-sqrt(a1R**2-4*a2R*a0R))/(2*a2R) )
     aR = max( &
          (-a1L+sqrt(a1L**2-4*a2L*a0L))/(2*a2L), &
          (-a1R+sqrt(a1R**2-4*a2R*a0R))/(2*a2R) )
     aLR = aL*aR

! ========== standard HLL ==========
!    F = F(L)
     if ( aL .gt. 0 ) then
        Fd(i)   = FL(de)
        Fm(i,1) = FL(mn)
        Fm(i,2) = FL(mt1)
        Fm(i,3) = FL(mt2)
        Fe(i)   = FL(en)
!        Fb(i,1) = FL(bn)
        Fb(i,2) = FL(bt1)
        Fb(i,3) = FL(bt2)
!    F = F(R)
     elseif ( aR .lt. 0 ) then
        Fd(i)   = FR(de)
        Fm(i,1) = FR(mn)
        Fm(i,2) = FR(mt1)
        Fm(i,3) = FR(mt2)
        Fe(i)   = FR(en)
!        Fb(i,1) = FR(bn)
        Fb(i,2) = FR(bt1)
        Fb(i,3) = FR(bt2)
!    F = F(HLL)
     else
        a = 1.d0 / ( aR - aL )
        Fd(i)  =a*(aR*FL(de) -aL*FR(de) +aLR*(UR(de) -UL(de)  ))
        Fm(i,1)=a*(aR*FL(mn) -aL*FR(mn) +aLR*(UR(mn) -UL(mn)))
        Fm(i,2)=a*(aR*FL(mt1)-aL*FR(mt1)+aLR*(UR(mt1)-UL(mt1)))
        Fm(i,3)=a*(aR*FL(mt2)-aL*FR(mt2)+aLR*(UR(mt2)-UL(mt2)))
        Fe(i)  =a*(aR*FL(en) -aL*FR(en) +aLR*(UR(en) -UL(en)  ))
!        Fb(i,1)=a*(aR*FL(bn) -aL*FR(bn) +aLR*(UR(bn) -UL(bn)))
        Fb(i,2)=a*(aR*FL(bt1)-aL*FR(bt1)+aLR*(UR(bt1)-UL(bt1)))
        Fb(i,3)=a*(aR*FL(bt2)-aL*FR(bt2)+aLR*(UR(bt2)-UL(bt2)))

! ========== HLLD module ==========
        if( hlld ) then
           ehll  = a*(aR*UR(en)-aL*UL(en)+FL(en)-FR(en))
           mxhll = a*(aR*UR(mn)-aL*UL(mn)+FL(mn)-FR(mn))
           
!       calculate covariant pressure
!       pt = p + 0.5 |b|^2
           a1 = ehll - Fm(i,1)
           a0 = mxhll * Fe(i) - Fm(i,1) * ehll
           pt = (-a1+sqrt(a1**2-4*a0)) / 2.d0

!       aM : speed of contact/tangential discon.
           aM = ( pt + aR*UR(mn) - FR(mn) ) &
                / ( aR*( UR(en)+pt ) - FR(en) )
           if ( ( aM .lt. aL).or.( aM .gt. aR )) then
              write(6,*) 'something is wrong in HLLD: ', aL, aM, aR
              stop
           endif
!          use HLLD only when |aR-aM| or |aM-aL| > threshold
           if( ( (aR-aM).gt.threshold ).and.( (aM-aL).gt.threshold ) ) then
!              write(6,*) 'HLLD'
!             HLLF flux : F_aL = FL + aL*(U_aL-UL)
!             Anyway, we are dealing with discrete numerical values
!             f2 = lambda * p + Re = E + p
              if ( aM .ge. 0 ) then
                 f1 = 1.d0 / ( aL - aM )
                 Rd   = ( aL*UL(de)  - FL(de)  )
                 Rmn  = ( aL*UL(mn)  - FL(mn)  )
                 Rmt1 = ( aL*UL(mt1) - FL(mt1) )
                 Rmt2 = ( aL*UL(mt2) - FL(mt2) )
                 Re   = ( aL*UL(en)  - FL(en)  )
                 Rbt1 = ( aL*UL(bt1) - FL(bt1) )
                 Rbt2 = ( aL*UL(bt2) - FL(bt2) )
                 mybt1 = f1*Rbt1
                 mybt2 = f1*Rbt2
                 f2 = aL*pt+Re
                 at1 = ( (f2-mybt2*Rbt2)*Rmt1+mybt1*Rbt2*Rmt2 ) &
                      / ( f2*( f2-mybt1*Rbt1-mybt2*Rbt2 ))
                 at2 = ( mybt2*Rbt1*Rmt1+(f2-mybt1*Rbt1)*Rmt2 ) &
                      / ( f2*( f2-mybt1*Rbt1-mybt2*Rbt2 ))
                 vB = at1*mybt1+at2*mybt2
                 Fd(i)   = f1*Rd*aM
                 Fm(i,1) = aL*( f1*f2*aM ) - Rmn
                 Fm(i,2) = aL*( f1*f2*at1 - vB*mybt1 ) - Rmt1
                 Fm(i,3) = aL*( f1*f2*at2 - vB*mybt2 ) - Rmt2
                 Fe(i)   = f1*f2*aM
                 Fb(i,1) = 0.d0
                 Fb(i,2) = mybt1*aM
                 Fb(i,3) = mybt2*aM
!             HLLF flux : F_aR = FR + aR*(U_aR-UR)
!             notice that F_aR != F( U_aR )
              else
                 f1 = 1.d0 / ( aR - aM )
                 Rd   = ( aL*UR(de)  - FR(de)  )
                 Rmn  = ( aL*UR(mn)  - FR(mn)  )
                 Rmt1 = ( aL*UR(mt1) - FR(mt1) )
                 Rmt2 = ( aL*UR(mt2) - FR(mt2) )
                 Re   = ( aL*UR(en)  - FR(en)  )
                 Rbt1 = ( aL*UR(bt1) - FR(bt1) )
                 Rbt2 = ( aL*UR(bt2) - FR(bt2) )
                 mybt1 = f1*Rbt1
                 mybt2 = f1*Rbt2
                 f2 = aR*pt+Re
                 at1 = ( (f2-mybt2*Rbt2)*Rmt1+mybt1*Rbt2*Rmt2 ) &
                      / ( f2*( f2-mybt1*Rbt1-mybt2*Rbt2 ))
                 at2 = ( mybt2*Rbt1*Rmt1+(f2-mybt1*Rbt1)*Rmt2 ) &
                      / ( f2*( f2-mybt1*Rbt1-mybt2*Rbt2 ))
                 vB = at1*mybt1+at2*mybt2
                 Fd(i)   = f1*Rd*aM
                 Fm(i,1) = aR*( f1*f2*aM ) - Rmn
                 Fm(i,2) = aR*( f1*f2*at1 - vB*mybt1 ) - Rmt1
                 Fm(i,3) = aR*( f1*f2*at2 - vB*mybt2 ) - Rmt2
                 Fe(i)   = f1*f2*aM
!                 Fb(i,1) = 0.d0
                 Fb(i,2) = mybt1*aM
                 Fb(i,3) = mybt2*aM
              endif
           endif
        endif
! ========== HLLD module ==========

     endif
! ========== standard HLL ==========

     F(i,j,mx:mz) = Fm(i,1:3)
     F(i,j,en)    = Fe(i)
     F(i,j,de)    = Fd(i)
     F(i,j,bt1)   = Fb(i,2)
     F(i,j,bt2)   = Fb(i,3)

  enddo

  return
end subroutine flux_solver
