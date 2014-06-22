subroutine hlld_g(G,VL,VR,ix,jx)
!-----------------------------------------------------------------------
!     HLLD flux solver in the Y direction
!       Ref: T. Miyoshi, K. Kusano, J. Comput. Phys., 208, 315 (2005)
!-----------------------------------------------------------------------
!     2010/05/11  S. Zenitani  HLLC-G solver
!     2010/05/13  S. Zenitani  HLLD solver
!     2010/09/18  S. Zenitani  fixed some bugs
!     2014/04/11  S. Zenitani  fixed a NaN problem in vfL and vfR
!-----------------------------------------------------------------------
  implicit none
  include 'param.h'
  integer, intent(in) :: ix, jx
! numerical flux (F) [output]
  real(8) :: G(ix,jx,var1)
! left/right states (VL,VR) [input]
  real(8) :: VL(ix,jx,var1), VR(ix,jx,var1)
! left/right conserved variables (UL & UR; local)
  real(8) :: UL(ix,jx,var2), UR(ix,jx,var2)
! numerical flux (GL & GR local)
  real(8) :: GL(ix,jx,var1), GR(ix,jx,var1)
  integer :: i, j

  real(8), parameter :: hllg_factor = 1.001   ! threshold to switch to HLLC-G

  real(8) :: B2, f1, f2
  real(8) :: aL, aL1, aM, aR1, aR, aVL, aVR, vBL, vBR
  real(8) :: vfL, vfR
  real(8) :: ptL, ptR
  real(8) :: UL1(var1), UR1(var1), U2(var1), U_hll(var1), pt
  real(8) :: ro_L1, ro_R1, vx_L1, vx_R1, vz_L1, vz_R1, vx_2, vz_2
  real(8) :: ro_Ls, ro_Rs
!-----------------------------------------------------------------------

  G(:,:,:) = 0.d0

  call v2u(VL,UL,ix,1,ix,jx,1,jx-1)
  call v2g(VL,GL,ix,1,ix,jx,1,jx-1)
  call v2u(VR,UR,ix,1,ix,jx,1,jx-1)
  call v2g(VR,GR,ix,1,ix,jx,1,jx-1)

  do j=1,jx-1
  do i=1,ix

!    VL -> coefficients
     B2 = dot_product( VL(i,j,bx:bz), VL(i,j,bx:bz) )
!    f1: Gamma p
!    f2: 4 gamma p B_n^2
     f1 = gamma * VL(i,j,pr)
     f2 = 4 * f1 * VL(i,j,by)**2
!    fast mode, total pressure, v dot B
!     vfL = sqrt( ( (f1+B2) + sqrt( (f1+B2)**2 - f2 )) / ( 2*VL(i,j,ro) ))
     vfL = sqrt( ( (f1+B2) + sqrt(max( (f1+B2)**2-f2, 0.d0 ))) / ( 2*VL(i,j,ro) ))
     ptL = VL(i,j,pr) + 0.5d0*B2
     vBL = dot_product( VL(i,j,vx:vz), VL(i,j,bx:bz) )

!    VR -> coefficients
     B2 = dot_product( VR(i,j,bx:bz), VR(i,j,bx:bz) )
!    f1: Gamma p
!    f2: 4 gamma p B_n^2
     f1 = gamma * VR(i,j,pr)
     f2 = 4 * f1 * VR(i,j,by)**2
!    fast mode, total pressure, v dot B
!     vfR = sqrt( ( (f1+B2) + sqrt( (f1+B2)**2 - f2 )) / ( 2*VR(i,j,ro) ))
     vfR = sqrt( ( (f1+B2) + sqrt(max( (f1+B2)**2-f2, 0.d0 ))) / ( 2*VR(i,j,ro) ))
     ptR = VR(i,j,pr) + 0.5d0*B2
     vBR = dot_product( VR(i,j,vx:vz), VR(i,j,bx:bz) )

!    Riemann fan speed (MK05 eq. 67)
!     aL = min( VL(i,j,vy) - vfL, VR(i,j,vy) - vfR )
!     aR = max( VL(i,j,vy) + vfL, VR(i,j,vy) + vfR )
     aL = min( VL(i,j,vy), VR(i,j,vy) ) - max( vfL, vfR )
     aR = max( VL(i,j,vy), VR(i,j,vy) ) + max( vfL, vfR )

!    G = G(L)
     if ( aL .ge. 0 ) then
        G(i,j,:) = GL(i,j,:)
!    G = G(R)
     elseif ( aR .le. 0 ) then
        G(i,j,:) = GR(i,j,:)

!    HLLC/HLLD flux
     else

!       HLL state
        f1 = 1.d0 / ( aR - aL )
        U_hll(mx:en) = f1*( aR*UR(i,j,mx:en) - aL*UL(i,j,mx:en) - GR(i,j,mx:en) + GL(i,j,mx:en) )
        U_hll(ro:bz) = f1*( aR*VR(i,j,ro:bz) - aL*VL(i,j,ro:bz) - GR(i,j,ro:bz) + GL(i,j,ro:bz) )

!       entropy wave
        aM = U_hll(my) / U_hll(ro)

!       Can we mathematically prove this?
        if ( aL .gt. aM .or. aR .lt. aM ) then
           write(6,*) 'error', aL, aM, aR
           write(6,*) ' fast mode: ', vfL, vfR
           stop
        endif

!       Total pressure
        pt    = ptL + VL(i,j,ro) * ( aL - VL(i,j,vy) ) * ( aM - VL(i,j,vy) )
        ro_L1 = VL(i,j,ro) * ( aL - VL(i,j,vy) ) / ( aL - aM ) ! ro(L*)
        ro_R1 = VR(i,j,ro) * ( aR - VR(i,j,vy) ) / ( aR - aM ) ! ro(R*)

!       For logical consistency, we employ By_hll as B_y in the intermediate states
        aVL = abs( U_hll(by) )/sqrt( ro_L1 )  ! Alfven wave (L)
        aVR = abs( U_hll(by) )/sqrt( ro_R1 )  ! Alfven wave (R)

! ========== revert to HLLC-G ==========
        if ( ( aL .ge. ( aM - hllg_factor*aVL ) ) .or. &
             ( ( aM + hllg_factor*aVR ) .ge. aR ) ) then

           vx_2 = U_hll(mx) / U_hll(ro) ! vx_HLL
           vz_2 = U_hll(mz) / U_hll(ro) ! vz_HLL

!          G = G(L*)
           if ( aM .ge. 0 ) then
              
              U2(en) = ( ( aL-VL(i,j,vy) )*UL(i,j,en) - ptL*VL(i,j,vy) + pt*aM + &
                   U_hll(by)*( vBL - vx_2*U_hll(bx) - aM*U_hll(by) - vz_2*U_hll(bz) ) ) / ( aL - aM )
              U2(mx) = ro_L1 * vx_2
              U2(my) = ro_L1 * aM
              U2(mz) = ro_L1 * vz_2

              G(i,j,mx:en) = aL *( U2(mx:en) - UL(i,j,mx:en) ) + GL(i,j,mx:en)
              G(i,j,ro)    = aL *( ro_L1     - VL(i,j,ro)    ) + GL(i,j,ro)
              G(i,j,bx:bz) = aL *( U_hll(bx:bz) - VL(i,j,bx:bz) ) + GL(i,j,bx:bz)

!          G = G(R*)
           else
              
              U2(en) = ( ( aR-VR(i,j,vy) )*UR(i,j,en) - ptR*VR(i,j,vy) + pt*aM + &
                   U_hll(by)*( vBR - vx_2*U_hll(bx) - aM*U_hll(by) - vz_2*U_hll(bz) ) ) / ( aR - aM )
              U2(mx) = ro_R1 * vx_2
              U2(my) = ro_R1 * aM
              U2(mz) = ro_R1 * vz_2

              G(i,j,mx:en) = aR *( U2(mx:en) - UR(i,j,mx:en) ) + GR(i,j,mx:en)
              G(i,j,ro)    = aR *( ro_R1     - VR(i,j,ro)    ) + GR(i,j,ro)
              G(i,j,bx:bz) = aR *( U_hll(bx:bz) - VR(i,j,bx:bz) ) + GR(i,j,bx:bz)

           endif

! ========== HLLD flux ==========
        else

!          Intermediate state (L*)
           f1 = 1.d0 / ( VL(i,j,ro)*(aL-VL(i,j,vy))*(aL-aM) - U_hll(by)**2 )
           UL1(bx) = VL(i,j,bx) * f1 * ( VL(i,j,ro)*(aL-VL(i,j,vy))**2 - U_hll(by)**2 )
           UL1(by) = U_hll(by)
           UL1(bz) = VL(i,j,bz) * f1 * ( VL(i,j,ro)*(aL-VL(i,j,vy))**2 - U_hll(by)**2 )
           vx_L1   = VL(i,j,vx) - f1 * U_hll(by)*VL(i,j,bx)*(aM-VL(i,j,vy))
           vz_L1   = VL(i,j,vz) - f1 * U_hll(by)*VL(i,j,bz)*(aM-VL(i,j,vy))
           
           UL1(ro) = ro_L1
           UL1(en) = ( ( aL-VL(i,j,vy) )*UL(i,j,en) - ptL*VL(i,j,vy) + pt*aM + &
                U_hll(by)*( vBL - vx_L1*UL1(bx) - aM*U_hll(by) - vz_L1*UL1(bz)) ) / ( aL - aM )
           UL1(mx) = ro_L1 * vx_L1
           UL1(my) = ro_L1 * aM
           UL1(mz) = ro_L1 * vz_L1

!          Intermediate state (R*)
           f1 = 1.d0 / ( VR(i,j,ro)*(aR-VR(i,j,vy))*(aR-aM) - U_hll(by)**2 )
           UR1(bx) = VR(i,j,bx) * f1 * ( VR(i,j,ro)*(aR-VR(i,j,vy))**2 - U_hll(by)**2 )
           UR1(by) = U_hll(by)
           UR1(bz) = VR(i,j,bz) * f1 * ( VR(i,j,ro)*(aR-VR(i,j,vy))**2 - U_hll(by)**2 )
           vx_R1   = VR(i,j,vx) - f1 * U_hll(by)*VR(i,j,bx)*(aM-VR(i,j,vy))
           vz_R1   = VR(i,j,vz) - f1 * U_hll(by)*VR(i,j,bz)*(aM-VR(i,j,vy))
           
           UR1(ro) = ro_R1
           UR1(en) = ( ( aR-VR(i,j,vy) )*UR(i,j,en) - ptR*VR(i,j,vy) + pt*aM + &
                U_hll(by)*( vBR - vx_R1*UR1(bx) - aM*U_hll(by) - vz_R1*UR1(bz)) ) / ( aR - aM )
           UR1(mx) = ro_R1 * vx_R1
           UR1(my) = ro_R1 * aM
           UR1(mz) = ro_R1 * vz_R1
           
!          rotational waves
           aL1 = aM - aVL
           aR1 = aM + aVR

!          G = G(L*)
           if( aL1 .ge. 0 ) then

              G(i,j,mx:en) = aL *( UL1(mx:en) - UL(i,j,mx:en) ) + GL(i,j,mx:en)
              G(i,j,ro)    = aL *( UL1(ro)    - VL(i,j,ro)    ) + GL(i,j,ro)
              G(i,j,bx:bz) = aL *( UL1(bx:bz) - VL(i,j,bx:bz) ) + GL(i,j,bx:bz)

!          G = G(R*)
           elseif( aR1 .le. 0 ) then

              G(i,j,mx:en) = aR *( UR1(mx:en) - UR(i,j,mx:en) ) + GR(i,j,mx:en)
              G(i,j,ro)    = aR *( UR1(ro)    - VR(i,j,ro) )    + GR(i,j,ro)
              G(i,j,bx:bz) = aR *( UR1(bx:bz) - VR(i,j,bx:bz) ) + GR(i,j,bx:bz)

!          Centeral states are tricky
!          Question: can we really rule out U_hll(by) == 0 ?
           else

              ro_Ls = sqrt( ro_L1 )  ! sqrt(ro(L*))
              ro_Rs = sqrt( ro_R1 )  ! sqrt(ro(R*))
              f1    = 1.d0 / ( ro_Ls + ro_Rs )
              f2    = dsign( 1.d0, U_hll(by) )
              vx_2   = f1 * ( ro_Ls*vx_L1 + ro_Rs*vx_R1 + ( UR1(bx)-UL1(bx) )*f2 )
              vz_2   = f1 * ( ro_Ls*vz_L1 + ro_Rs*vz_R1 + ( UR1(bz)-UL1(bz) )*f2 )
              U2(bx) = f1 * ( ro_Ls*UR1(bx) + ro_Rs*UL1(bx) + ro_Ls*ro_Rs*( vx_R1-vx_L1 )*f2 )
              U2(by) = U_hll(by)
              U2(bz) = f1 * ( ro_Ls*UR1(bz) + ro_Rs*UL1(bz) + ro_Ls*ro_Rs*( vz_R1-vz_L1 )*f2 )

!             G = G(L**)
              if( aM .ge. 0 ) then

                 U2(ro) = ro_L1
                 U2(en) = UL1(en) - ro_Ls * ( vx_L1*UL1(bx) + vz_L1*UL1(bz) &
                      - vx_2*U2(bx) - vz_2*U2(bz) ) * f2
                 U2(mx) = ro_L1 * vx_2
                 U2(my) = ro_L1 * aM
                 U2(mz) = ro_L1 * vz_2

                 G(i,j,mx:en) = aL1*U2(mx:en) - (aL1-aL)*UL1(mx:en) - aL*UL(i,j,mx:en) + GL(i,j,mx:en)
                 G(i,j,ro)    = aL1*U2(ro)    - (aL1-aL)*UL1(ro)    - aL*VL(i,j,ro)    + GL(i,j,ro)
                 G(i,j,bx:bz) = aL1*U2(bx:bz) - (aL1-aL)*UL1(bx:bz) - aL*VL(i,j,bx:bz) + GL(i,j,bx:bz)

!             G = G(R**)
              else

                 U2(ro) = ro_R1
                 U2(en) = UR1(en) + ro_Rs * ( vx_R1*UR1(bx) + vz_R1*UR1(bz) &
                      - vx_2*U2(bx) - vz_2*U2(bz) ) * f2
                 U2(mx) = ro_R1 * vx_2
                 U2(my) = ro_R1 * aM
                 U2(mz) = ro_R1 * vz_2

                 G(i,j,mx:en) = aR1*U2(mx:en) - (aR1-aR)*UR1(mx:en) - aR*UR(i,j,mx:en) + GR(i,j,mx:en)
                 G(i,j,ro)    = aR1*U2(ro)    - (aR1-aR)*UR1(ro)    - aR*VR(i,j,ro)    + GR(i,j,ro)
                 G(i,j,bx:bz) = aR1*U2(bx:bz) - (aR1-aR)*UR1(bx:bz) - aR*VR(i,j,bx:bz) + GR(i,j,bx:bz)

              endif

           endif

        endif

     endif

  enddo
  enddo

  return
end subroutine hlld_g
