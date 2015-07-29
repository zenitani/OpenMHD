subroutine hllc_g(G,VL,VR,ix,jx)
!-----------------------------------------------------------------------
!     HLLC-G solver in the Y direction
!       Ref: K. F. Gurski, SIAM J. Sci. Comput., 25, 2165 (2004)
!-----------------------------------------------------------------------
!     2010/05/11  S. Zenitani  HLLC-G solver
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

  real(8) :: vB, B2, f1, f2
  real(8) :: aL, aR, a, aM
  real(8) :: vf, vfL2, vfR2
  real(8) :: ptL, ptR
  real(8) :: U_tmp(var1), ro_hll, pt_tmp

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
     vf = sqrt( max( vfL2, vfR2 ) )
     aL = min( VL(i,j,vy), VR(i,j,vy) ) - vf
     aR = max( VL(i,j,vy), VR(i,j,vy) ) + vf

!    entropy wave
     aM = ( &
          ( (aR-VR(i,j,vy))*VR(i,j,ro)*VR(i,j,vy) - ptR ) - &
          ( (aL-VL(i,j,vy))*VL(i,j,ro)*VL(i,j,vy) - ptL ) &
          ) / ( (aR-VR(i,j,vy))*VR(i,j,ro) - (aL-VL(i,j,vy))*VL(i,j,ro) )

!     if ( aL .gt. aM .or. aR .lt. aM ) then
!        write(6,*) 'error', aL, aM, aR
!        write(6,*) ' fast mode: ', vfL, vfR
!        stop
!     endif


!    G = G(L)
     if ( aL .gt. 0 ) then
        G(i,j,:) = GL(i,j,:)
!    G = G(R)
     elseif ( aR .lt. 0 ) then
        G(i,j,:) = GR(i,j,:)
!    G = G(HLLC)
     else
        a = 1.d0 / ( aR - aL )

        U_tmp(mx) = a*( aR*UR(i,j,mx) - aL*UL(i,j,mx) - GR(i,j,mx) + GL(i,j,mx) ) ! ( ro*vx )_HLL
        U_tmp(mz) = a*( aR*UR(i,j,mz) - aL*UL(i,j,mz) - GR(i,j,mz) + GL(i,j,mz) ) ! ( ro*vz )_HLL
        ro_hll    = a*( aR*VR(i,j,ro) - aL*VL(i,j,ro) - GR(i,j,ro) + GL(i,j,ro) )
        U_tmp(bx) = a*( aR*VR(i,j,bx) - aL*VL(i,j,bx) - GR(i,j,bx) + GL(i,j,bx) ) ! Bx_hll
        U_tmp(bz) = a*( aR*VR(i,j,bz) - aL*VL(i,j,bz) - GR(i,j,bz) + GL(i,j,bz) ) ! Bz_hll

!       G = G(L*)
        if ( aM .ge. 0 ) then
           vB        = dot_product( VL(i,j,vx:vz), VL(i,j,bx:bz) )
           U_tmp(ro) = VL(i,j,ro) * ( aL - VL(i,j,vy) ) / ( aL - aM )
           pt_tmp    = ptL + VL(i,j,ro) * ( aL - VL(i,j,vy) ) * ( aM - VL(i,j,vy) )
           U_tmp(en) = ( ( aL - VL(i,j,vy) )*UL(i,j,en) - ptL*VL(i,j,vy) + pt_tmp * aM + &
                VL(i,j,by)*( vB - aM*VL(i,j,by) - (U_tmp(mx)*U_tmp(bx)+U_tmp(mz)*U_tmp(bz))/ro_hll ) ) /  ( aL - aM )
           U_tmp(mx) = U_tmp(mx) * U_tmp(ro) / ro_hll
           U_tmp(my) = U_tmp(ro) * aM
           U_tmp(mz) = U_tmp(mz) * U_tmp(ro) / ro_hll
           G(i,j,mx:en) = aL *( U_tmp(mx:en) - UL(i,j,mx:en) ) + GL(i,j,mx:en)
           G(i,j,ro) = aL *( U_tmp(ro) - VL(i,j,ro) ) + GL(i,j,ro)
           G(i,j,bx) = aL *( U_tmp(bx) - VL(i,j,bx) ) + GL(i,j,bx)
           G(i,j,bz) = aL *( U_tmp(bz) - VL(i,j,bz) ) + GL(i,j,bz)

!       G = G(R*)
        else
           vB        = dot_product( VR(i,j,vx:vz), VR(i,j,bx:bz) )
           U_tmp(ro) = VR(i,j,ro) * ( aR - VR(i,j,vy) ) / ( aR - aM )
           pt_tmp    = ptR + VR(i,j,ro) * ( aR - VR(i,j,vy) ) * ( aM - VR(i,j,vy) )
           U_tmp(en) = ( ( aR - VR(i,j,vy) )*UR(i,j,en) - ptR*VR(i,j,vy) + pt_tmp * aM + &
                VR(i,j,by)*( vB - aM*VR(i,j,by) - (U_tmp(mx)*U_tmp(bx)+U_tmp(mz)*U_tmp(bz))/ro_hll ) ) /  ( aR - aM )
           U_tmp(mx) = U_tmp(mx) * U_tmp(ro) / ro_hll
           U_tmp(my) = U_tmp(ro) * aM
           U_tmp(mz) = U_tmp(mz) * U_tmp(ro) / ro_hll
           G(i,j,mx:en) = aR *( U_tmp(mx:en) - UR(i,j,mx:en) ) + GR(i,j,mx:en)
           G(i,j,ro) = aR *( U_tmp(ro) - VR(i,j,ro) ) + GR(i,j,ro)
           G(i,j,bx) = aR *( U_tmp(bx) - VR(i,j,bx) ) + GR(i,j,bx)
           G(i,j,bz) = aR *( U_tmp(bz) - VR(i,j,bz) ) + GR(i,j,bz)

        endif

     endif

  enddo
  enddo

  return
end subroutine hllc_g
