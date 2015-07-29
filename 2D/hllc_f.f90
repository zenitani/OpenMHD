subroutine hllc_f(F,VL,VR,ix,jx)
!-----------------------------------------------------------------------
!     HLLC-G solver in the X direction
!       Ref: K. F. Gurski, SIAM J. Sci. Comput., 25, 2165 (2004)
!-----------------------------------------------------------------------
!     2010/05/11  S. Zenitani  HLLC-G solver
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

  real(8) :: vB, B2, f1, f2
  real(8) :: aL, aR, a, aM
  real(8) :: vf, vfL2, vfR2
  real(8) :: ptL, ptR
  real(8) :: U_tmp(var1), ro_hll, pt_tmp

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
     ptL = VL(i,j,pr) + 0.5d0*B2

!    VR -> coefficients
     B2 = dot_product( VR(i,j,bx:bz), VR(i,j,bx:bz) )
!    f1: Gamma p
!    f2: 4 gamma p B_n^2
     f1 = gamma * VR(i,j,pr)
     f2 = 4 * f1 * VR(i,j,bx)**2
!    fast mode^2
!     vfR = sqrt( ( (f1+B2) + sqrt( (f1+B2)**2 - f2 )) / ( 2*VR(i,j,ro) ))
     vfR2 = ( (f1+B2) + sqrt(max( (f1+B2)**2-f2, 0.d0 ))) / ( 2*VR(i,j,ro) )
     ptR = VR(i,j,pr) + 0.5d0*B2

!    Riemann fan speed (MK05 eq. 67)
!     aL = min( VL(i,j,vx) - vfL, VR(i,j,vx) - vfR )
!     aR = max( VL(i,j,vx) + vfL, VR(i,j,vx) + vfR )
!     aL = min( VL(i,j,vx), VR(i,j,vx) ) - max( vfL, vfR )
!     aR = max( VL(i,j,vx), VR(i,j,vx) ) + max( vfL, vfR )
     vf = sqrt( max( vfL2, vfR2 ) )
     aL = min( VL(i,j,vx), VR(i,j,vx) ) - vf
     aR = max( VL(i,j,vx), VR(i,j,vx) ) + vf

!    entropy wave
     aM = ( &
          ( (aR-VR(i,j,vx))*VR(i,j,ro)*VR(i,j,vx) - ptR ) - &
          ( (aL-VL(i,j,vx))*VL(i,j,ro)*VL(i,j,vx) - ptL ) &
          ) / ( (aR-VR(i,j,vx))*VR(i,j,ro) - (aL-VL(i,j,vx))*VL(i,j,ro) )

!     if ( aL .gt. aM .or. aR .lt. aM ) then
!        write(6,*) 'error', aL, aM, aR
!        write(6,*) ' fast mode: ', vfL, vfR
!        stop
!     endif

!    F = F(L)
     if ( aL .ge. 0 ) then
        F(i,j,:) = FL(i,j,:)
!    F = F(R)
     elseif ( aR .le. 0 ) then
        F(i,j,:) = FR(i,j,:)
!    F = F(HLLC)
     else
        a = 1.d0 / ( aR - aL )

        U_tmp(my) = a*( aR*UR(i,j,my) - aL*UL(i,j,my) - FR(i,j,my) + FL(i,j,my) ) ! ( ro*vy )_HLL
        U_tmp(mz) = a*( aR*UR(i,j,mz) - aL*UL(i,j,mz) - FR(i,j,mz) + FL(i,j,mz) ) ! ( ro*vz )_HLL
        ro_hll    = a*( aR*VR(i,j,ro) - aL*VL(i,j,ro) - FR(i,j,ro) + FL(i,j,ro) )
        U_tmp(by) = a*( aR*VR(i,j,by) - aL*VL(i,j,by) - FR(i,j,by) + FL(i,j,by) ) ! By_hll
        U_tmp(bz) = a*( aR*VR(i,j,bz) - aL*VL(i,j,bz) - FR(i,j,bz) + FL(i,j,bz) ) ! Bz_hll

!       F = F(L*)
        if ( aM .ge. 0 ) then
           vB        = dot_product( VL(i,j,vx:vz), VL(i,j,bx:bz) )
           U_tmp(ro) = VL(i,j,ro) * ( aL - VL(i,j,vx) ) / ( aL - aM )
           pt_tmp    = ptL + VL(i,j,ro) * ( aL - VL(i,j,vx) ) * ( aM - VL(i,j,vx) )
           U_tmp(en) = ( ( aL - VL(i,j,vx) )*UL(i,j,en) - ptL*VL(i,j,vx) + pt_tmp * aM + &
                U_tmp(bx)*( vB - aM*U_tmp(bx) - (U_tmp(my)*U_tmp(by)+U_tmp(mz)*U_tmp(bz))/ro_hll ) ) /  ( aL - aM )
           U_tmp(mx) = U_tmp(ro) * aM
           U_tmp(my) = U_tmp(my) * U_tmp(ro) / ro_hll
           U_tmp(mz) = U_tmp(mz) * U_tmp(ro) / ro_hll
           F(i,j,mx:en) = aL *( U_tmp(mx:en) - UL(i,j,mx:en) ) + FL(i,j,mx:en)
           F(i,j,ro) = aL *( U_tmp(ro) - VL(i,j,ro) ) + FL(i,j,ro)
           F(i,j,by) = aL *( U_tmp(by) - VL(i,j,by) ) + FL(i,j,by)
           F(i,j,bz) = aL *( U_tmp(bz) - VL(i,j,bz) ) + FL(i,j,bz)
!       F = F(R*)
        else
           vB        = dot_product( VR(i,j,vx:vz), VR(i,j,bx:bz) )
           U_tmp(ro) = VR(i,j,ro) * ( aR - VR(i,j,vx) ) / ( aR - aM )
           pt_tmp    = ptR + VR(i,j,ro) * ( aR - VR(i,j,vx) ) * ( aM - VR(i,j,vx) )
           U_tmp(en) = ( ( aR - VR(i,j,vx) )*UR(i,j,en) - ptR*VR(i,j,vx) + pt_tmp * aM + &
                U_tmp(bx)*( vB - aM*U_tmp(bx) - (U_tmp(my)*U_tmp(by)+U_tmp(mz)*U_tmp(bz))/ro_hll ) ) /  ( aR - aM )
           U_tmp(mx) = U_tmp(ro) * aM
           U_tmp(my) = U_tmp(my) * U_tmp(ro) / ro_hll
           U_tmp(mz) = U_tmp(mz) * U_tmp(ro) / ro_hll
           F(i,j,mx:en) = aR *( U_tmp(mx:en) - UR(i,j,mx:en) ) + FR(i,j,mx:en)
           F(i,j,ro) = aR *( U_tmp(ro) - VR(i,j,ro) ) + FR(i,j,ro)
           F(i,j,by) = aR *( U_tmp(by) - VR(i,j,by) ) + FR(i,j,by)
           F(i,j,bz) = aR *( U_tmp(bz) - VR(i,j,bz) ) + FR(i,j,bz)
        endif

     endif

  enddo
  enddo

  return
end subroutine hllc_f
