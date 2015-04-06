subroutine llf_g(G,VL,VR,ix,jx)
!-----------------------------------------------------------------------
!     Local Lax-Friedrich (Rusanov) solver in the Y direction
!-----------------------------------------------------------------------
!     2010/05/13  S. Zenitani  LLF solver
!     2010/09/18  S. Zenitani  fixed some bugs
!-----------------------------------------------------------------------
  implicit none
  include 'param.h'
  integer, intent(in) :: ix, jx
! numerical flux (G) [output]
  real(8), intent(out) :: G(ix,jx,var1)
! left/right states (VL,VR) [input]
  real(8) :: VL(ix,jx,var1), VR(ix,jx,var1)
! left/right conserved variables (UL & UR; local)
  real(8) :: UL(ix,jx,var2), UR(ix,jx,var2)
! numerical flux (GL & GR local)
  real(8) :: GL(ix,jx,var1), GR(ix,jx,var1)
  integer :: i, j

  real(8) :: B2, f1, f2
  real(8) :: aL, aR
  real(8) :: vfL, vfR

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
!    fast mode
!     vfL = sqrt( ( (f1+B2) + sqrt( (f1+B2)**2 - f2 )) / ( 2*VL(i,j,ro) ))
     vfL = sqrt( ( (f1+B2) + sqrt(max( (f1+B2)**2-f2, 0.d0 ))) / ( 2*VL(i,j,ro) ))

!    VR -> coefficients
     B2 = dot_product( VR(i,j,bx:bz), VR(i,j,bx:bz) )
!    f1: Gamma p
!    f2: 4 gamma p B_n^2
     f1 = gamma * VR(i,j,pr)
     f2 = 4 * f1 * VR(i,j,by)**2
!    fast mode
!     vfR = sqrt( ( (f1+B2) + sqrt( (f1+B2)**2 - f2 )) / ( 2*VR(i,j,ro) ))
     vfR = sqrt( ( (f1+B2) + sqrt(max( (f1+B2)**2-f2, 0.d0 ))) / ( 2*VR(i,j,ro) ))

!    Local Lax Friedrich
     aR = max( abs(VL(i,j,vy))+vfL, abs(VR(i,j,vy))+vfR )
     aL = -aR

!    G = G(L)
     if ( aL .ge. 0 ) then
        G(i,j,:) = GL(i,j,:)
!    G = G(R)
     elseif ( aR .le. 0 ) then
        G(i,j,:) = GR(i,j,:)
!    G = G(LLF)
     else

        f1 = 1.d0 / ( aR - aL )
        f2 = aL * aR

        G(i,j,mx:en) = f1*( aR*GL(i,j,mx:en) - aL*GR(i,j,mx:en) &
             + f2 *(UR(i,j,mx:en)-UL(i,j,mx:en)) )
        G(i,j,ro:bz) = f1*( aR*GL(i,j,ro:bz) - aL*GR(i,j,ro:bz) &
             + f2 *(VR(i,j,ro:bz)-VL(i,j,ro:bz)) )

     endif

  enddo
  enddo

  return
end subroutine llf_g
