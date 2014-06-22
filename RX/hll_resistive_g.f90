subroutine hll_resistive_g(G,U,VL,VR,EtS,dx,ix,jx)
!-----------------------------------------------------------------------
!     Resistive HLL solver in the Y direction
!-----------------------------------------------------------------------
!     2010/09/23  S. Zenitani  resistive HLL solver
!     2010/11/30  S. Zenitani  fixed bug in J_x
!     2014/04/11  S. Zenitani  fixed a NaN problem in vfL and vfR
!-----------------------------------------------------------------------
  implicit none
  include 'param.h'
!-----------------------------------------------------------------------
  integer, intent(in) :: ix, jx
  real(8), intent(in) :: dx
! numerical flux (G) [output]
  real(8), intent(out) :: G(ix,jx,var1)
! left/right states (VL,VR) [input]
  real(8), intent(in) :: VL(ix,jx,var1), VR(ix,jx,var1)
! conserved variables (U) [input]
  real(8), intent(in) :: U(ix,jx,var1)
!-----------------------------------------------------------------------
! left/right conserved variables (UL & UR; local)
  real(8) :: UL(ix,jx,var2), UR(ix,jx,var2)
! numerical flux (GL & GR local)
  real(8) :: GL(ix,jx,var1), GR(ix,jx,var1)
! electric current, magnetic fields, and resistivity at the cell surface
! J and B are local, while EtS comes from the input
  real(8) :: JxS(ix,jx), JzS(ix,jx), EtS(ix,jx)
  integer :: i, j

  real(8) :: B2, f1, f2
  real(8) :: aL, aR
  real(8) :: vfL, vfR

  G(:,:,:) = 0.d0
  JxS(:,:) = 0.d0
  JzS(:,:) = 0.d0

  call v2u(VL,UL,ix,1,ix,jx,1,jx-1)
  call v2g(VL,GL,ix,1,ix,jx,1,jx-1)
  call v2u(VR,UR,ix,1,ix,jx,1,jx-1)
  call v2g(VR,GR,ix,1,ix,jx,1,jx-1)

! surface current (Toth+ 2008, JCP)
  f1 = 1.d0 / dx
  do j=1,jx-1
     do i=2,ix-1
        JxS(i,j) = f1*( U(i,j+1,bz)-U(i,j,bz) )
!        JyS(i,j) = -0.25d0*f1*( U(i+1,j+1,bz)+U(i+1,j,bz)-U(i-1,j+1,bz)-U(i-1,j,bz) )
        JzS(i,j) =  0.25d0*f1*( U(i+1,j+1,by)+U(i+1,j,by)-U(i-1,j+1,by)-U(i-1,j,by) ) &
             - f1*( U(i,j+1,bx)-U(i,j,bx) )
     enddo
  enddo

  do j=1,jx-1
     GL(:,j,en) = GL(:,j,en) + EtS(:,j) * ( JzS(:,j)*VL(:,j,bx) - JxS(:,j)*VL(:,j,bz) )
     GL(:,j,bx) = GL(:,j,bx) + EtS(:,j) * JzS(:,j)
     GL(:,j,bz) = GL(:,j,bz) - EtS(:,j) * JxS(:,j)
     GR(:,j,en) = GR(:,j,en) + EtS(:,j) * ( JzS(:,j)*VR(:,j,bx) - JxS(:,j)*VR(:,j,bz) )
     GR(:,j,bx) = GR(:,j,bx) + EtS(:,j) * JzS(:,j)
     GR(:,j,bz) = GR(:,j,bz) - EtS(:,j) * JxS(:,j)
  enddo

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

!    Riemann fan speed (MK05 eq. 67)
!     aL = min( VL(i,j,vy) - vfL, VR(i,j,vy) - vfR )
!     aR = max( VL(i,j,vy) + vfL, VR(i,j,vy) + vfR )
     aL = min( VL(i,j,vy), VR(i,j,vy) ) - max( vfL, vfR )
     aR = max( VL(i,j,vy), VR(i,j,vy) ) + max( vfL, vfR )

!    F = F(L)
     if ( aL .ge. 0 ) then
        G(i,j,:) = GL(i,j,:)
!    F = F(R)
     elseif ( aR .le. 0 ) then
        G(i,j,:) = GR(i,j,:)
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
end subroutine hll_resistive_g
