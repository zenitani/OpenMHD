subroutine u2v(U,V,ix,jx)
!-----------------------------------------------------------------------
!     Gamma-law version of Mignone & McKinney (2007)'s inversion scheme
!     See also PLUTO (Mignone et al. 2007 ApJS)
!-----------------------------------------------------------------------
!     2009/07/27  S. Zenitani
!-----------------------------------------------------------------------
  implicit none
  include 'param_rela.h'
  integer, parameter :: loop_max = 50        ! Usually 10-20 is OK.
  real(8), parameter :: threshold = 1.0d-13  ! May be computer-dependent
!-----------------------------------------------------------------------
  integer, intent(in)  :: ix, jx
  real(8), intent(in)  :: U(ix,jx,var1)  ! conserved variables  [input]
  real(8), intent(out) :: V(ix,jx,var2)  ! primitive variables  [output]
!-----------------------------------------------------------------------
  real(8) :: Ud     ! density
  real(8) :: Um(3)  ! momentum density
  real(8) :: Ue     ! energy density
  real(8) :: Ub(3)  ! magnetic field
  real(8) :: Vv(3)  ! 3-velocity (tmp)
!-----------------------------------------------------------------------
  integer :: i, j, k, izero
  real(8) :: B2, M2, S, S2, W, W2, W3, S2_W2
  real(8) :: f1, f2, v2, u0, u02
  real(8) :: dv2_dW, chi, dchi_dW, dfW, dW, fW
!  real(8) :: rro
  real(8) :: pre

  do j=1,jx
  do i=1,ix
     Ue = U(i,j,en)
     Ud = U(i,j,de)
     Um(1:3) = U(i,j,mx:mz)
     Ub(1:3) = U(i,j,bx:bz)
     B2 = dot_product( U(i,j,bx:bz), U(i,j,bx:bz) )
     M2 = dot_product( U(i,j,mx:mz), U(i,j,mx:mz) )
     S  = dot_product( U(i,j,bx:bz), U(i,j,mx:mz) )
     S2 = S*S
     
     f1 = - 4 * ( Ue - B2 )
     f2 = M2 + ( B2 - 2*Ue )*B2
     W = ( - f1 + sqrt(f1**2 - 12*f2))/6.d0
     if( W.lt.Ud ) then
        W = Ud
     elseif( W.gt.(4.d0/3.d0)*(Ue-B2/2.d0) ) then
        W = (4.d0/3.d0)*(Ue-B2/2.d0)
     endif
     
! ---- iteration loop ----
     do k=1,loop_max

        W2 = W*W
        W3 = W2*W
        S2_W2 = S2/W2
        f1 = W + B2

        v2 = (M2 + S2_W2*(f1 + W))/(f1*f1)
        if( v2.gt.1 )then
           write(6,*) 'stop: v2 > 1', k, 'th iteration'
           stop
        endif

        u02 = 1.d0 / (1.d0 - v2)
        u0  = sqrt(u02)
        
        dv2_dW  = ( S2*(3*W*f1 + B2*B2) + M2*W3 ) * (-2) / W3 / (f1**3)
        chi = (W - Ud*u0) * (1.d0 - v2)
!        rro = Ud/u0 ! unused
        dchi_dW =  1.d0 - v2 - 0.5d0*u0*(Ud + 2*chi*u0)*dv2_dW
        pre = 0.25d0*chi

        fW  = W - pre + 0.5d0*(1.d0 + v2)*B2 - 0.5d0*S2_W2 - Ue
        dfW = 1.d0 - 0.25d0*dchi_dW + 0.5d0*B2*dv2_dW + S2/W3
        dW = fW/dfW

        W = W - dW
!        write(6,*) k, ' W, dW, dW/W', W, dW, dW/W
        if (abs(dW).lt.(threshold*W)) then
           if( pre.gt.0.d0 ) then
              goto 1000
           endif
        endif
        if( k.eq.loop_max ) then
           write(6,*) 'No convergence'
           write(6,*) 'U: ', i, Ud, Ue, Um(1:3)
           stop
        endif

     enddo
1000 continue
! ---- iteration loop ----

!     write(6,*) 'error: ', abs(dW)/W

     Vv(1:3) = (1.d0/(W+B2))*( Um(1:3) + (S/W)*Ub(1:3) )
     v2 = dot_product( Vv(1:3), Vv(1:3) )
     if( v2 .ge. 1 ) then
        write(6,*) 'error: superluminous velocity', i, v2
        stop
     endif
     V(i,j,ux:uz) = Vv(1:3) / sqrt( 1.d0 - v2 )
     V(i,j,pr) = pre
     V(i,j,ro) = Ud*sqrt( 1.d0 - v2 )
     
  enddo
  enddo

  return
end subroutine u2v
