subroutine set_dt(U,V,vmax,dt,dx,cfl,ix,jx,kx)
!-----------------------------------------------------------------------
!     CFL condition
!-----------------------------------------------------------------------
  implicit none
  include 'param.h'
  real(8), parameter :: dtmin = 1.0d-7
!-----------------------------------------------------------------------
  integer, intent(in) :: ix, jx, kx
  real(8), intent(in) :: U(ix,jx,kx,var1)
  real(8), intent(in) :: V(ix,jx,kx,var2)
  real(8), intent(in) :: dx, cfl
  real(8), intent(out) :: vmax, dt
!-----------------------------------------------------------------------
  integer :: i, j, k, is, ie, js, je, ks, ke
  integer :: imax, jmax, kmax, mymax(3)
  real(8) :: vtmp(ix,jx,kx)
  real(8) :: B2, f1, f2, vfx, vfy, vfz


  is = min(2,ix); ie = max(1,ix-1)
  js = min(2,jx); je = max(1,jx-1)
  ks = min(2,kx); ke = max(1,kx-1)
  vmax = -1.d0
! vtmp = 0.d0

!$omp parallel do private(i,j,k,B2,f1,f2,vfx,vfy,vfz) reduction(max:vmax)
  do k=ks,ke
  do j=js,je
  do i=is,ie

     B2 = dot_product( U(i,j,k,bx:bz), U(i,j,k,bx:bz) )
     f1 = gamma * V(i,j,k,pr)

!    fast mode in the X direction
     f2 = 4 * f1 * U(i,j,k,bx)**2
!     vfx = sqrt( ( (f1+B2) + sqrt( (f1+B2)**2 - f2 )) / ( 2*U(i,j,k,ro) ))
     vfx = sqrt( ( (f1+B2) + sqrt(max( (f1+B2)**2-f2, 0.d0 ))) / ( 2*U(i,j,k,ro) ))

!    fast mode in the Y direction
     f2 = 4 * f1 * U(i,j,k,by)**2
     vfy = sqrt( ( (f1+B2) + sqrt(max( (f1+B2)**2-f2, 0.d0 ))) / ( 2*U(i,j,k,ro) ))

!    fast mode in the Z direction
     f2 = 4 * f1 * U(i,j,k,bz)**2
     vfz = sqrt( ( (f1+B2) + sqrt(max( (f1+B2)**2-f2, 0.d0 ))) / ( 2*U(i,j,k,ro) ))

!    max speed of MHD waves
     vtmp(i,j,k) = max( abs(V(i,j,k,vx))+vfx, abs(V(i,j,k,vy))+vfy, abs(V(i,j,k,vz))+vfz )
     vmax = max( vtmp(i,j,k), vmax )

  enddo
  enddo
  enddo
!$omp end parallel do

! dt
  dt = cfl * dx / vmax

! error check
  if( dt < dtmin ) then
     if( vmax < 0.d0 ) then
        write(6,*) 'Potential error in OpenMP/reduction ...'
     endif
     mymax = maxloc(vtmp)
     imax = mymax(1); jmax = mymax(2); kmax = mymax(3)
     write(6,*) ' dt is too small : ', dt, ' < ', dtmin
     write(6,*) '     velocity is : ', vmax
     write(6,*) imax, jmax, kmax, U(imax,jmax,kmax,ro), V(imax,jmax,kmax,pr)
     stop
  endif

  return
end subroutine set_dt
