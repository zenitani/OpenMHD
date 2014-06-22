subroutine input(filename,ix,jx,t,x,y,U,V)
  implicit none
  include 'param.h'
  character*256 :: filename      
  integer, intent(in) :: ix, jx
  real(8), intent(out) :: t
  real(8), intent(out) :: x(ix), y(jx)
  real(8), intent(out) :: U(ix,jx,var1)
  real(8), intent(out) :: V(ix,jx,var2)
  integer :: i, j, ix0, jx0
  real(8) :: f1, v2, B2

  U(:,:,:) = 0.d0
  V(:,:,:) = 0.d0

  open(15,file=filename,form='unformatted')
  read(15) t
  read(15) ix0
  read(15) jx0
  if(( ix0.ne.ix ) .or. ( jx0.ne.jx )) then
     write(6,*) 'parameter mismatch'
     write(6,*) ' ix= ',ix,' ix0= ',ix0
     write(6,*) ' jx= ',jx,' jx0= ',jx0
     stop
  endif
  read(15) (x(i),i=1,ix)
  read(15) (y(j),j=1,jx)
  read(15) ((U(i,j,mx),i=1,ix),j=1,jx)
  read(15) ((U(i,j,my),i=1,ix),j=1,jx)
  read(15) ((U(i,j,mz),i=1,ix),j=1,jx)
  read(15) ((U(i,j,en),i=1,ix),j=1,jx)
  read(15) ((U(i,j,ro),i=1,ix),j=1,jx)
  read(15) ((U(i,j,bx),i=1,ix),j=1,jx)
  read(15) ((U(i,j,by),i=1,ix),j=1,jx)
  read(15) ((U(i,j,bz),i=1,ix),j=1,jx)
  read(15) ((U(i,j,ps),i=1,ix),j=1,jx)
  read(15) ((V(i,j,vx),i=1,ix),j=1,jx)
  read(15) ((V(i,j,vy),i=1,ix),j=1,jx)
  read(15) ((V(i,j,vz),i=1,ix),j=1,jx)
  read(15) ((V(i,j,pr),i=1,ix),j=1,jx)
  close(15)

! ----------------------------------------------------------------------
! reconstructing V --> U
! I don't use the v2u routine due to the different martirx size
!  f1 = 1.d0 / ( gamma - 1 )
!  do j=1,jx
!  do i=1,ix
!     v2 = dot_product( V(i,j,vx:vz), V(i,j,vx:vz) )
!     B2 = dot_product( U(i,j,bx:bz), U(i,j,bx:bz) )
!     U(i,j,mx:mz) = U(i,j,ro) * V(i,j,vx:vz)
!     U(i,j,en)    = 0.5d0 * ( U(i,j,ro)*v2 + B2 ) + f1*V(i,j,pr)
!  enddo
!  enddo
! ----------------------------------------------------------------------

  return
end subroutine input
