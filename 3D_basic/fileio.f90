!-----------------------------------------------------------------------
!    File IO routines
!-----------------------------------------------------------------------
! This file contains the following two routines:
!
!  * fileio_output(filename,ix,jx,kx,t,x,y,z,U,V)
!  * fileio_input(filename,ix,jx,kx,t,x,y,z,U)
!
! The fileio_output routine outputs primitive variables for analysis.
! One can comment out the V-related part, in order to save disk space.
!-----------------------------------------------------------------------

! ----------------------------------------------------------------------
!     Output routine
! ----------------------------------------------------------------------
subroutine fileio_output(filename,ix,jx,kx,t,x,y,z,U,V)
  implicit none
  include 'param.h'
  integer, intent(in) :: ix, jx, kx
  real(8), intent(in) :: x(ix), y(jx), z(kx)
  real(8), intent(in) :: U(ix,jx,kx,var1)
  real(8), intent(in) :: V(ix,jx,kx,var2)
  real(8), intent(in) :: t
  integer :: i, j, k
  character*256 :: filename      

  open(16,file=filename,form='unformatted',access='stream')
  write(16) t
  write(16) ix
  write(16) jx
  write(16) kx
  write(16) (x(i),i=1,ix)
  write(16) (y(j),j=1,jx)
  write(16) (z(k),k=1,kx)
! -------------------------------------
  write(16) (((U(i,j,k,mx),i=1,ix),j=1,jx),k=1,kx)
  write(16) (((U(i,j,k,my),i=1,ix),j=1,jx),k=1,kx)
  write(16) (((U(i,j,k,mz),i=1,ix),j=1,jx),k=1,kx)
  write(16) (((U(i,j,k,en),i=1,ix),j=1,jx),k=1,kx)
  write(16) (((U(i,j,k,ro),i=1,ix),j=1,jx),k=1,kx)
  write(16) (((U(i,j,k,bx),i=1,ix),j=1,jx),k=1,kx)
  write(16) (((U(i,j,k,by),i=1,ix),j=1,jx),k=1,kx)
  write(16) (((U(i,j,k,bz),i=1,ix),j=1,jx),k=1,kx)
  write(16) (((U(i,j,k,ps),i=1,ix),j=1,jx),k=1,kx)
! -------------------------------------
  write(16) (((V(i,j,k,vx),i=1,ix),j=1,jx),k=1,kx)
  write(16) (((V(i,j,k,vy),i=1,ix),j=1,jx),k=1,kx)
  write(16) (((V(i,j,k,vz),i=1,ix),j=1,jx),k=1,kx)
  write(16) (((V(i,j,k,pr),i=1,ix),j=1,jx),k=1,kx)
! -------------------------------------
  close(16)

  return
end subroutine fileio_output


! ----------------------------------------------------------------------
!     Input routine
! ----------------------------------------------------------------------
subroutine fileio_input(filename,ix,jx,kx,t,x,y,z,U)
  implicit none
  include 'param.h'
  character*256 :: filename      
  integer, intent(in) :: ix, jx, kx
  real(8), intent(out) :: t
  real(8), intent(out) :: x(ix), y(jx), z(kx)
  real(8), intent(out) :: U(ix,jx,kx,var1)
!  real(8), intent(out) :: V(ix,jx,kx,var2)
  integer :: i, j, k, ix0, jx0, kx0
!  real(8) :: f1, v2, B2

  U(:,:,:,:) = 0.d0
!  V(:,:,:,:) = 0.d0

  open(15,file=filename,form='unformatted',access='stream')
  read(15) t
  read(15) ix0
  read(15) jx0
  read(15) kx0
  if(( ix0 /= ix ) .or. ( jx0 /= jx ) .or. ( kx0 /= kx )) then
     write(6,*) 'parameter mismatch'
     write(6,*) ' ix= ',ix,' ix0= ',ix0
     write(6,*) ' jx= ',jx,' jx0= ',jx0
     write(6,*) ' kx= ',kx,' kx0= ',kx0
     stop
  endif
  read(15) (x(i),i=1,ix)
  read(15) (y(j),j=1,jx)
  read(15) (z(k),k=1,kx)
  read(15) (((U(i,j,k,mx),i=1,ix),j=1,jx),k=1,kx)
  read(15) (((U(i,j,k,my),i=1,ix),j=1,jx),k=1,kx)
  read(15) (((U(i,j,k,mz),i=1,ix),j=1,jx),k=1,kx)
  read(15) (((U(i,j,k,en),i=1,ix),j=1,jx),k=1,kx)
  read(15) (((U(i,j,k,ro),i=1,ix),j=1,jx),k=1,kx)
  read(15) (((U(i,j,k,bx),i=1,ix),j=1,jx),k=1,kx)
  read(15) (((U(i,j,k,by),i=1,ix),j=1,jx),k=1,kx)
  read(15) (((U(i,j,k,bz),i=1,ix),j=1,jx),k=1,kx)
  read(15) (((U(i,j,k,ps),i=1,ix),j=1,jx),k=1,kx)
!  read(15) (((V(i,j,k,vx),i=1,ix),j=1,jx),k=1,kx)
!  read(15) (((V(i,j,k,vy),i=1,ix),j=1,jx),k=1,kx)
!  read(15) (((V(i,j,k,vz),i=1,ix),j=1,jx),k=1,kx)
!  read(15) (((V(i,j,k,pr),i=1,ix),j=1,jx),k=1,kx)
  close(15)

! ----------------------------------------------------------------------
! reconstructing V --> U
! I don't use the v2u routine due to the different martirx size
!  f1 = 1.d0 / ( gamma - 1 )
!  do k=1,kx
!  do j=1,jx
!  do i=1,ix
!     v2 = dot_product( V(i,j,k,vx:vz), V(i,j,k,vx:vz) )
!     B2 = dot_product( U(i,j,k,bx:bz), U(i,j,k,bx:bz) )
!     U(i,j,k,mx:mz) = U(i,j,k,ro) * V(i,j,k,vx:vz)
!     U(i,j,k,en)    = 0.5d0 * ( U(i,j,k,ro)*v2 + B2 ) + f1*V(i,j,k,pr)
!  enddo
!  enddo
!  enddo
! ----------------------------------------------------------------------

  return
end subroutine fileio_input

! ----------------------------------------------------------------------
