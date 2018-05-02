subroutine mpiinput(filename,ix,jx,t,x,y,U)
!-----------------------------------------------------------------------
!     MPI-IO file input routine
!-----------------------------------------------------------------------
!     2015/04/06  S. Zenitani  MPI-IO
!     2017/04/14  S. Zenitani  no longer use record markers
!     2018/05/02  S. Zenitani  2-D decomposition
!-----------------------------------------------------------------------
  use parallel
  implicit none
  include 'param.h'
  character*256 :: filename
  integer, intent(in) :: ix, jx
  real(8), intent(out) :: t
  real(8), intent(out) :: x(ix), y(jx)
  real(8), intent(out) :: U(ix,jx,var1)
!  real(8), intent(out) :: V(ix,jx,var2)
  integer :: iix, iix0, is, ie
  integer :: jjx, jjx0, js, je
  integer :: fh, ierr, ftype
  integer :: gsizes(2), subsizes(2), starts(2)
  integer(kind=mpi_offset_kind) :: disp, mk, mkiix, mkjjx


  iix = cart2d%sizes(1)*(ix-2) + 2
  jjx = cart2d%sizes(2)*(jx-2) + 2

  call mpi_file_open(mpi_comm_world, filename, & 
       mpi_mode_rdonly, &
       mpi_info_null, fh, ierr )

  disp = 0
  call mpi_file_read_all( fh,  t, 1,   mpi_real8, mpi_status_ignore, ierr); disp = disp+8
  call mpi_file_read_all( fh,iix0,1, mpi_integer, mpi_status_ignore, ierr); disp = disp+4
  call mpi_file_read_all( fh,jjx0,1, mpi_integer, mpi_status_ignore, ierr); disp = disp+4
  if(( iix0 /= iix ) .or. ( jjx0 /= jjx )) then
     write(6,*) 'parameter mismatch'
     write(6,*) ' ix= ',iix,' ix0= ',iix0
     write(6,*) ' jx= ',jjx,' jx0= ',jjx0
     stop
  endif

  is = 2
  ie = ix-1
  if( cart2d%coords(1) == 0                 )  is = 1
  if( cart2d%coords(1) == cart2d%sizes(1)-1 )  ie = ix
  js = 2
  je = jx-1
  if( cart2d%coords(2) == 0                 )  js = 1
  if( cart2d%coords(2) == cart2d%sizes(2)-1 )  je = jx

  call mpi_file_set_view(fh, disp+8*cart2d%coords(1)*(ix-2), mpi_real8, mpi_real8, "native", mpi_info_null, ierr)
  call mpi_file_read_all( fh, x, ix, mpi_real8, mpi_status_ignore, ierr); disp = disp+8*iix

!  call mpi_file_set_view(fh, disp, mpi_byte, mpi_byte, "native", mpi_info_null, ierr)
  call mpi_file_set_view(fh, disp+8*cart2d%coords(2)*(jx-2), mpi_real8, mpi_real8, "native", mpi_info_null, ierr)
  call mpi_file_read_all( fh, y, jx,   mpi_real8, mpi_status_ignore, ierr); disp = disp+8*jjx

! --------- 2D matrix of iix * jjx -------------------------------------------------
  gsizes = (/iix,jjx/);  subsizes = (/ix,jx/);
  starts = (/0+cart2d%coords(1)*(ix-2), 0+cart2d%coords(2)*(jx-2)/)
  call mpi_type_create_subarray(2,gsizes,subsizes,starts,mpi_order_fortran,mpi_real8,ftype,ierr)
  call mpi_type_commit( ftype,ierr )

  mkiix = iix
  mkjjx = jjx
  mk = 8*mkiix*mkjjx
  call mpi_file_set_view(fh, disp, mpi_real8, ftype, "native", mpi_info_null, ierr)
  call mpi_file_read_all( fh, U(1,1,mx), ix*jx, mpi_real8, mpi_status_ignore, ierr); disp = disp+mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype, "native", mpi_info_null, ierr)
  call mpi_file_read_all( fh, U(1,1,my), ix*jx, mpi_real8, mpi_status_ignore, ierr); disp = disp+mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype, "native", mpi_info_null, ierr)
  call mpi_file_read_all( fh, U(1,1,mz), ix*jx, mpi_real8, mpi_status_ignore, ierr); disp = disp+mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype, "native", mpi_info_null, ierr)
  call mpi_file_read_all( fh, U(1,1,en), ix*jx, mpi_real8, mpi_status_ignore, ierr); disp = disp+mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype, "native", mpi_info_null, ierr)
  call mpi_file_read_all( fh, U(1,1,ro), ix*jx, mpi_real8, mpi_status_ignore, ierr); disp = disp+mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype, "native", mpi_info_null, ierr)
  call mpi_file_read_all( fh, U(1,1,bx), ix*jx, mpi_real8, mpi_status_ignore, ierr); disp = disp+mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype, "native", mpi_info_null, ierr)
  call mpi_file_read_all( fh, U(1,1,by), ix*jx, mpi_real8, mpi_status_ignore, ierr); disp = disp+mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype, "native", mpi_info_null, ierr)
  call mpi_file_read_all( fh, U(1,1,bz), ix*jx, mpi_real8, mpi_status_ignore, ierr); disp = disp+mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype, "native", mpi_info_null, ierr)
  call mpi_file_read_all( fh, U(1,1,ps), ix*jx, mpi_real8, mpi_status_ignore, ierr); disp = disp+mk

!  call mpi_file_set_view(fh, disp, mpi_real8, ftype, "native", mpi_info_null, ierr)
!  call mpi_file_read_all( fh, V(1,1,vx), ix*jx, mpi_real8, mpi_status_ignore, ierr); disp = disp+mk
!  call mpi_file_set_view(fh, disp, mpi_real8, ftype, "native", mpi_info_null, ierr)
!  call mpi_file_read_all( fh, V(1,1,vy), ix*jx, mpi_real8, mpi_status_ignore, ierr); disp = disp+mk
!  call mpi_file_set_view(fh, disp, mpi_real8, ftype, "native", mpi_info_null, ierr)
!  call mpi_file_read_all( fh, V(1,1,vz), ix*jx, mpi_real8, mpi_status_ignore, ierr); disp = disp+mk
!  call mpi_file_set_view(fh, disp, mpi_real8, ftype, "native", mpi_info_null, ierr)
!  call mpi_file_read_all( fh, V(1,1,pr), ix*jx, mpi_real8, mpi_status_ignore, ierr); disp = disp+mk
! --------- 2D matrix of iix * jx -------------------------------------------------

  call mpi_file_close( fh, ierr )

  return
end subroutine mpiinput
