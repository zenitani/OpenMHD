subroutine mpioutput(filename,ix,jx,t,x,y,U,V)
!-----------------------------------------------------------------------
!     MPI-IO file output routine
!-----------------------------------------------------------------------
!     2015/04/06  S. Zenitani  MPI-IO
!     2017/04/14  S. Zenitani  no longer use record markers
!     2018/01/17  S. Zenitani  subarray for mpi_file_write_all
!-----------------------------------------------------------------------
!  use mpi
  use parallel
  implicit none
  include 'param.h'
  character*256 :: filename
  integer, intent(in) :: ix, jx
  real(8), intent(in) :: x(ix), y(jx)
  real(8), intent(in) :: U(ix,jx,var1)
  real(8), intent(in) :: V(ix,jx,var2)
  real(8), intent(in) :: t
  integer :: iix, is, ie
  integer :: jjx, js, je
  integer :: fh, ierr, ftype1, ftype2
  integer :: gsize(1), subsize(1), start(1)
  integer :: gsizes(2), subsizes(2), starts(2)
  integer(kind=mpi_offset_kind) :: disp, mk, mkiix, mkjjx

  iix = cart2d%sizes(1)*(ix-2) + 2
  jjx = cart2d%sizes(2)*(jx-2) + 2

  call mpi_file_open(cart2d%comm, filename, & 
       mpi_mode_wronly + mpi_mode_create, & 
       mpi_info_null, fh, ierr )

  disp = 0
  call mpi_file_set_size( fh,disp,ierr)
  if( ranks%myrank==0 ) then
     call mpi_file_write( fh,  t, 1,   mpi_real8, mpi_status_ignore, ierr)
     call mpi_file_write( fh,iix, 1, mpi_integer, mpi_status_ignore, ierr)
     call mpi_file_write( fh,jjx, 1, mpi_integer, mpi_status_ignore, ierr)
  endif
  disp = disp + 8 + 4 + 4

  is = 2
  ie = ix-1
  if( cart2d%coords(1) == 0                 )  is = 1
  if( cart2d%coords(1) == cart2d%sizes(1)-1 )  ie = ix
  js = 2
  je = jx-1
  if( cart2d%coords(2) == 0                 )  js = 1
  if( cart2d%coords(2) == cart2d%sizes(2)-1 )  je = jx

! --------- 1D matrix of iix -------------------------------------------------
!  if( cart2d%coords(1) == 0 ) then
!     write(6,*) 'Hello?'
!  endif
     gsize(1)   = iix
     subsize(1) = ie-is+1
     start(1)   = cart2d%coords(1)*(ix-2)+(is-1)
     call mpi_type_create_subarray(1,gsize,subsize,start,mpi_order_fortran,mpi_real8,ftype1,ierr)
     call mpi_type_commit( ftype1,ierr )

     gsize(1)   = ix
     subsize(1) = ie-is+1
     start(1)   = is-1
     call mpi_type_create_subarray(1,gsize,subsize,start,mpi_order_fortran,mpi_real8,ftype2,ierr)
     call mpi_type_commit( ftype2,ierr )

     call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null, ierr)
     call mpi_file_write( fh, x, 1, ftype2, mpi_status_ignore, ierr)
!     call mpi_file_write_all( fh, x, 1, ftype2, mpi_status_ignore, ierr)
!  endif
  disp = disp+8*iix

! --------- 1D matrix of jjx -------------------------------------------------
!  if( cart2d%coords(2) == 0 ) then
!     write(6,*) 'Hello?'
!  endif
     gsize(1)   = jjx
     subsize(1) = je-js+1
     start(1)   = cart2d%coords(2)*(jx-2)+(js-1)
     call mpi_type_create_subarray(1,gsize,subsize,start,mpi_order_fortran,mpi_real8,ftype1,ierr)
     call mpi_type_commit( ftype1,ierr )

     gsize(1)   = jx
     subsize(1) = je-js+1
     start(1)   = js-1
     call mpi_type_create_subarray(1,gsize,subsize,start,mpi_order_fortran,mpi_real8,ftype2,ierr)
     call mpi_type_commit( ftype2,ierr )

     call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null, ierr)
     call mpi_file_write( fh, y, 1, ftype2, mpi_status_ignore, ierr)
!     call mpi_file_write_all( fh, y, 1, ftype2, mpi_status_ignore, ierr)
!     call mpi_file_set_view(fh, disp, mpi_byte, mpi_byte, "native", mpi_info_null, ierr)
!     call mpi_file_write( fh, y, jx, mpi_real8, mpi_status_ignore, ierr)
!  endif
  disp = disp+8*jjx

! --------- 2D matrix of iix * jjx -------------------------------------------------
  gsizes = (/iix,jjx/);  subsizes = (/(ie-is+1),(je-js+1)/);
  starts = (/cart2d%coords(1)*(ix-2)+(is-1), cart2d%coords(2)*(jx-2)+(js-1)/)
  call mpi_type_create_subarray(2,gsizes,subsizes,starts,mpi_order_fortran,mpi_real8,ftype1,ierr)
  call mpi_type_commit( ftype1,ierr )
  gsizes = (/ix,jx/);  subsizes = (/(ie-is+1),(je-js+1)/);  starts = (/(is-1), (js-1)/)
  call mpi_type_create_subarray(2,gsizes,subsizes,starts,mpi_order_fortran,mpi_real8,ftype2,ierr)
  call mpi_type_commit( ftype2,ierr )

  mkiix = iix
  mkjjx = jjx
  mk = 8*mkiix*mkjjx
! mk = 8*iix*jx
  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, U(1,1,mx), 1, ftype2, mpi_status_ignore, ierr); disp = disp+mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, U(1,1,my), 1, ftype2, mpi_status_ignore, ierr); disp = disp+mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, U(1,1,mz), 1, ftype2, mpi_status_ignore, ierr); disp = disp+mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, U(1,1,en), 1, ftype2, mpi_status_ignore, ierr); disp = disp+mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, U(1,1,ro), 1, ftype2, mpi_status_ignore, ierr); disp = disp+mk

  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, U(1,1,bx), 1, ftype2, mpi_status_ignore, ierr); disp = disp+mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, U(1,1,by), 1, ftype2, mpi_status_ignore, ierr); disp = disp+mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, U(1,1,bz), 1, ftype2, mpi_status_ignore, ierr); disp = disp+mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, U(1,1,ps), 1, ftype2, mpi_status_ignore, ierr); disp = disp+mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, V(1,1,vx), 1, ftype2, mpi_status_ignore, ierr); disp = disp+mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, V(1,1,vy), 1, ftype2, mpi_status_ignore, ierr); disp = disp+mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, V(1,1,vz), 1, ftype2, mpi_status_ignore, ierr); disp = disp+mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, V(1,1,pr), 1, ftype2, mpi_status_ignore, ierr); disp = disp+mk
! --------- 2D matrix of iix * jx -------------------------------------------------

  call mpi_file_close( fh, ierr )

  return
end subroutine mpioutput
