subroutine mpioutput(filename,ix,jx,t,x,y,U,V,myrank,npe)
  implicit none
  include 'param.h'
  include 'mpif.h'
  character*256 :: filename
  integer, intent(in) :: ix, jx, myrank, npe
  real(8), intent(in) :: x(ix), y(jx)
  real(8), intent(in) :: U(ix,jx,var1)
  real(8), intent(in) :: V(ix,jx,var2)
  real(8), intent(in) :: t
  integer :: iix, is, ie
  integer :: fh, ierr, ftype1, ftype2
  integer :: gsize(2), subsize(2), start(2)
  integer(kind=mpi_offset_kind) :: disp, mk, mkiix, mkjjx

  iix = npe*(ix-2)+2
  call mpi_file_open(mpi_comm_world, filename, & 
       mpi_mode_wronly + mpi_mode_create, & 
       mpi_info_null, fh, ierr )

  disp = 0
  call mpi_file_set_size( fh,disp,ierr)
  if( myrank==0 ) then
     call mpi_file_write( fh,  t, 1,   mpi_real8, mpi_status_ignore, ierr)
     call mpi_file_write( fh,iix, 1, mpi_integer, mpi_status_ignore, ierr)
     call mpi_file_write( fh, jx, 1, mpi_integer, mpi_status_ignore, ierr)
  endif
  disp = disp + 8 + 4 + 4

  if( myrank==0 ) then
     is = 1
     ie = ix-1
  elseif( myrank == (npe-1) ) then
     is = 2
     ie = ix
  else
     is = 2
     ie = ix-1
  endif

! --------- 1D matrix of iix -------------------------------------------------
  call mpi_type_create_subarray(1,iix,(ie-is+1),myrank*(ix-2)+(is-1),mpi_order_fortran,mpi_real8,ftype1,ierr)
  call mpi_type_commit( ftype1,ierr )
  call mpi_type_create_subarray(1,ix,(ie-is+1),(is-1),mpi_order_fortran,mpi_real8,ftype2,ierr)
  call mpi_type_commit( ftype2,ierr )

  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, x, 1, ftype2, mpi_status_ignore, ierr)
  disp = disp+8*iix

  call mpi_file_set_view(fh, disp, mpi_byte, mpi_byte, "native", mpi_info_null, ierr)
  if( myrank==0 ) then
    call mpi_file_write( fh, y, jx, mpi_real8, mpi_status_ignore, ierr)
  endif
  disp = disp+8*jx

! --------- 2D matrix of iix * jx -------------------------------------------------
  gsize = (/iix,jx/) ;  subsize = (/(ie-is+1),jx/) ;  start = (/myrank*(ix-2)+(is-1), 0/)
  call mpi_type_create_subarray(2,gsize,subsize,start,mpi_order_fortran,mpi_real8,ftype1,ierr)
  call mpi_type_commit( ftype1,ierr )
  gsize = (/ix,jx/) ;  subsize = (/(ie-is+1),jx/) ;  start = (/(is-1), 0/)
  call mpi_type_create_subarray(2,gsize,subsize,start,mpi_order_fortran,mpi_real8,ftype2,ierr)
  call mpi_type_commit( ftype2,ierr )

  mkiix = iix
  mkjjx = jx
  mk = 8*mkiix*mkjjx
! mk = 8*iix*jx
  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, U(1,1,mx), 1, ftype2, mpi_status_ignore, ierr) ; disp = disp + mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, U(1,1,my), 1, ftype2, mpi_status_ignore, ierr) ; disp = disp + mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, U(1,1,mz), 1, ftype2, mpi_status_ignore, ierr) ; disp = disp + mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, U(1,1,en), 1, ftype2, mpi_status_ignore, ierr) ; disp = disp + mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, U(1,1,ro), 1, ftype2, mpi_status_ignore, ierr) ; disp = disp + mk

  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, U(1,1,bx), 1, ftype2, mpi_status_ignore, ierr) ; disp = disp + mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, U(1,1,by), 1, ftype2, mpi_status_ignore, ierr) ; disp = disp + mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, U(1,1,bz), 1, ftype2, mpi_status_ignore, ierr) ; disp = disp + mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, U(1,1,ps), 1, ftype2, mpi_status_ignore, ierr) ; disp = disp + mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, V(1,1,vx), 1, ftype2, mpi_status_ignore, ierr) ; disp = disp + mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, V(1,1,vy), 1, ftype2, mpi_status_ignore, ierr) ; disp = disp + mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, V(1,1,vz), 1, ftype2, mpi_status_ignore, ierr) ; disp = disp + mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, V(1,1,pr), 1, ftype2, mpi_status_ignore, ierr) ; disp = disp + mk
! --------- 2D matrix of iix * jx -------------------------------------------------

  call mpi_file_close( fh, ierr )

  return
end subroutine mpioutput
