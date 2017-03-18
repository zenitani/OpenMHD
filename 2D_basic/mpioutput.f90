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
  integer :: mk
  integer :: iix
  integer :: fh, ierr, ftype
  integer :: gsize(2), subsize(2), start(2)
  integer(kind=mpi_offset_kind) :: disp
!-----------------------------------------------------------------------
! record marker
  integer, parameter :: mk_type = mpi_integer
  integer, parameter :: mk_size = 4
!  integer, parameter :: mk_type = mpi_integer8
!  integer, parameter :: mk_size = 8
!-----------------------------------------------------------------------

  iix = npe*(ix-2)+2
  call mpi_file_open(mpi_comm_world, filename, & 
       mpi_mode_wronly + mpi_mode_create, & 
       mpi_info_null, fh, ierr )

  disp = 0
  call mpi_file_set_size( fh,disp,ierr)

  mk = 8
  call mpi_file_write_all( fh, mk, 1, mk_type, mpi_status_ignore, ierr) ; disp = disp+mk_size
  call mpi_file_write_all( fh,  t, 1, mpi_real8, mpi_status_ignore, ierr) ; disp = disp+mk
  call mpi_file_write_all( fh, mk, 1, mk_type, mpi_status_ignore, ierr) ; disp = disp+mk_size
  mk = 4
  call mpi_file_write_all( fh, mk, 1, mk_type, mpi_status_ignore, ierr) ; disp = disp+mk_size
  call mpi_file_write_all( fh,iix, 1, mpi_integer, mpi_status_ignore, ierr) ; disp = disp+mk
  call mpi_file_write_all( fh, mk, 1, mk_type, mpi_status_ignore, ierr) ; disp = disp+mk_size
  call mpi_file_write_all( fh, mk, 1, mk_type, mpi_status_ignore, ierr) ; disp = disp+mk_size
  call mpi_file_write_all( fh, jx, 1, mpi_integer, mpi_status_ignore, ierr) ; disp = disp+mk
  call mpi_file_write_all( fh, mk, 1, mk_type, mpi_status_ignore, ierr) ; disp = disp+mk_size

  mk = 8*iix
  call mpi_file_write_all( fh, mk, 1, mk_type, mpi_status_ignore, ierr) ; disp = disp+mk_size
  call mpi_file_set_view(fh, disp+8*myrank*(ix-2), mpi_real8, mpi_real8, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, x, ix, mpi_real8, mpi_status_ignore, ierr) ; disp = disp+mk
  call mpi_file_set_view(fh, disp, mpi_byte, mpi_byte, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, mk, 1, mk_type, mpi_status_ignore, ierr) ; disp = disp+mk_size

  mk = 8*jx
  call mpi_file_write_all( fh, mk, 1, mk_type, mpi_status_ignore, ierr) ; disp = disp+mk_size
  call mpi_file_write_all( fh, y, jx, mpi_real8, mpi_status_ignore, ierr) ; disp = disp+mk
  call mpi_file_write_all( fh, mk, 1, mk_type, mpi_status_ignore, ierr) ; disp = disp+mk_size

! --------- 2D matrix of iix * jx -------------------------------------------------
  gsize = (/iix,jx/) ;  subsize = (/ix,jx/) ;  start = (/0+myrank*(ix-2), 0/)
  call mpi_type_create_subarray(2,gsize,subsize,start,mpi_order_fortran,mpi_real8,ftype,ierr)
  call mpi_type_commit( ftype,ierr )
  mk = 8*iix*jx

  call mpi_file_set_view(fh, disp, mpi_byte, mpi_byte, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, mk, 1, mk_type, mpi_status_ignore, ierr) ; disp = disp+mk_size
  call mpi_file_set_view(fh, disp, mpi_real8, ftype, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, U(1,1,mx), ix*jx, mpi_real8, mpi_status_ignore, ierr) ; disp = disp + mk
  call mpi_file_set_view(fh, disp, mpi_byte, mpi_byte, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, mk, 1, mk_type, mpi_status_ignore, ierr) ; disp = disp+mk_size

  call mpi_file_set_view(fh, disp, mpi_byte, mpi_byte, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, mk, 1, mk_type, mpi_status_ignore, ierr) ; disp = disp+mk_size
  call mpi_file_set_view(fh, disp, mpi_real8, ftype, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, U(1,1,my), ix*jx, mpi_real8, mpi_status_ignore, ierr) ; disp = disp + mk
  call mpi_file_set_view(fh, disp, mpi_byte, mpi_byte, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, mk, 1, mk_type, mpi_status_ignore, ierr) ; disp = disp+mk_size

  call mpi_file_set_view(fh, disp, mpi_byte, mpi_byte, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, mk, 1, mk_type, mpi_status_ignore, ierr) ; disp = disp+mk_size
  call mpi_file_set_view(fh, disp, mpi_real8, ftype, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, U(1,1,mz), ix*jx, mpi_real8, mpi_status_ignore, ierr) ; disp = disp + mk
  call mpi_file_set_view(fh, disp, mpi_byte, mpi_byte, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, mk, 1, mk_type, mpi_status_ignore, ierr) ; disp = disp+mk_size

  call mpi_file_write_all( fh, mk, 1, mk_type, mpi_status_ignore, ierr) ; disp = disp+mk_size
  call mpi_file_set_view(fh, disp, mpi_real8, ftype, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, U(1,1,en), ix*jx, mpi_real8, mpi_status_ignore, ierr) ; disp = disp + mk
  call mpi_file_set_view(fh, disp, mpi_byte, mpi_byte, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, mk, 1, mk_type, mpi_status_ignore, ierr) ; disp = disp+mk_size

  call mpi_file_write_all( fh, mk, 1, mk_type, mpi_status_ignore, ierr) ; disp = disp+mk_size
  call mpi_file_set_view(fh, disp, mpi_real8, ftype, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, U(1,1,ro), ix*jx, mpi_real8, mpi_status_ignore, ierr) ; disp = disp + mk
  call mpi_file_set_view(fh, disp, mpi_byte, mpi_byte, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, mk, 1, mk_type, mpi_status_ignore, ierr) ; disp = disp+mk_size

  call mpi_file_write_all( fh, mk, 1, mk_type, mpi_status_ignore, ierr) ; disp = disp+mk_size
  call mpi_file_set_view(fh, disp, mpi_real8, ftype, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, U(1,1,bx), ix*jx, mpi_real8, mpi_status_ignore, ierr) ; disp = disp + mk
  call mpi_file_set_view(fh, disp, mpi_byte, mpi_byte, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, mk, 1, mk_type, mpi_status_ignore, ierr) ; disp = disp+mk_size

  call mpi_file_write_all( fh, mk, 1, mk_type, mpi_status_ignore, ierr) ; disp = disp+mk_size
  call mpi_file_set_view(fh, disp, mpi_real8, ftype, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, U(1,1,by), ix*jx, mpi_real8, mpi_status_ignore, ierr) ; disp = disp + mk
  call mpi_file_set_view(fh, disp, mpi_byte, mpi_byte, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, mk, 1, mk_type, mpi_status_ignore, ierr) ; disp = disp+mk_size

  call mpi_file_write_all( fh, mk, 1, mk_type, mpi_status_ignore, ierr) ; disp = disp+mk_size
  call mpi_file_set_view(fh, disp, mpi_real8, ftype, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, U(1,1,bz), ix*jx, mpi_real8, mpi_status_ignore, ierr) ; disp = disp + mk
  call mpi_file_set_view(fh, disp, mpi_byte, mpi_byte, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, mk, 1, mk_type, mpi_status_ignore, ierr) ; disp = disp+mk_size

  call mpi_file_write_all( fh, mk, 1, mk_type, mpi_status_ignore, ierr) ; disp = disp+mk_size
  call mpi_file_set_view(fh, disp, mpi_real8, ftype, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, U(1,1,ps), ix*jx, mpi_real8, mpi_status_ignore, ierr) ; disp = disp + mk
  call mpi_file_set_view(fh, disp, mpi_byte, mpi_byte, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, mk, 1, mk_type, mpi_status_ignore, ierr) ; disp = disp+mk_size

  call mpi_file_write_all( fh, mk, 1, mk_type, mpi_status_ignore, ierr) ; disp = disp+mk_size
  call mpi_file_set_view(fh, disp, mpi_real8, ftype, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, V(1,1,vx), ix*jx, mpi_real8, mpi_status_ignore, ierr) ; disp = disp + mk
  call mpi_file_set_view(fh, disp, mpi_byte, mpi_byte, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, mk, 1, mk_type, mpi_status_ignore, ierr) ; disp = disp+mk_size

  call mpi_file_write_all( fh, mk, 1, mk_type, mpi_status_ignore, ierr) ; disp = disp+mk_size
  call mpi_file_set_view(fh, disp, mpi_real8, ftype, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, V(1,1,vy), ix*jx, mpi_real8, mpi_status_ignore, ierr) ; disp = disp + mk
  call mpi_file_set_view(fh, disp, mpi_byte, mpi_byte, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, mk, 1, mk_type, mpi_status_ignore, ierr) ; disp = disp+mk_size

  call mpi_file_write_all( fh, mk, 1, mk_type, mpi_status_ignore, ierr) ; disp = disp+mk_size
  call mpi_file_set_view(fh, disp, mpi_real8, ftype, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, V(1,1,vz), ix*jx, mpi_real8, mpi_status_ignore, ierr) ; disp = disp + mk
  call mpi_file_set_view(fh, disp, mpi_byte, mpi_byte, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, mk, 1, mk_type, mpi_status_ignore, ierr) ; disp = disp+mk_size

  call mpi_file_write_all( fh, mk, 1, mk_type, mpi_status_ignore, ierr) ; disp = disp+mk_size
  call mpi_file_set_view(fh, disp, mpi_real8, ftype, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, V(1,1,pr), ix*jx, mpi_real8, mpi_status_ignore, ierr) ; disp = disp + mk
  call mpi_file_set_view(fh, disp, mpi_byte, mpi_byte, "native", mpi_info_null, ierr)
  call mpi_file_write_all( fh, mk, 1, mk_type, mpi_status_ignore, ierr) ; disp = disp+mk_size
! --------- 2D matrix of iix * jx -------------------------------------------------

  call mpi_file_close( fh, ierr )

  return
end subroutine mpioutput
