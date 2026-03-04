!-----------------------------------------------------------------------
!    MPI file IO routines
!-----------------------------------------------------------------------
!     2015/04/06  S. Zenitani  MPI-IO
!     2017/04/14  S. Zenitani  no longer use record markers
!     2018/05/02  S. Zenitani  2-D decomposition
!     2023/12/24  S. Zenitani  3-D decomposition
!     2026/03/04  S. Zenitani  rewritten for mpi_f08
!-----------------------------------------------------------------------
! This file contains the following two routines:
!
!  * mpiio_output(filename,ix,jx,kx,t,x,y,z,U,V)
!  * mpiio_input(filename,ix,jx,kx,t,x,y,z,U)
!
! The mpiio_output routine outputs primitive variables for analysis.
! One can comment out the V-related part, in order to save disk space.
! These routines depend on the parallel module.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     MPI-IO file output routine
!-----------------------------------------------------------------------
subroutine mpiio_output(filename,ix,jx,kx,t,x,y,z,U,V)
  use parallel
  implicit none
  include 'param.h'
  character*256 :: filename
  integer, intent(in) :: ix, jx, kx
  real(8), intent(in) :: x(ix), y(jx), z(kx)
  real(8), intent(in) :: U(ix,jx,kx,var1)
  real(8), intent(in) :: V(ix,jx,kx,var2)
  real(8), intent(in) :: t
  integer :: iix, is, ie
  integer :: jjx, js, je
  integer :: kkx, ks, ke
  type(mpi_file) :: fh
  type(mpi_datatype) :: ftype1, ftype2
  integer :: gsize(1), subsize(1), start(1)
  integer :: gsizes(3), subsizes(3), starts(3)
  integer(kind=mpi_offset_kind) :: disp, mk, mkiix, mkjjx, mkkkx

  iix = cart3d%sizes(1)*(ix-2) + 2
  jjx = cart3d%sizes(2)*(jx-2) + 2
  kkx = cart3d%sizes(3)*(kx-2) + 2

  call mpi_file_open(cart3d%comm, filename, & 
       mpi_mode_wronly + mpi_mode_create, & 
       mpi_info_null, fh )

  disp = 0
  call mpi_file_set_size(fh,disp)
  if( ranks%myrank==0 ) then
     call mpi_file_write(fh,  t, 1,   mpi_real8, mpi_status_ignore)
     call mpi_file_write(fh,iix, 1, mpi_integer, mpi_status_ignore)
     call mpi_file_write(fh,jjx, 1, mpi_integer, mpi_status_ignore)
     call mpi_file_write(fh,kkx, 1, mpi_integer, mpi_status_ignore)
  endif
  disp = disp + 8 + 4 + 4 + 4

  is = 2
  ie = ix-1
  if( cart3d%coords(1) == 0                 )  is = 1
  if( cart3d%coords(1) == cart3d%sizes(1)-1 )  ie = ix
  js = 2
  je = jx-1
  if( cart3d%coords(2) == 0                 )  js = 1
  if( cart3d%coords(2) == cart3d%sizes(2)-1 )  je = jx
  ks = 2
  ke = kx-1
  if( cart3d%coords(3) == 0                 )  ks = 1
  if( cart3d%coords(3) == cart3d%sizes(3)-1 )  ke = kx

! --------- 1D array of iix -------------------------------------------------
  gsize(1)   = iix
  subsize(1) = ie-is+1
  start(1)   = cart3d%coords(1)*(ix-2)+(is-1)
  call mpi_type_create_subarray(1,gsize,subsize,start,mpi_order_fortran,mpi_real8,ftype1)
  call mpi_type_commit(ftype1)

  gsize(1)   = ix
  subsize(1) = ie-is+1
  start(1)   = is-1
  call mpi_type_create_subarray(1,gsize,subsize,start,mpi_order_fortran,mpi_real8,ftype2)
  call mpi_type_commit(ftype2)

  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null)
  if( cart3d%coords(2) == 0 .and. cart3d%coords(3) == 0 ) then
     call mpi_file_write( fh, x, 1, ftype2, mpi_status_ignore )
  endif
! call mpi_file_write_all( fh, x, 1, ftype2, mpi_status_ignore )
  call mpi_type_free(ftype1)
  call mpi_type_free(ftype2)
  disp = disp+8*iix

  call mpi_barrier( cart3d%comm )

! --------- 1D array of jjx -------------------------------------------------
  gsize(1)   = jjx
  subsize(1) = je-js+1
  start(1)   = cart3d%coords(2)*(jx-2)+(js-1)
  call mpi_type_create_subarray(1,gsize,subsize,start,mpi_order_fortran,mpi_real8,ftype1)
  call mpi_type_commit(ftype1)

  gsize(1)   = jx
  subsize(1) = je-js+1
  start(1)   = js-1
  call mpi_type_create_subarray(1,gsize,subsize,start,mpi_order_fortran,mpi_real8,ftype2)
  call mpi_type_commit(ftype2)

  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null)
  if( cart3d%coords(1) == 0 .and. cart3d%coords(3) == 0 ) then
     call mpi_file_write(fh, y, 1, ftype2, mpi_status_ignore)
  endif
! call mpi_file_write_all( fh, y, 1, ftype2, mpi_status_ignore )
  call mpi_type_free(ftype1)
  call mpi_type_free(ftype2)
  disp = disp+8*jjx

  call mpi_barrier( cart3d%comm )

! --------- 1D array of kkx -------------------------------------------------
  gsize(1)   = kkx
  subsize(1) = ke-ks+1
  start(1)   = cart3d%coords(3)*(kx-2)+(ks-1)
  call mpi_type_create_subarray(1,gsize,subsize,start,mpi_order_fortran,mpi_real8,ftype1)
  call mpi_type_commit(ftype1)

  gsize(1)   = kx
  subsize(1) = ke-ks+1
  start(1)   = ks-1
  call mpi_type_create_subarray(1,gsize,subsize,start,mpi_order_fortran,mpi_real8,ftype2)
  call mpi_type_commit(ftype2)

  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null)
  if( cart3d%coords(1) == 0 .and. cart3d%coords(2) == 0 ) then
     call mpi_file_write( fh, z, 1, ftype2, mpi_status_ignore)
  endif
! call mpi_file_write_all( fh, z, 1, ftype2, mpi_status_ignore)
  call mpi_type_free(ftype1)
  call mpi_type_free(ftype2)
  disp = disp+8*jjx

  call mpi_barrier( cart3d%comm )

! --------- 3D matrix of iix * jjx * kkx -------------------------------------------------
  gsizes = (/iix,jjx,kkx/);  subsizes = (/(ie-is+1),(je-js+1),(ke-ks+1)/);
  starts = (/cart3d%coords(1)*(ix-2)+(is-1), cart3d%coords(2)*(jx-2)+(js-1), cart3d%coords(3)*(kx-2)+(ks-1)/)
  call mpi_type_create_subarray(3,gsizes,subsizes,starts,mpi_order_fortran,mpi_real8,ftype1)
  call mpi_type_commit(ftype1)
  gsizes = (/ix,jx,kx/);  subsizes = (/(ie-is+1),(je-js+1),(ke-ks+1)/);  starts = (/(is-1),(js-1),(ks-1)/)
  call mpi_type_create_subarray(3,gsizes,subsizes,starts,mpi_order_fortran,mpi_real8,ftype2)
  call mpi_type_commit(ftype2)

  mkiix = iix;  mkjjx = jjx;  mkkkx = kkx;  mk = 8*mkiix*mkjjx*mkkkx
! mk = 8*iix*jx
  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null)
  call mpi_file_write_all(fh, U(1,1,1,mx), 1, ftype2, mpi_status_ignore); disp = disp+mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null)
  call mpi_file_write_all(fh, U(1,1,1,my), 1, ftype2, mpi_status_ignore); disp = disp+mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null)
  call mpi_file_write_all(fh, U(1,1,1,mz), 1, ftype2, mpi_status_ignore); disp = disp+mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null)
  call mpi_file_write_all(fh, U(1,1,1,en), 1, ftype2, mpi_status_ignore); disp = disp+mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null)
  call mpi_file_write_all(fh, U(1,1,1,ro), 1, ftype2, mpi_status_ignore); disp = disp+mk

  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null)
  call mpi_file_write_all(fh, U(1,1,1,bx), 1, ftype2, mpi_status_ignore); disp = disp+mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null)
  call mpi_file_write_all(fh, U(1,1,1,by), 1, ftype2, mpi_status_ignore); disp = disp+mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null)
  call mpi_file_write_all(fh, U(1,1,1,bz), 1, ftype2, mpi_status_ignore); disp = disp+mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null)
  call mpi_file_write_all(fh, U(1,1,1,ps), 1, ftype2, mpi_status_ignore); disp = disp+mk
! -------------------------------------
  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null)
  call mpi_file_write_all(fh, V(1,1,1,vx), 1, ftype2, mpi_status_ignore); disp = disp+mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null)
  call mpi_file_write_all(fh, V(1,1,1,vy), 1, ftype2, mpi_status_ignore); disp = disp+mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null)
  call mpi_file_write_all(fh, V(1,1,1,vz), 1, ftype2, mpi_status_ignore); disp = disp+mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype1, "native", mpi_info_null)
  call mpi_file_write_all(fh, V(1,1,1,pr), 1, ftype2, mpi_status_ignore); disp = disp+mk
! -------------------------------------
  call mpi_type_free(ftype1)
  call mpi_type_free(ftype2)
! --------- 2D matrix of iix * jjx * kkx -------------------------------------------------

  call mpi_file_close( fh )

  return
end subroutine mpiio_output


!-----------------------------------------------------------------------
!     MPI-IO file input routine
!-----------------------------------------------------------------------
subroutine mpiio_input(filename,ix,jx,kx,t,x,y,z,U)
  use parallel
  implicit none
  include 'param.h'
  character*256 :: filename
  integer, intent(in) :: ix, jx, kx
  real(8), intent(out) :: t
  real(8), intent(out) :: x(ix), y(jx), z(kx)
  real(8), intent(out) :: U(ix,jx,kx,var1)
!  real(8), intent(out) :: V(ix,jx,kx,var2)
  integer :: iix, iix0, is, ie
  integer :: jjx, jjx0, js, je
  integer :: kkx, kkx0, ks, ke
  integer :: ijk
  type(mpi_file) :: fh
  type(mpi_datatype) :: ftype
  integer :: gsizes(3), subsizes(3), starts(3)
  integer(kind=mpi_offset_kind) :: disp, mk, mkiix, mkjjx, mkkkx


  iix = cart3d%sizes(1)*(ix-2) + 2
  jjx = cart3d%sizes(2)*(jx-2) + 2
  kkx = cart3d%sizes(3)*(kx-2) + 2

  call mpi_file_open(mpi_comm_world, filename, & 
       mpi_mode_rdonly, &
       mpi_info_null, fh )

  disp = 0
  call mpi_file_read_all(fh,  t, 1,   mpi_real8, mpi_status_ignore); disp = disp+8
  call mpi_file_read_all(fh,iix0,1, mpi_integer, mpi_status_ignore); disp = disp+4
  call mpi_file_read_all(fh,jjx0,1, mpi_integer, mpi_status_ignore); disp = disp+4
  call mpi_file_read_all(fh,kkx0,1, mpi_integer, mpi_status_ignore); disp = disp+4
  if(( iix0 /= iix ) .or. ( jjx0 /= jjx ) .or. ( kkx0 /= kkx )) then
     write(6,*) 'parameter mismatch'
     write(6,*) ' ix= ',iix,' ix0= ',iix0
     write(6,*) ' jx= ',jjx,' jx0= ',jjx0
     write(6,*) ' kx= ',kkx,' kx0= ',kkx0
     stop
  endif

  is = 2
  ie = ix-1
  if( cart3d%coords(1) == 0                 )  is = 1
  if( cart3d%coords(1) == cart3d%sizes(1)-1 )  ie = ix
  js = 2
  je = jx-1
  if( cart3d%coords(2) == 0                 )  js = 1
  if( cart3d%coords(2) == cart3d%sizes(2)-1 )  je = jx
  ks = 2
  ke = kx-1
  if( cart3d%coords(3) == 0                 )  ks = 1
  if( cart3d%coords(3) == cart3d%sizes(3)-1 )  ke = kx

  call mpi_file_set_view(fh, disp+8*cart3d%coords(1)*(ix-2), mpi_real8, mpi_real8, "native", mpi_info_null)
  call mpi_file_read_all(fh, x, ix, mpi_real8, mpi_status_ignore); disp = disp+8*iix
!  call mpi_file_set_view(fh, disp, mpi_byte, mpi_byte, "native", mpi_info_null)
  call mpi_file_set_view(fh, disp+8*cart3d%coords(2)*(jx-2), mpi_real8, mpi_real8, "native", mpi_info_null)
  call mpi_file_read_all(fh, y, jx,   mpi_real8, mpi_status_ignore); disp = disp+8*jjx
  call mpi_file_set_view(fh, disp+8*cart3d%coords(3)*(kx-2), mpi_real8, mpi_real8, "native", mpi_info_null)
  call mpi_file_read_all(fh, z, kx,   mpi_real8, mpi_status_ignore); disp = disp+8*kkx

! --------- 3D matrix of iix * jjx * kkx -------------------------------------------------
  gsizes = (/iix,jjx,kkx/);  subsizes = (/ix,jx,kx/);
  starts = (/0+cart3d%coords(1)*(ix-2), 0+cart3d%coords(2)*(jx-2), 0+cart3d%coords(3)*(kx-2)/)
  call mpi_type_create_subarray(3,gsizes,subsizes,starts,mpi_order_fortran,mpi_real8,ftype)
  call mpi_type_commit(ftype)

  mkiix = iix;  mkjjx = jjx;  mkkkx = kkx;  mk = 8*mkiix*mkjjx*mkkkx
  ijk = ix*jx*kx
  call mpi_file_set_view(fh, disp, mpi_real8, ftype, "native", mpi_info_null)
  call mpi_file_read_all(fh, U(1,1,1,mx), ijk, mpi_real8, mpi_status_ignore); disp = disp+mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype, "native", mpi_info_null)
  call mpi_file_read_all(fh, U(1,1,1,my), ijk, mpi_real8, mpi_status_ignore); disp = disp+mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype, "native", mpi_info_null)
  call mpi_file_read_all(fh, U(1,1,1,mz), ijk, mpi_real8, mpi_status_ignore); disp = disp+mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype, "native", mpi_info_null)
  call mpi_file_read_all(fh, U(1,1,1,en), ijk, mpi_real8, mpi_status_ignore); disp = disp+mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype, "native", mpi_info_null)
  call mpi_file_read_all(fh, U(1,1,1,ro), ijk, mpi_real8, mpi_status_ignore); disp = disp+mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype, "native", mpi_info_null)
  call mpi_file_read_all(fh, U(1,1,1,bx), ijk, mpi_real8, mpi_status_ignore); disp = disp+mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype, "native", mpi_info_null)
  call mpi_file_read_all(fh, U(1,1,1,by), ijk, mpi_real8, mpi_status_ignore); disp = disp+mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype, "native", mpi_info_null)
  call mpi_file_read_all(fh, U(1,1,1,bz), ijk, mpi_real8, mpi_status_ignore); disp = disp+mk
  call mpi_file_set_view(fh, disp, mpi_real8, ftype, "native", mpi_info_null)
  call mpi_file_read_all(fh, U(1,1,1,ps), ijk, mpi_real8, mpi_status_ignore); disp = disp+mk

!  call mpi_file_set_view(fh, disp, mpi_real8, ftype, "native", mpi_info_null)
!  call mpi_file_read_all(fh, V(1,1,1,vx), ijk, mpi_real8, mpi_status_ignore); disp = disp+mk
!  call mpi_file_set_view(fh, disp, mpi_real8, ftype, "native", mpi_info_null)
!  call mpi_file_read_all(fh, V(1,1,1,vy), ijk, mpi_real8, mpi_status_ignore); disp = disp+mk
!  call mpi_file_set_view(fh, disp, mpi_real8, ftype, "native", mpi_info_null)
!  call mpi_file_read_all(fh, V(1,1,1,vz), ijk, mpi_real8, mpi_status_ignore); disp = disp+mk
!  call mpi_file_set_view(fh, disp, mpi_real8, ftype, "native", mpi_info_null)
!  call mpi_file_read_all(fh, V(1,1,1,pr), ijk, mpi_real8, mpi_status_ignore); disp = disp+mk

  call mpi_type_free(ftype)
! --------- 2D matrix of iix * jjx -------------------------------------------------

  call mpi_file_close( fh )

  return
end subroutine mpiio_input

! ----------------------------------------------------------------------
