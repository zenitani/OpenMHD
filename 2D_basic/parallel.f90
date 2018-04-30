module parallel
  use mpi
  implicit none

  ! 2-D domain
  integer, parameter :: ndims = 2

  ! rank
  type myranks
     integer :: size = 0
     integer :: myrank
     integer :: north, east, south, west
  end type myranks

  type(myranks), save :: ranks         ! inter-node communication
  type(myranks), save :: ranks_local   ! local communication
  integer, save, private :: comm_local ! local communication

  ! topology
  type mycoords
     integer :: comm
     integer :: sizes(ndims)
     integer :: coords(ndims)
     logical :: periods(ndims)
  end type mycoords
  type(mycoords), save :: cart2d

  ! MPI-3 shared-memory communication
  logical, parameter, private :: use_shm = .false.


contains

  subroutine parallel_init(my_sizes,my_periods)
    use, intrinsic :: iso_c_binding, only : c_f_pointer, c_ptr
    include 'param.h'
    integer, intent(in) :: my_sizes(2)
    logical, intent(in) :: my_periods(2)
    integer :: merr
    integer :: group_world, group_local
    integer :: tmpA(4), tmpB(4)

    call mpi_init(merr)
!   call mpi_comm_size(mpi_comm_world, ranks%size, merr)
!   call mpi_comm_rank(mpi_comm_world, ranks%myrank, merr)

    cart2d%sizes(:) = my_sizes(:)
    cart2d%coords(:)  = 0
    cart2d%periods(:) = my_periods(:)
    call mpi_cart_create(mpi_comm_world, ndims, cart2d%sizes, cart2d%periods, .true., cart2d%comm, merr)
    call mpi_comm_size (cart2d%comm, ranks%size, merr)
    call mpi_comm_rank (cart2d%comm, ranks%myrank, merr)
    call mpi_cart_shift(cart2d%comm, 0, 1,ranks%west, ranks%east, merr)
    call mpi_cart_shift(cart2d%comm, 1, 1,ranks%south,ranks%north,merr)
    call mpi_cart_coords(cart2d%comm, ranks%myrank, 2, cart2d%coords, merr)

    if( ranks%size /= my_sizes(1)*my_sizes(2) ) then
       if( ranks%myrank == 0 ) then
          write(6,*) 'MPI process numbers mismatch.'
       endif
       stop
    endif

    tmpA = (/ranks%north, ranks%east, ranks%south, ranks%west/)
    call mpi_comm_split_type(cart2d%comm, mpi_comm_type_shared, 0, mpi_info_null, comm_local, merr)
    call mpi_comm_group(cart2d%comm, group_world, merr)
    call mpi_comm_group(comm_local, group_local, merr)

    call mpi_group_translate_ranks(group_world, 4, tmpA, group_local, tmpB, merr)
    call mpi_comm_size(comm_local, ranks_local%size, merr)
    ranks_local%north = tmpB(1) ;   ranks_local%east = tmpB(2)
    ranks_local%south = tmpB(3) ;   ranks_local%west = tmpB(4)

!    write(6,*) 'I am ', ranks%myrank, cart2d%coords(1), cart2d%coords(2)
!    write(6,*) 'My (', ranks%myrank, ') neighbours are ', ranks_local%north, ranks_local%east, ranks_local%south, ranks_local%west
    return
  end subroutine parallel_init


  subroutine parallel_finilize
    integer :: merr

    call mpi_finalize(merr)

    return
  end subroutine parallel_finilize


  subroutine parallel_exchange(U,ix,jx,dir)
    include 'param.h'
    integer, intent(in) :: ix, jx
    real(8), intent(inout) :: U(ix,jx,var1)
    ! direction [input]: 1 (X), 2 (Y), 3 (Z)
    integer, intent(in)  :: dir
!----------------------------------------------------------------------
    real(8), allocatable :: bufsnd1(:,:), bufrcv1(:,:)
    real(8), allocatable :: bufsnd2(:,:), bufrcv2(:,:)
    integer :: mreq1(2), mreq2(2)
    integer :: msize, merr
    integer :: mleft, mright
!----------------------------------------------------------------------

    select case(dir)
    case(1)

       allocate( bufsnd1(jx,var1), bufrcv1(jx,var1) )
       allocate( bufsnd2(jx,var1), bufrcv2(jx,var1) )
       mleft  = ranks%west
       mright = ranks%east
       msize  = jx*var1

       if( cart2d%sizes(1) == 1 ) then

          ! left/right
          if( mleft /= mpi_proc_null ) then
             U(1,:,:)  = U(ix-1,:,:)
             U(ix,:,:) = U(2,:,:)
          endif

       else

          if( .not. use_shm ) then

             ! nonblocking communication (mreq1)
             call mpi_irecv(bufrcv1,msize,mpi_real8,mright,0,cart2d%comm,mreq1(1),merr)
             bufsnd1(:,:) = U(2,:,:)
             call mpi_isend(bufsnd1,msize,mpi_real8,mleft,0,cart2d%comm,mreq1(2),merr)
             ! nonblocking communication (mreq2)
             call mpi_irecv(bufrcv2,msize,mpi_real8,mleft,1,cart2d%comm,mreq2(1),merr)
             bufsnd2(:,:) = U(ix-1,:,:)
             call mpi_isend(bufsnd2,msize,mpi_real8,mright,1,cart2d%comm,mreq2(2),merr)

             call mpi_waitall(2,mreq1,mpi_statuses_ignore,merr)
             if( mright /= mpi_proc_null ) then
                U(ix,:,:) = bufrcv1(:,:)
             endif
             call mpi_waitall(2,mreq2,mpi_statuses_ignore,merr)
             if( mleft /= mpi_proc_null ) then
                U(1,:,:) = bufrcv2(:,:)
             endif

          else
          
          endif

       endif
       deallocate( bufsnd1, bufrcv1 )
       deallocate( bufsnd2, bufrcv2 )


    case(2)

       allocate( bufsnd1(ix,var1), bufrcv1(ix,var1) )
       allocate( bufsnd2(ix,var1), bufrcv2(ix,var1) )
       mleft  = ranks%south
       mright = ranks%north
       msize  = ix*var1

       if( cart2d%sizes(2) == 1 ) then

          ! top/bottom boundaries
          if( mleft /= mpi_proc_null ) then
             U(:,jx,:) = U(:,2,:)
             U(:,1,:)  = U(:,jx-1,:)
          endif

       else

          if( .not. use_shm ) then

             ! nonblocking communication (mreq1)
             call mpi_irecv(bufrcv1,msize,mpi_real8,mright,0,cart2d%comm,mreq1(1),merr)
             bufsnd1(:,:) = U(:,2,:)
             call mpi_isend(bufsnd1,msize,mpi_real8,mleft,0,cart2d%comm,mreq1(2),merr)
             ! nonblocking communication (mreq2)
             call mpi_irecv(bufrcv2,msize,mpi_real8,mleft,1,cart2d%comm,mreq2(1),merr)
             bufsnd2(:,:) = U(:,jx-1,:)
             call mpi_isend(bufsnd2,msize,mpi_real8,mright,1,cart2d%comm,mreq2(2),merr)

             call mpi_waitall(2,mreq1,mpi_statuses_ignore,merr)
             if( mright /= mpi_proc_null ) then
                U(:,jx,:) = bufrcv1(:,:)
             endif
             call mpi_waitall(2,mreq2,mpi_statuses_ignore,merr)
             if( mleft /= mpi_proc_null ) then
                U(:,1,:) = bufrcv2(:,:)
             endif

          else

          endif

       endif
       deallocate( bufsnd1, bufrcv1 )
       deallocate( bufsnd2, bufrcv2 )
   
    endselect

  end subroutine parallel_exchange


  subroutine parallel_exchange2(VL,VR,ix,jx,dir)
    include 'param.h'
    integer, intent(in) :: ix, jx
    ! left flux (VL) [input]
    real(8), intent(inout) :: VL(ix,jx,var1)
    ! right flux (VR) [input]
    real(8), intent(inout) :: VR(ix,jx,var1)
    ! direction [input]: 1 (X), 2 (Y), 3 (Z)
    integer, intent(in)  :: dir
!----------------------------------------------------------------------
    real(8), allocatable :: bufsnd1(:,:), bufrcv1(:,:)
    real(8), allocatable :: bufsnd2(:,:), bufrcv2(:,:)
    integer :: mreq1(2), mreq2(2)
    integer :: msize, merr
    integer :: mleft, mright
    integer, parameter :: use_mpi3_shm = 0
!----------------------------------------------------------------------

    select case(dir)
    case(1)

       allocate( bufsnd1(jx,var1), bufrcv1(jx,var1) )
       allocate( bufsnd2(jx,var1), bufrcv2(jx,var1) )
       mleft  = ranks%west
       mright = ranks%east
       msize  = jx*var1

       if( cart2d%sizes(1) == 1 ) then

          ! left/right
          if( mleft /= mpi_proc_null ) then
             VR(ix-1,:,:) = VR(1,:,:)
             VL(1,:,:)    = VL(ix-1,:,:)
          endif

       else

          if( .not. use_shm ) then

             ! nonblocking communication (mreq1)
             call mpi_irecv(bufrcv1,msize,mpi_real8,mright,0,cart2d%comm,mreq1(1),merr)
             bufsnd1(:,:) = VR(1,:,:)
             call mpi_isend(bufsnd1,msize,mpi_real8,mleft,0,cart2d%comm,mreq1(2),merr)
             ! nonblocking communication (mreq2)
             call mpi_irecv(bufrcv2,msize,mpi_real8,mleft,1,cart2d%comm,mreq2(1),merr)
             bufsnd2(:,:) = VL(ix-1,:,:)
             call mpi_isend(bufsnd2,msize,mpi_real8,mright,1,cart2d%comm,mreq2(2),merr)

             call mpi_waitall(2,mreq1,mpi_statuses_ignore,merr)
             if( mright /= mpi_proc_null ) then
                VR(ix-1,:,:) = bufrcv1(:,:)
             endif
             call mpi_waitall(2,mreq2,mpi_statuses_ignore,merr)
             if( mleft /= mpi_proc_null ) then
                VL(1,:,:) = bufrcv2(:,:)
             endif

          else
          
          endif

       endif
       deallocate( bufsnd1, bufrcv1 )
       deallocate( bufsnd2, bufrcv2 )

    case(2)

       allocate( bufsnd1(ix,var1), bufrcv1(ix,var1) )
       allocate( bufsnd2(ix,var1), bufrcv2(ix,var1) )
       mleft  = ranks%south
       mright = ranks%north
       msize  = ix*var1

       if( cart2d%sizes(2) == 1 ) then

          ! top/bottom boundaries
          if( mleft /= mpi_proc_null ) then
             VR(:,jx-1,:) = VR(:,1,:)
             VL(:,   1,:) = VL(:,jx-1,:)
          endif

       else

          if( .not. use_shm ) then

             ! nonblocking communication (mreq1)
             call mpi_irecv(bufrcv1,msize,mpi_real8,mright,0,cart2d%comm,mreq1(1),merr)
             bufsnd1(:,:) = VR(:,1,:)
             call mpi_isend(bufsnd1,msize,mpi_real8,mleft,0,cart2d%comm,mreq1(2),merr)
             ! nonblocking communication (mreq2)
             call mpi_irecv(bufrcv2,msize,mpi_real8,mleft,1,cart2d%comm,mreq2(1),merr)
             bufsnd2(:,:) = VL(:,jx-1,:)
             call mpi_isend(bufsnd2,msize,mpi_real8,mright,1,cart2d%comm,mreq2(2),merr)

             call mpi_waitall(2,mreq1,mpi_statuses_ignore,merr)
             if( mright /= mpi_proc_null ) then
                VR(:,jx-1,:) = bufrcv1(:,:)
             endif
             call mpi_waitall(2,mreq2,mpi_statuses_ignore,merr)
             if( mleft /= mpi_proc_null ) then
                VL(:,1,:) = bufrcv2(:,:)
             endif

          else

          endif

       endif
       deallocate( bufsnd1, bufrcv1 )
       deallocate( bufsnd2, bufrcv2 )
   
    endselect

  end subroutine parallel_exchange2

end module parallel
