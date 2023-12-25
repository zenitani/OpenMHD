module parallel
  use mpi
  implicit none

  ! 3-D domain
  integer, parameter :: ndims = 3

  ! topology
  type mycoords
     integer :: comm
     integer :: sizes(ndims)
     integer :: coords(ndims)
     logical :: periods(ndims)
  end type mycoords
  type(mycoords) :: cart3d

  ! rank info
  type myranks
     integer :: size = 0
     integer :: myrank
     integer :: north, east, south, west, up, down
  end type myranks

  type(myranks) :: ranks         ! global communication
  type(myranks) :: ranks_local   ! intra-node communication
  integer, private :: comm_local ! intra-node communication

  ! ----- MPI-3 shared memory communication ------------------------
  logical, parameter, private :: use_shm = .false.  ! MPI communication
! logical, parameter, private :: use_shm = .true.   ! MPI-3 SHM model (experimental)

  ! Please do not edit here
  integer, private :: mpi_mode(ndims) = (/1, 1, 1/)  ! 0: no MPI, 1: MPI-1, 3: MPI-3
  integer, private :: mwin1, mwin2, mwin3
  real(8), private, dimension(:,:,:,:), pointer, contiguous :: fptr1 => null()
  real(8), private, dimension(:,:,:,:), pointer, contiguous :: fwest => null()
  real(8), private, dimension(:,:,:,:), pointer, contiguous :: feast => null()
  real(8), private, dimension(:,:,:,:), pointer, contiguous :: fptr2  => null()
  real(8), private, dimension(:,:,:,:), pointer, contiguous :: fsouth => null()
  real(8), private, dimension(:,:,:,:), pointer, contiguous :: fnorth => null()
  real(8), private, dimension(:,:,:,:), pointer, contiguous :: fptr3  => null()
  real(8), private, dimension(:,:,:,:), pointer, contiguous :: fdown => null()
  real(8), private, dimension(:,:,:,:), pointer, contiguous :: fup => null()
  ! ----- MPI-3 shared memory communication ------------------------

contains

  subroutine parallel_init(my_sizes,my_periods,ix,jx,kx)
    use, intrinsic :: iso_c_binding, only : c_f_pointer, c_ptr
    include 'param.h'
    integer, intent(in) :: my_sizes(ndims)
    logical, intent(in) :: my_periods(ndims)
    integer, intent(in) :: ix, jx, kx
    integer :: i
    integer :: tmpA(6), tmpB(6)
    integer :: group_world, group_local
    integer :: mdisp, merr, merrcode
    integer(kind=mpi_address_kind) :: msize
    type(c_ptr) :: baseptr1, baseptr2, baseptr3

    call mpi_init(merr)

    cart3d%sizes(:) = my_sizes(:)
    cart3d%coords(:)  = 0
    cart3d%periods(:) = my_periods(:)

    ! After this, one should use cart3d%comm instead of mpi_comm_world
    call mpi_cart_create(mpi_comm_world, ndims, cart3d%sizes, cart3d%periods, .true., cart3d%comm, merr)
    call mpi_comm_size (cart3d%comm, ranks%size, merr)
    call mpi_comm_rank (cart3d%comm, ranks%myrank, merr)
    call mpi_cart_shift(cart3d%comm, 0, 1, ranks%west, ranks%east, merr)
    call mpi_cart_shift(cart3d%comm, 1, 1, ranks%south,ranks%north,merr)
    call mpi_cart_shift(cart3d%comm, 2, 1, ranks%down, ranks%up, merr)
    call mpi_cart_coords(cart3d%comm, ranks%myrank, 3, cart3d%coords, merr)

    if( ranks%size /= my_sizes(1)*my_sizes(2)*my_sizes(3) ) then
       if( ranks%myrank == 0 ) then
          write(6,*) 'MPI process numbers mismatch.'
       endif
       merrcode = -1
       call mpi_abort(cart3d%comm, merrcode, merr)
    endif

    ! MPI mode
    ! No MPI transport in the case of 1
    if( cart3d%sizes(1) == 1 )  mpi_mode(1) = 0
    if( cart3d%sizes(2) == 1 )  mpi_mode(2) = 0
    if( cart3d%sizes(3) == 1 )  mpi_mode(3) = 0

    ! ----- MPI-3 shared memory communication ------------------------
    if( use_shm ) then
       ! node-local mapping
       tmpA = (/ranks%north, ranks%east, ranks%south, ranks%west, ranks%up, ranks%down/)
       call mpi_comm_split_type(cart3d%comm, mpi_comm_type_shared, 0, mpi_info_null, comm_local, merr)
       call mpi_comm_group(cart3d%comm, group_world, merr)
       call mpi_comm_group(comm_local,  group_local, merr)
       call mpi_group_translate_ranks(group_world, 6, tmpA, group_local, tmpB, merr)
       call mpi_comm_size(comm_local, ranks_local%size, merr)
       call mpi_comm_rank(comm_local, ranks_local%myrank, merr)
       do i=1,6
          if( tmpB(i) == mpi_proc_null )  tmpB(i) = mpi_undefined
       enddo
       ranks_local%north = tmpB(1) ;   ranks_local%east = tmpB(2)
       ranks_local%south = tmpB(3) ;   ranks_local%west = tmpB(4)
       ranks_local%up    = tmpB(5) ;   ranks_local%down = tmpB(6)

       ! preparing shared memories for west <--> east communication
       if(( ranks_local%size > 1 ).and.( mpi_mode(1) == 1 )) then

          ! NOTE: This value should be shared in the local node, because
          ! the local mpi_barrier is critical in the communication part.
          mpi_mode(1) = 3

          if(( ranks_local%west /= mpi_undefined ).or.( ranks_local%east /= mpi_undefined )) then
             msize = 2*8*jx*kx*var1
             mdisp = 0
             call mpi_win_allocate_shared(msize,8,mpi_info_null,comm_local,baseptr1,mwin1,merr)
             call c_f_pointer(baseptr1,fptr1,(/jx,kx,var1,2/))
             if( ranks_local%west /= mpi_undefined ) then
                call mpi_win_shared_query(mwin1,ranks_local%west,msize,mdisp,baseptr1,merr)
                call c_f_pointer(baseptr1,fwest,(/jx,kx,var1,2/))
             endif
             if( ranks_local%east /= mpi_undefined ) then
                call mpi_win_shared_query(mwin1,ranks_local%east,msize,mdisp,baseptr1,merr)
                call c_f_pointer(baseptr1,feast,(/jx,kx,var1,2/))
             endif
          else
             ! dummy pointer
             msize = 0
             call mpi_win_allocate_shared(msize,8,mpi_info_null,comm_local,baseptr1,mwin1,merr)
          endif

       endif

       ! preparing shared memories for south <--> north communication
       if(( ranks_local%size > 1 ).and.( mpi_mode(2) == 1 )) then

          ! NOTE: This value should be shared in the local node, because
          ! the local mpi_barrier is critical in the communication part.
          mpi_mode(2) = 3

          if(( ranks_local%north /= mpi_undefined ).or.( ranks_local%south /= mpi_undefined )) then
             msize = 2*8*ix*kx*var1
             mdisp = 0
             call mpi_win_allocate_shared(msize,8,mpi_info_null,comm_local,baseptr2,mwin2,merr)
             call c_f_pointer(baseptr2,fptr2,(/ix,kx,var1,2/))
             if( ranks_local%south /= mpi_undefined ) then
                call mpi_win_shared_query(mwin2,ranks_local%south,msize,mdisp,baseptr2,merr)
                call c_f_pointer(baseptr2,fsouth,(/ix,kx,var1,2/))
             endif
             if( ranks_local%north /= mpi_undefined ) then
                call mpi_win_shared_query(mwin2,ranks_local%north,msize,mdisp,baseptr2,merr)
                call c_f_pointer(baseptr2,fnorth,(/ix,kx,var1,2/))
             endif
          else
             ! dummy pointer
             msize = 0
             call mpi_win_allocate_shared(msize,8,mpi_info_null,comm_local,baseptr2,mwin2,merr)
          endif

       endif
       
       ! preparing shared memories for downstairs <--> upstairs communication
       if(( ranks_local%size > 1 ).and.( mpi_mode(3) == 1 )) then

          ! NOTE: This value should be shared in the local node, because
          ! the local mpi_barrier is critical in the communication part.
          mpi_mode(3) = 3

          if(( ranks_local%down /= mpi_undefined ).or.( ranks_local%up /= mpi_undefined )) then
             msize = 2*8*ix*jx*var1
             mdisp = 0
             call mpi_win_allocate_shared(msize,8,mpi_info_null,comm_local,baseptr3,mwin3,merr)
             call c_f_pointer(baseptr3,fptr3,(/ix,jx,var1,2/))
             if( ranks_local%down /= mpi_undefined ) then
                call mpi_win_shared_query(mwin3,ranks_local%down,msize,mdisp,baseptr3,merr)
                call c_f_pointer(baseptr3,fdown,(/ix,jx,var1,2/))
             endif
             if( ranks_local%up /= mpi_undefined ) then
                call mpi_win_shared_query(mwin3,ranks_local%up,msize,mdisp,baseptr3,merr)
                call c_f_pointer(baseptr3,fup,(/ix,jx,var1,2/))
             endif
          else
             ! dummy pointer
             msize = 0
             call mpi_win_allocate_shared(msize,8,mpi_info_null,comm_local,baseptr3,mwin3,merr)
          endif

       endif
    endif
    ! ----- MPI-3 shared memory communication ------------------------

    return
  end subroutine parallel_init


  subroutine parallel_finalize
    integer :: merr

    if( mpi_mode(1) == 3 )  call mpi_win_free(mwin1,merr)
    if( mpi_mode(2) == 3 )  call mpi_win_free(mwin2,merr)
    if( mpi_mode(3) == 3 )  call mpi_win_free(mwin3,merr)
    call mpi_finalize(merr)
  
  end subroutine parallel_finalize


  subroutine parallel_exchange(U,ix,jx,kx,dir)
    include 'param.h'
    integer, intent(in) :: ix, jx, kx
    real(8), intent(inout) :: U(ix,jx,kx,var1)
    ! direction [input]: 1 (X), 2 (Y), 3 (Z)
    integer, intent(in)  :: dir
!----------------------------------------------------------------------
    real(8), allocatable :: bufsnd1(:,:,:), bufrcv1(:,:,:)
    real(8), allocatable :: bufsnd2(:,:,:), bufrcv2(:,:,:)
    integer :: mreq1(2), mreq2(2)
    integer :: msize, merr
!----------------------------------------------------------------------

    select case(dir)
    case(1)  ! west <--> east

       allocate( bufsnd1(jx,kx,var1), bufrcv1(jx,kx,var1) )
       allocate( bufsnd2(jx,kx,var1), bufrcv2(jx,kx,var1) )
       msize  = jx*kx*var1

       ! MPI mode switch
       select case(mpi_mode(1))
       case(0)  ! no MPI

          ! periodic
          if( ranks%west /= mpi_proc_null ) then
             U(1,:,:,:)  = U(ix-1,:,:,:)
             U(ix,:,:,:) = U(2,:,:,:)
          endif

       case(1)  ! MPI-1

          ! nonblocking communication (mreq1)
          call mpi_irecv(bufrcv1,msize,mpi_real8,ranks%east,0,cart3d%comm,mreq1(1),merr)
          bufsnd1(:,:,:) = U(2,:,:,:)
          call mpi_isend(bufsnd1,msize,mpi_real8,ranks%west,0,cart3d%comm,mreq1(2),merr)
          ! nonblocking communication (mreq2)
          call mpi_irecv(bufrcv2,msize,mpi_real8,ranks%west,1,cart3d%comm,mreq2(1),merr)
          bufsnd2(:,:,:) = U(ix-1,:,:,:)
          call mpi_isend(bufsnd2,msize,mpi_real8,ranks%east,1,cart3d%comm,mreq2(2),merr)

          call mpi_waitall(2,mreq1,mpi_statuses_ignore,merr)
          if( ranks%east /= mpi_proc_null ) then
             U(ix,:,:,:) = bufrcv1(:,:,:)
          endif
          call mpi_waitall(2,mreq2,mpi_statuses_ignore,merr)
          if( ranks%west /= mpi_proc_null ) then
             U(1,:,:,:) = bufrcv2(:,:,:)
          endif

       case(3)  ! MPI-3

          call mpi_win_lock_all(0,mwin1,merr)
          if( ranks_local%west /= mpi_undefined ) then
!            call mpi_win_lock(mpi_lock_shared,ranks_local%west,0,mwin1,merr)
             fwest(:,:,:,2) = U(2,:,:,:)
!            call mpi_win_unlock(ranks_local%west,mwin1,merr)
          else
             call mpi_irecv(bufrcv2,msize,mpi_real8,ranks%west,1,cart3d%comm,mreq1(1),merr)
             bufsnd1(:,:,:) = U(2,:,:,:)
             call mpi_isend(bufsnd1,msize,mpi_real8,ranks%west,0,cart3d%comm,mreq1(2),merr)
          endif
          if( ranks_local%east /= mpi_undefined ) then
!            call mpi_win_lock(mpi_lock_shared,ranks_local%east,0,mwin1,merr)
             feast(:,:,:,1) = U(ix-1,:,:,:)
!            call mpi_win_unlock(ranks_local%east,mwin1,merr)
          else
             call mpi_irecv(bufrcv1,msize,mpi_real8,ranks%east,0,cart3d%comm,mreq2(1),merr)
             bufsnd2(:,:,:) = U(ix-1,:,:,:)
             call mpi_isend(bufsnd2,msize,mpi_real8,ranks%east,1,cart3d%comm,mreq2(2),merr)
          endif

          call mpi_win_unlock_all(mwin1,merr)
          call mpi_barrier(comm_local,merr)  ! local barrier
!         call mpi_win_sync(mwin1,merr)
!         call mpi_win_fence(mpi_mode_noput,mwin1,merr)

          if( ranks_local%west /= mpi_undefined ) then
             U(1,:,:,:) = fptr1(:,:,:,1)
          else
             call mpi_waitall(2,mreq1,mpi_statuses_ignore,merr)
             if( ranks%west /= mpi_proc_null ) then
                U(1,:,:,:) = bufrcv2(:,:,:)
             endif
          endif
          if( ranks_local%east /= mpi_undefined ) then
             U(ix,:,:,:) = fptr1(:,:,:,2)
          else
             call mpi_waitall(2,mreq2,mpi_statuses_ignore,merr)
             if( ranks%east /= mpi_proc_null ) then
                U(ix,:,:,:) = bufrcv1(:,:,:)
             endif
          endif
!         call mpi_win_fence(0,mwin1,merr)

       end select
       ! MPI mode switch

       deallocate( bufsnd1, bufrcv1 )
       deallocate( bufsnd2, bufrcv2 )

    !----------------------------------------------------------
    case(2)  ! south <--> north

       allocate( bufsnd1(ix,kx,var1), bufrcv1(ix,kx,var1) )
       allocate( bufsnd2(ix,kx,var1), bufrcv2(ix,kx,var1) )
       msize  = ix*kx*var1

       ! MPI mode switch
       select case(mpi_mode(2))
       case(0)  ! no MPI

          ! periodic
          if( ranks%south /= mpi_proc_null ) then
             U(:,jx,:,:) = U(:,   2,:,:)
             U(:, 1,:,:) = U(:,jx-1,:,:)
          endif

       case(1)  ! MPI-1

          ! nonblocking communication (mreq1)
          call mpi_irecv(bufrcv1,msize,mpi_real8,ranks%north,0,cart3d%comm,mreq1(1),merr)
          bufsnd1(:,:,:) = U(:,2,:,:)
          call mpi_isend(bufsnd1,msize,mpi_real8,ranks%south,0,cart3d%comm,mreq1(2),merr)
          ! nonblocking communication (mreq2)
          call mpi_irecv(bufrcv2,msize,mpi_real8,ranks%south,1,cart3d%comm,mreq2(1),merr)
          bufsnd2(:,:,:) = U(:,jx-1,:,:)
          call mpi_isend(bufsnd2,msize,mpi_real8,ranks%north,1,cart3d%comm,mreq2(2),merr)

          call mpi_waitall(2,mreq1,mpi_statuses_ignore,merr)
          if( ranks%north /= mpi_proc_null ) then
             U(:,jx,:,:) = bufrcv1(:,:,:)
          endif
          call mpi_waitall(2,mreq2,mpi_statuses_ignore,merr)
          if( ranks%south /= mpi_proc_null ) then
             U(:,1,:,:) = bufrcv2(:,:,:)
          endif

       case(3)  ! MPI-3

          call mpi_win_lock_all(0,mwin2,merr)
          if( ranks_local%south /= mpi_undefined ) then
             fsouth(:,:,:,2) = U(:,2,:,:)
          else
             call mpi_irecv(bufrcv2,msize,mpi_real8,ranks%south,1,cart3d%comm,mreq1(1),merr)
             bufsnd1(:,:,:) = U(:,2,:,:)
             call mpi_isend(bufsnd1,msize,mpi_real8,ranks%south,0,cart3d%comm,mreq1(2),merr)
          endif
          if( ranks_local%north /= mpi_undefined ) then
             fnorth(:,:,:,1) = U(:,jx-1,:,:)
          else
             call mpi_irecv(bufrcv1,msize,mpi_real8,ranks%north,0,cart3d%comm,mreq2(1),merr)
             bufsnd2(:,:,:) = U(:,jx-1,:,:)
             call mpi_isend(bufsnd2,msize,mpi_real8,ranks%north,1,cart3d%comm,mreq2(2),merr)
          endif

          call mpi_win_unlock_all(mwin2,merr)
          call mpi_barrier(comm_local,merr)  ! local barrier

          if( ranks_local%south /= mpi_undefined ) then
             U(:,1,:,:) = fptr2(:,:,:,1)
          else
             call mpi_waitall(2,mreq1,mpi_statuses_ignore,merr)
             if( ranks%south /= mpi_proc_null ) then
                U(:,1,:,:) = bufrcv2(:,:,:)
             endif
          endif
          if( ranks_local%north /= mpi_undefined ) then
             U(:,jx,:,:) = fptr2(:,:,:,2)
          else
             call mpi_waitall(2,mreq2,mpi_statuses_ignore,merr)
             if( ranks%north /= mpi_proc_null ) then
                U(:,jx,:,:) = bufrcv1(:,:,:)
             endif
          endif

       end select
       ! MPI mode switch

       deallocate( bufsnd1, bufrcv1 )
       deallocate( bufsnd2, bufrcv2 )

    !----------------------------------------------------------
    case(3)  ! downstairs <--> upstairs

       allocate( bufsnd1(ix,jx,var1), bufrcv1(ix,jx,var1) )
       allocate( bufsnd2(ix,jx,var1), bufrcv2(ix,jx,var1) )
       msize  = ix*jx*var1

       ! MPI mode switch
       select case(mpi_mode(3))
       case(0)  ! no MPI

          ! periodic
          if( ranks%down /= mpi_proc_null ) then
             U(:,:,kx,:) = U(:,:,2,:)
             U(:,:,1,:)  = U(:,:,kx-1,:)
          endif

       case(1)  ! MPI-1

          ! nonblocking communication (mreq1)
          call mpi_irecv(bufrcv1,msize,mpi_real8,ranks%up,0,cart3d%comm,mreq1(1),merr)
          bufsnd1(:,:,:) = U(:,:,2,:)
          call mpi_isend(bufsnd1,msize,mpi_real8,ranks%down,0,cart3d%comm,mreq1(2),merr)
          ! nonblocking communication (mreq2)
          call mpi_irecv(bufrcv2,msize,mpi_real8,ranks%down,1,cart3d%comm,mreq2(1),merr)
          bufsnd2(:,:,:) = U(:,:,jx-1,:)
          call mpi_isend(bufsnd2,msize,mpi_real8,ranks%up,1,cart3d%comm,mreq2(2),merr)

          call mpi_waitall(2,mreq1,mpi_statuses_ignore,merr)
          if( ranks%up /= mpi_proc_null ) then
             U(:,:,kx,:) = bufrcv1(:,:,:)
          endif
          call mpi_waitall(2,mreq2,mpi_statuses_ignore,merr)
          if( ranks%down /= mpi_proc_null ) then
             U(:,:,1,:) = bufrcv2(:,:,:)
          endif

       case(3)  ! MPI-3

          call mpi_win_lock_all(0,mwin3,merr)
          if( ranks_local%down /= mpi_undefined ) then
             fdown(:,:,:,2) = U(:,:,2,:)
          else
             call mpi_irecv(bufrcv2,msize,mpi_real8,ranks%down,1,cart3d%comm,mreq1(1),merr)
             bufsnd1(:,:,:) = U(:,:,2,:)
             call mpi_isend(bufsnd1,msize,mpi_real8,ranks%down,0,cart3d%comm,mreq1(2),merr)
          endif
          if( ranks_local%up /= mpi_undefined ) then
             fup(:,:,:,1) = U(:,:,kx-1,:)
          else
             call mpi_irecv(bufrcv1,msize,mpi_real8,ranks%up,0,cart3d%comm,mreq2(1),merr)
             bufsnd2(:,:,:) = U(:,:,kx-1,:)
             call mpi_isend(bufsnd2,msize,mpi_real8,ranks%up,1,cart3d%comm,mreq2(2),merr)
          endif

          call mpi_win_unlock_all(mwin3,merr)
          call mpi_barrier(comm_local,merr)  ! local barrier

          if( ranks_local%down /= mpi_undefined ) then
             U(:,:,1,:) = fptr3(:,:,:,1)
          else
             call mpi_waitall(2,mreq1,mpi_statuses_ignore,merr)
             if( ranks%down /= mpi_proc_null ) then
                U(:,:,1,:) = bufrcv2(:,:,:)
             endif
          endif
          if( ranks_local%up /= mpi_undefined ) then
             U(:,:,kx,:) = fptr3(:,:,:,2)
          else
             call mpi_waitall(2,mreq2,mpi_statuses_ignore,merr)
             if( ranks%up /= mpi_proc_null ) then
                U(:,:,kx,:) = bufrcv1(:,:,:)
             endif
          endif

       end select
       ! MPI mode switch

       deallocate( bufsnd1, bufrcv1 )
       deallocate( bufsnd2, bufrcv2 )
   
    end select

  end subroutine parallel_exchange


  subroutine parallel_exchange2(VL,VR,ix,jx,kx,dir)
    include 'param.h'
    integer, intent(in) :: ix, jx, kx
    real(8), intent(inout) :: VL(ix,jx,kx,var1)  ! left values (VL)
    real(8), intent(inout) :: VR(ix,jx,kx,var1)  ! right values (VR)
    ! direction [input]: 1 (X), 2 (Y), 3 (Z)
    integer, intent(in)  :: dir
!----------------------------------------------------------------------
    real(8), allocatable :: bufsnd1(:,:,:), bufrcv1(:,:,:)
    real(8), allocatable :: bufsnd2(:,:,:), bufrcv2(:,:,:)
    integer :: mreq1(2), mreq2(2)
    integer :: msize, merr
!----------------------------------------------------------------------

    select case(dir)
    case(1)  ! west <--> east

       allocate( bufsnd1(jx,kx,var1), bufrcv1(jx,kx,var1) )
       allocate( bufsnd2(jx,kx,var1), bufrcv2(jx,kx,var1) )
       msize  = jx*kx*var1

       ! MPI mode switch
       select case(mpi_mode(1))
       case(0)  ! no MPI

          ! periodic
          if( ranks%west /= mpi_proc_null ) then
             VR(ix-1,:,:,:) = VR(1,:,:,:)
             VL(1,:,:,:)    = VL(ix-1,:,:,:)
          endif

       case(1)  ! MPI-1

          ! nonblocking communication (mreq1)
          call mpi_irecv(bufrcv1,msize,mpi_real8,ranks%east,0,cart3d%comm,mreq1(1),merr)
          bufsnd1(:,:,:) = VR(1,:,:,:)
          call mpi_isend(bufsnd1,msize,mpi_real8,ranks%west,0,cart3d%comm,mreq1(2),merr)
          ! nonblocking communication (mreq2)
          call mpi_irecv(bufrcv2,msize,mpi_real8,ranks%west,1,cart3d%comm,mreq2(1),merr)
          bufsnd2(:,:,:) = VL(ix-1,:,:,:)
          call mpi_isend(bufsnd2,msize,mpi_real8,ranks%east,1,cart3d%comm,mreq2(2),merr)

          call mpi_waitall(2,mreq1,mpi_statuses_ignore,merr)
          if( ranks%east /= mpi_proc_null ) then
             VR(ix-1,:,:,:) = bufrcv1(:,:,:)
          endif
          call mpi_waitall(2,mreq2,mpi_statuses_ignore,merr)
          if( ranks%west /= mpi_proc_null ) then
             VL(1,:,:,:) = bufrcv2(:,:,:)
          endif

       case(3)  ! MPI-3

          call mpi_win_lock_all(0,mwin1,merr)
          if( ranks_local%west /= mpi_undefined ) then
             fwest(:,:,:,2) = VR(1,:,:,:)
          else
             call mpi_irecv(bufrcv2,msize,mpi_real8,ranks%west,1,cart3d%comm,mreq1(1),merr)
             bufsnd1(:,:,:) = VR(1,:,:,:)
             call mpi_isend(bufsnd1,msize,mpi_real8,ranks%west,0,cart3d%comm,mreq1(2),merr)
          endif
          if( ranks_local%east /= mpi_undefined ) then
             feast(:,:,:,1) = VL(ix-1,:,:,:)
          else
             call mpi_irecv(bufrcv1,msize,mpi_real8,ranks%east,0,cart3d%comm,mreq2(1),merr)
             bufsnd2(:,:,:) = VL(ix-1,:,:,:)
             call mpi_isend(bufsnd2,msize,mpi_real8,ranks%east,1,cart3d%comm,mreq2(2),merr)
          endif

          call mpi_win_unlock_all(mwin1,merr)
          call mpi_barrier(comm_local,merr)  ! local barrier

          if( ranks_local%west /= mpi_undefined ) then
             VL(1,:,:,:) = fptr1(:,:,:,1)
          else
             call mpi_waitall(2,mreq1,mpi_statuses_ignore,merr)
             if( ranks%west /= mpi_proc_null ) then
                VL(1,:,:,:) = bufrcv2(:,:,:)
             endif
          endif
          if( ranks_local%east /= mpi_undefined ) then
             VR(ix-1,:,:,:) = fptr1(:,:,:,2)
          else
             call mpi_waitall(2,mreq2,mpi_statuses_ignore,merr)
             if( ranks%east /= mpi_proc_null ) then
                VR(ix-1,:,:,:) = bufrcv1(:,:,:)
             endif
          endif

       end select
       ! MPI mode switch

       deallocate( bufsnd1, bufrcv1 )
       deallocate( bufsnd2, bufrcv2 )

    !----------------------------------------------------------
    case(2)  ! south <--> north

       allocate( bufsnd1(ix,kx,var1), bufrcv1(ix,kx,var1) )
       allocate( bufsnd2(ix,kx,var1), bufrcv2(ix,kx,var1) )
       msize  = ix*kx*var1

       ! MPI mode switch
       select case(mpi_mode(2))
       case(0)  ! no MPI

          ! periodic
          if( ranks%south /= mpi_proc_null ) then
             VR(:,jx-1,:,:) = VR(:,1,:,:)
             VL(:,   1,:,:) = VL(:,jx-1,:,:)
          endif

       case(1)  ! MPI-1

          ! nonblocking communication (mreq1)
          call mpi_irecv(bufrcv1,msize,mpi_real8,ranks%north,0,cart3d%comm,mreq1(1),merr)
          bufsnd1(:,:,:) = VR(:,1,:,:)
          call mpi_isend(bufsnd1,msize,mpi_real8,ranks%south,0,cart3d%comm,mreq1(2),merr)
          ! nonblocking communication (mreq2)
          call mpi_irecv(bufrcv2,msize,mpi_real8,ranks%south,1,cart3d%comm,mreq2(1),merr)
          bufsnd2(:,:,:) = VL(:,jx-1,:,:)
          call mpi_isend(bufsnd2,msize,mpi_real8,ranks%north,1,cart3d%comm,mreq2(2),merr)

          call mpi_waitall(2,mreq1,mpi_statuses_ignore,merr)
          if( ranks%north /= mpi_proc_null ) then
             VR(:,jx-1,:,:) = bufrcv1(:,:,:)
          endif
          call mpi_waitall(2,mreq2,mpi_statuses_ignore,merr)
          if( ranks%south /= mpi_proc_null ) then
             VL(:,1,:,:) = bufrcv2(:,:,:)
          endif

       case(3)  ! MPI-3

          call mpi_win_lock_all(0,mwin2,merr)
          if( ranks_local%south /= mpi_undefined ) then
             fsouth(:,:,:,2) = VR(:,1,:,:)
          else
             call mpi_irecv(bufrcv2,msize,mpi_real8,ranks%south,1,cart3d%comm,mreq1(1),merr)
             bufsnd1(:,:,:) = VR(:,1,:,:)
             call mpi_isend(bufsnd1,msize,mpi_real8,ranks%south,0,cart3d%comm,mreq1(2),merr)
          endif
          if( ranks_local%north /= mpi_undefined ) then
             fnorth(:,:,:,1) = VL(:,jx-1,:,:)
          else
             call mpi_irecv(bufrcv1,msize,mpi_real8,ranks%north,0,cart3d%comm,mreq2(1),merr)
             bufsnd2(:,:,:) = VL(:,jx-1,:,:)
             call mpi_isend(bufsnd2,msize,mpi_real8,ranks%north,1,cart3d%comm,mreq2(2),merr)
          endif

          call mpi_win_unlock_all(mwin2,merr)
          call mpi_barrier(comm_local,merr)  ! local barrier

          if( ranks_local%south /= mpi_undefined ) then
             VL(:,1,:,:) = fptr2(:,:,:,1)
          else
             call mpi_waitall(2,mreq1,mpi_statuses_ignore,merr)
             if( ranks%south /= mpi_proc_null ) then
                VL(:,1,:,:) = bufrcv2(:,:,:)
             endif
          endif
          if( ranks_local%north /= mpi_undefined ) then
             VR(:,jx-1,:,:) = fptr2(:,:,:,2)
          else
             call mpi_waitall(2,mreq2,mpi_statuses_ignore,merr)
             if( ranks%north /= mpi_proc_null ) then
                VR(:,jx-1,:,:) = bufrcv1(:,:,:)
             endif
          endif

       end select
       ! MPI mode switch

       deallocate( bufsnd1, bufrcv1 )
       deallocate( bufsnd2, bufrcv2 )

    !----------------------------------------------------------
    case(3)  ! downstairs <--> upstairs

       allocate( bufsnd1(ix,jx,var1), bufrcv1(ix,jx,var1) )
       allocate( bufsnd2(ix,jx,var1), bufrcv2(ix,jx,var1) )
       msize  = ix*jx*var1

       ! MPI mode switch
       select case(mpi_mode(3))
       case(0)  ! no MPI

          ! periodic
          if( ranks%down /= mpi_proc_null ) then
             VR(:,:,kx-1,:) = VR(:,:,   1,:)
             VL(:,:,   1,:) = VL(:,:,kx-1,:)
          endif

       case(1)  ! MPI-1

          ! nonblocking communication (mreq1)
          call mpi_irecv(bufrcv1,msize,mpi_real8,ranks%up,0,cart3d%comm,mreq1(1),merr)
          bufsnd1(:,:,:) = VR(:,:,1,:)
          call mpi_isend(bufsnd1,msize,mpi_real8,ranks%down,0,cart3d%comm,mreq1(2),merr)
          ! nonblocking communication (mreq2)
          call mpi_irecv(bufrcv2,msize,mpi_real8,ranks%down,1,cart3d%comm,mreq2(1),merr)
          bufsnd2(:,:,:) = VL(:,:,kx-1,:)
          call mpi_isend(bufsnd2,msize,mpi_real8,ranks%up,1,cart3d%comm,mreq2(2),merr)

          call mpi_waitall(2,mreq1,mpi_statuses_ignore,merr)
          if( ranks%up /= mpi_proc_null ) then
             VR(:,:,kx-1,:) = bufrcv1(:,:,:)
          endif
          call mpi_waitall(2,mreq2,mpi_statuses_ignore,merr)
          if( ranks%down /= mpi_proc_null ) then
             VL(:,:,1,:) = bufrcv2(:,:,:)
          endif

       case(3)  ! MPI-3

          call mpi_win_lock_all(0,mwin3,merr)
          if( ranks_local%down /= mpi_undefined ) then
             fdown(:,:,:,2) = VR(:,:,1,:)
          else
             call mpi_irecv(bufrcv2,msize,mpi_real8,ranks%down,1,cart3d%comm,mreq1(1),merr)
             bufsnd1(:,:,:) = VR(:,:,1,:)
             call mpi_isend(bufsnd1,msize,mpi_real8,ranks%down,0,cart3d%comm,mreq1(2),merr)
          endif
          if( ranks_local%up /= mpi_undefined ) then
             fup(:,:,:,1) = VL(:,:,kx-1,:)
          else
             call mpi_irecv(bufrcv1,msize,mpi_real8,ranks%up,0,cart3d%comm,mreq2(1),merr)
             bufsnd2(:,:,:) = VL(:,:,kx-1,:)
             call mpi_isend(bufsnd2,msize,mpi_real8,ranks%up,1,cart3d%comm,mreq2(2),merr)
          endif

          call mpi_win_unlock_all(mwin3,merr)
          call mpi_barrier(comm_local,merr)  ! local barrier

          if( ranks_local%down /= mpi_undefined ) then
             VL(:,:,1,:) = fptr3(:,:,:,1)
          else
             call mpi_waitall(2,mreq1,mpi_statuses_ignore,merr)
             if( ranks%down /= mpi_proc_null ) then
                VL(:,:,1,:) = bufrcv2(:,:,:)
             endif
          endif
          if( ranks_local%up /= mpi_undefined ) then
             VR(:,:,kx-1,:) = fptr3(:,:,:,2)
          else
             call mpi_waitall(2,mreq2,mpi_statuses_ignore,merr)
             if( ranks%up /= mpi_proc_null ) then
                VR(:,:,kx-1,:) = bufrcv1(:,:,:)
             endif
          endif

       end select
       ! MPI mode switch

       deallocate( bufsnd1, bufrcv1 )
       deallocate( bufsnd2, bufrcv2 )
   
    end select

  end subroutine parallel_exchange2

end module parallel
