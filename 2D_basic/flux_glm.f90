subroutine flux_glm(F,VL,VR,ch,ix,jx,dir)
!-----------------------------------------------------------------------
!     GLM flux terms in the X/Y/Z directions
!-----------------------------------------------------------------------
!     This needs to be called right after LLF/HLL/HLLC/HLLD flux solver,
!     because these solvers may set F(:,:,bx) and F(:,:,ps) to zero
!-----------------------------------------------------------------------
  implicit none
  include 'param.h'
  integer, intent(in) :: ix, jx
! numerical flux (F) [input/output]
  real(8), intent(inout) :: F(ix,jx,var1)
! left/right states (VL,VR) [input]
  real(8), intent(in) :: VL(ix,jx,var1), VR(ix,jx,var1)
! direction [input]: 1 (X), 2 (Y), 3 (Z)
  integer, intent(in)  :: dir
! div cleaning wave speed
  real(8), intent(in) :: ch
  real(8) :: ch2
  integer :: i, j, is=0, ie=0, js=0, je=0
  integer :: bn

  bn = bx + mod(dir-1,3)
  select case(dir)
  case(1)
     is = 1; ie = ix-1
     js = 2; je = jx-1
  case(2)
     is = 2; ie = ix-1
     js = 1; je = jx-1
  case(3)
  endselect

  ch2 = ch**2

!$omp parallel do private(i,j)
  do j=js,je
  do i=is,ie

     F(i,j,bn) = 0.5d0*( VL(i,j,ps) + VR(i,j,ps) &
          - ch*( VR(i,j,bn)-VL(i,j,bn) ) )
     F(i,j,ps) = 0.5d0*( ch2*( VL(i,j,bn)+VR(i,j,bn) ) &
          - ch*( VR(i,j,ps)-VL(i,j,ps) ) )

  enddo
  enddo
!$omp end parallel do

  return
end subroutine flux_glm
