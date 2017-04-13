subroutine output(filename,ix,jx,t,x,y,U,V)
  implicit none
  include 'param.h'
  integer, intent(in) :: ix, jx
  real(8), intent(in) :: x(ix), y(jx)
  real(8), intent(in) :: U(ix,jx,var1)
  real(8), intent(in) :: V(ix,jx,var2)
  real(8), intent(in) :: t
  integer :: i, j
  character*256 :: filename      

  open(16,file=filename,form='unformatted',access='stream')
  write(16) t
  write(16) ix
  write(16) jx
  write(16) (x(i),i=1,ix)
  write(16) (y(j),j=1,jx)
  write(16) ((U(i,j,mx),i=1,ix),j=1,jx)
  write(16) ((U(i,j,my),i=1,ix),j=1,jx)
  write(16) ((U(i,j,mz),i=1,ix),j=1,jx)
  write(16) ((U(i,j,en),i=1,ix),j=1,jx)
  write(16) ((U(i,j,ro),i=1,ix),j=1,jx)
  write(16) ((U(i,j,bx),i=1,ix),j=1,jx)
  write(16) ((U(i,j,by),i=1,ix),j=1,jx)
  write(16) ((U(i,j,bz),i=1,ix),j=1,jx)
  write(16) ((U(i,j,ps),i=1,ix),j=1,jx)
  write(16) ((V(i,j,vx),i=1,ix),j=1,jx)
  write(16) ((V(i,j,vy),i=1,ix),j=1,jx)
  write(16) ((V(i,j,vz),i=1,ix),j=1,jx)
  write(16) ((V(i,j,pr),i=1,ix),j=1,jx)
  close(16)

  return
end subroutine output
