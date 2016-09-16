subroutine output(filename,ix,jx,t,x,y,U,V)
  implicit none
  include 'param.h'
  real(8) :: x(ix), y(jx)
  real(8) :: U(ix,jx,var1)
  real(8) :: V(ix,jx,var2)
  real(8) :: t
  integer :: ix, jx
  integer :: i, j
  character*256 :: filename      

!  open(16,file=filename,form='unformatted')
  open(16,file=filename)
  write(16,*) '# x(1),y(2),ro(3),pr(4),v(5-7),B(8-10), t=', t
  do j=1,jx
  do i=1,ix
     write(16,999) x(i), y(j), U(i,j,ro), V(i,j,pr), V(i,j,vx), V(i,j,vy), V(i,j,vz), &
          U(i,j,bx),U(i,j,by),U(i,j,bz)
999  format(10(1x,1p,e12.4e3))
  enddo
  enddo
  close(16)

  return
end subroutine output
