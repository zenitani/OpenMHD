subroutine output(filename,ix,jx,t,x,y,U,V)
  implicit none
  include 'param_rela.h'
  real(8) :: x(ix), y(jx)
  real(8) :: U(ix,jx,var1)
  real(8) :: V(ix,jx,var2)
  real(8) :: B2,u2,u02,u0,vB,ex,ey,ez,pt
  real(8) :: t
  integer :: ix, jx, i, j
  character*256 :: filename      

!  open(16,file=filename,form='unformatted')
  open(16,file=filename)
  write(16,*) '# x(1),y(2),ro(3),pr(4),v(5-7),u(8-10),u0(11),',&
       'B(12-14),E(15-17),pr_t(18)'
  do j=1,jx
  do i=1,ix
     B2 = dot_product( U(i,j,bx:bz), U(i,j,bx:bz) )
     u2 = dot_product( V(i,j,ux:uz), V(i,j,ux:uz) )
     u02= 1.d0 + u2
     u0 = sqrt( u02 )
     vB = dot_product( U(i,j,bx:bz), V(i,j,ux:uz) ) / u0
     ex = ( U(i,j,by)*V(i,j,uz)-U(i,j,bz)*V(i,j,uy) )/u0
     ey = ( U(i,j,bz)*V(i,j,ux)-U(i,j,bx)*V(i,j,uz) )/u0
     ez = ( U(i,j,bx)*V(i,j,uy)-U(i,j,by)*V(i,j,ux) )/u0
     pt = V(i,j,pr) + 0.5d0 * ( b2 / u02 + vB*vB )
     write(16,999) x(i), y(j), V(i,j,ro), V(i,j,pr), &
          V(i,j,ux:uz)/u0, V(i,j,ux:uz), u0, U(i,j,bx:bz),ex,ey,ez,pt
999  format(18(1x,1p,e12.4e3))
  enddo
  enddo
  close(16)

  return
end subroutine output
