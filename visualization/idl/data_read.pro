pro data_read,data,x,y,t,it,ix1=ix1,ix2=ix2,jx1=jx1,jx2=jx2

ix=0
jx=0
if not keyword_set(ix1) then ix1=0
if not keyword_set(jx1) then jx1=0
t0=0.d0
ix0=0
jx0=0

spawn, 'printf "data/field-%05d.dat" ' + string(it), retvars
filename = retvars(0)

openr,/get_lun,/f77,unit,filename ;,/swap_endian
print,'reading '+filename+'...'
readu,unit,t0
readu,unit,ix0
readu,unit,jx0
print,' t = ', t0
print," size = (",ix0," x ",jx0,")"

ix=ix0
jx=jx0
t =t0
x2=dblarr(ix)
y2=dblarr(jx)
tmp=dblarr(ix,jx)

if not keyword_set(ix2) then ix2=ix-1
if not keyword_set(jx2) then jx2=jx-1
ix=ix2-ix1+1
jx=jx2-jx1+1

x=dblarr(ix)
y=dblarr(jx)
data=dblarr(ix,jx,9)

readu,unit, x2 & x=x2[ix1:ix2]
readu,unit, y2 & y=y2[jx1:jx2]
readu,unit,tmp & readu,unit,tmp & readu,unit,tmp & readu,unit,tmp
readu,unit,tmp & data[*,*,4]=tmp[ix1:ix2,jx1:jx2]
readu,unit,tmp & data[*,*,5]=tmp[ix1:ix2,jx1:jx2]
readu,unit,tmp & data[*,*,6]=tmp[ix1:ix2,jx1:jx2]
readu,unit,tmp & data[*,*,7]=tmp[ix1:ix2,jx1:jx2]
readu,unit,tmp & data[*,*,8]=tmp[ix1:ix2,jx1:jx2]
readu,unit,tmp & data[*,*,0]=tmp[ix1:ix2,jx1:jx2]
readu,unit,tmp & data[*,*,1]=tmp[ix1:ix2,jx1:jx2]
readu,unit,tmp & data[*,*,2]=tmp[ix1:ix2,jx1:jx2]
readu,unit,tmp & data[*,*,3]=tmp[ix1:ix2,jx1:jx2]
close,unit
free_lun,unit

return
end
