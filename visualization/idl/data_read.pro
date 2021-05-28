pro data_read,data,x,y,t,arg1,ix1=ix1,ix2=ix2,jx1=jx1,jx2=jx2,xrange=xrange,yrange=yrange

ix=0
jx=0
t0=0.d0
ix0=0L
jx0=0L

if isa(arg1, 'int') then begin
   filename = 'data/field-'+string(arg1,format='(i05)')+'.dat'
endif else if isa(arg1, 'string') then begin
   filename = arg1
endif

; no record marker
openr,/get_lun,unit,filename ;,/swap_endian
; with record marker (/f77)
;openr,/get_lun,/f77,unit,filename ;,/swap_endian

print,'reading '+filename+'...'
readu,unit,t0
readu,unit,ix0
readu,unit,jx0
print,' t = ', t0
print," size = (",ix0," x ",jx0,")"

ix=ix0
jx=jx0
t =t0

tmpx=dblarr(ix)
tmpy=dblarr(jx)
tmp =dblarr(ix,jx)
readu,unit,tmpx
readu,unit,tmpy

if not keyword_set(ix1) then ix1=0
if not keyword_set(jx1) then jx1=0
if not keyword_set(ix2) then ix2=ix-1
if not keyword_set(jx2) then jx2=jx-1
if keyword_set(xrange) then begin
   for i=ix0-1,0,-1 do begin
      if( tmpx[i] lt xrange[0] ) then begin
         ix1 = i
         break
      endif
   endfor
   for i=0,ix0-1 do begin
      if( tmpx[i] gt xrange[1] ) then begin
         ix2 = i
         break
      endif
   endfor
endif
if keyword_set(yrange) then begin
   for j=jx0-1,0,-1 do begin
      if( tmpy[j] lt yrange[0] ) then begin
         jx1 = j
         break
      endif
   endfor
   for j=0,jx0-1 do begin
      if( tmpy[j] gt yrange[1] ) then begin
         jx2 = j
         break
      endif
   endfor
endif
;; print, ix1, tmpx[ix1], ix2, tmpx[ix2]
;; print, jx1, tmpy[jx1], jx2, tmpx[jx2]

ix = ix2-ix1+1
jx = jx2-jx1+1
data = dblarr(ix,jx,9)

x = dblarr(ix) & x = tmpx[ix1:ix2]
y = dblarr(jx) & y = tmpy[jx1:jx2]
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
