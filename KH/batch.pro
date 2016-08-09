;; dummy index
vx=0 & vy=1 & vz=2 & pr=3 & ro=4 & bx=5 & by=6 & bz=7 & ps=8

;; find and compile data_read routine
resolve_routine, "data_read"
;; reading data ...
data_read,data,x,y,t,10
;data_read,data,x,y,t,8,ix1=0,ix2=1501,jx1=0,jx2=201

;; 2D image
myimg = image(data[*,*,ro],x,y,axis_style=2,xtitle='$X$',ytitle='$Y$',xtickdir=1,xticklen=0.02,ytickdir=1,yticklen=0.01,font_size=16,rgb_table=13,dimensions=[600,600])
;; ,renderer=1) ;; use software rendering over a remote connection.


;; options
;myimg.font_name = 'Times'
;myimg.font_name = 'Helivetica'
;myimg.rgb_table = 4 ;; color table
myimg.title = 'Density (t=' + string(format='(f6.1)',t) + ')'
;myimg.title.font_name = 'Times'
myimg.title.font_size = 16
;myimg.max_value = 6
;myimg.scale,1.0,2.0 ;; ARRANGE ASPECT RATIO
;myimg.yrange = [0,20]
myimg.position = [0.2,0.15,0.8,0.85]

;; preparing Vector potential
ix = (size(x))[1]
jx = (size(y))[1]
az = dblarr(ix,jx)
az[0,0] = 0.0d0
for j=1,jx-1 do begin
   az[0,j] = az[0,j-1] + 0.5*(data[0,j-1,bx]+data[0,j,bx])
endfor
for i=1,ix-1 do begin
   az[i,0] = az[i-1,0] - 0.5*(data[i-1,0,by]+data[i,0,by])
   for j=1,jx-1 do begin
      az[i,j] = az[i,j-1] + 0.5*(data[i,j-1,bx]+data[i,j,bx])
   endfor
endfor

;; contour plot
myct = contour(az,x,y,/overplot,color='white',c_label_show=0)
;myct.color='gray'

;; colorbar
mycb = colorbar(target=myimg,orientation=1,tickdir=1,minor=4,border=1,font_size='14',position=[0.925,0.15,0.95,0.85])
mycb.font_size = '14'

;; image file
;myimg.save,'output.png',resolution=200

end
