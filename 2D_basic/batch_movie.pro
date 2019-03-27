;; dummy index
vx=0 & vy=1 & vz=2 & pr=3 & ro=4 & bx=5 & by=6 & bz=7 & ps=8
;; find and compile data_read routine
resolve_routine, "data_read"

;; ---- for movies --------------------
print, 'This routine outputs image files in "./movie/".'
file_mkdir,'movie' ;; if 'movie' exists, just proceed

;; ---- loop for movies ---------------
for ii=0,40,1 do begin

   ;; reading the data ...
   data_read,data,x,y,t,ii

   ;; 2D image
   myimg = image(data[*,*,pr],x,y,axis_style=2,xtitle='$X$',ytitle='$Y$', $
                 xtickdir=1,xticklen=0.02,ytickdir=1,yticklen=0.01, $
                 min_value=0,font_size=16,rgb_table=13,dimensions=[600,600])
   ;; ,renderer=1) ;; use software rendering over a remote connection.

   myimg.title = 'Pressure (t=' + string(format='(f6.1)',t) + ')'
   myimg.title.font_size = 16
   myimg.position = [0.15,0.15,0.85,0.85]

   ;; preparing Vector potential (Az)
   ix = (size(x))[1]
   jx = (size(y))[1]
   az = dblarr(ix,jx)
   ;az[0,0] = 0.0d0
   az[0,-1] = 0.5d0*(data[0,-1,bx] - data[0,-1,by]) ; Set the top-left corner (az[0.5,-1.5]) to zero
   for j=jx-1,1,-1 do az[0,j-1] = az[0,j]   - 0.5d0*(data[0,j-1,bx]+data[0,j,bx])
   for i=1,ix-1    do az[i,*]   = az[i-1,*] - 0.5d0*(data[i-1,*,by]+data[i,*,by])

   ;; contour of Az = magnetic field lines
   myct = contour(az,x,y,/overplot,color='white',c_label_show=0)

   ;; colorbar
   mycb = colorbar(target=myimg,orientation=1,tickdir=1,minor=4,border=1,font_size='14',position=[0.925,0.15,0.95,0.85])
   mycb.font_size = '14'

   ;; image file
   filename = "movie/output-" + string(ii,format='(i05)') + ".png"
   myimg.save,filename,resolution=96
   myimg.close

endfor
;; ---- loop for movies ---------------

;; ---- for movies --------------------
print, 'The image files will be found in "./movie/".'

end
