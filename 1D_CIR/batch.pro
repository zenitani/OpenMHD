;; dummy index
x=0 & y=1 & ro=2 & pr=3 & vx=4 & vy=5 & vz=6 & bx=7 & by=8 & bz=9

;; reading data ... (skipping the header)
d1 = (read_ascii('data/x-00000.dat',data_start=1)).FIELD01
d2 = (read_ascii('data/x-00006.dat',data_start=1)).FIELD01

;; Y axis (density)
;; mypl1 = plot(d1[x,*],d1[ro,*],xtitle='x',ytitle='$\rho$',name='t=000')
;; ;; ,renderer=1) ;; use software rendering over a remote connection.
;; mypl2 = plot(d2[x,*],d2[ro,*],/overplot,name='t=600')
;; mylg = legend(target=[mypl1,mypl2],position=[2000,6],/data,/auto_text_color)

;; Y axis (|B|)
mypl1 = plot(d1[x,*],sqrt(d1[bx,*]^2+d1[by,*]^2+d1[bz,*]^2),xtitle='x',ytitle='$|B|$',name='t=000')
;; ,renderer=1) ;; use software rendering over a remote connection.
mypl2 = plot(d2[x,*],sqrt(d2[bx,*]^2+d2[by,*]^2+d2[bz,*]^2),/overplot,name='t=600')
mylg = legend(target=[mypl1,mypl2],position=[2000,4],/data,/auto_text_color)

;; mypl.symbol = 'circle'
;; mypl.color  = 'red'
mypl1.font_size = 16
mypl1.xrange = [0,10000]
mypl1.linestyle = 'dash'
mypl2.linestyle = 'solid'
;; mypl1.font_name = 'Times'
;; mypl1.font_name = 'Helivetica'
;mylg.font_name = 'Times'
;; mylg.shadow = 0
mylg.font_size = '14'
mylg.horizontal_spacing = 0.05

;; image file
;mypl1.save,'output.png'
;mypl1.save,'output.png',resolution=150

end
