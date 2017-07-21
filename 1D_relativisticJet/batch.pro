;; dummy index
x=0 & y=1 & ro=2 & pr=3 & vx=4 & vy=5 & vz=6 & ux=7 & uy=8 & uz=9 & u0=10
bx=11 & by=12 & bz=13 & ex=14 & ey=15 & ez=16 & pt=17

;; reading the data ... (skipping the header)
d1 = (read_ascii('data/x-00000.dat',data_start=1)).FIELD01
d2 = (read_ascii('data/x-00002.dat',data_start=1)).FIELD01

;; Y axis (density)
;; mypl1 = plot(d1[x,*],d1[ro,*],xtitle='x',ytitle='$\rho$',name='t=000')
;; ;; ,renderer=1) ;; use software rendering over a remote connection.
;; mypl2 = plot(d2[x,*],d2[ro,*],/overplot,name='t=600')
;; mylg = legend(target=[mypl1,mypl2],position=[2000,6],/data,/auto_text_color)

;; Y axis (|B|)
mypl1 = plot(d1[x,*],d1[u0,*],xrange=[-0.1,0.1],xtitle='x',ytitle='$\gamma$',name='t=0.0')
;; ,renderer=1) ;; use software rendering over a remote connection.
mypl2 = plot(d2[x,*],d2[u0,*],/overplot,name='t=0.2')
;; legend
mylg = legend(target=[mypl1,mypl2],position=[-0.05,18.0],/data,/auto_text_color)

;; mypl.symbol = 'circle'
;; mypl.color  = 'red'
mypl1.font_size = 16
mypl1.linestyle = 'dash'
mypl2.linestyle = 'solid'
;; mypl1.font_name = 'Times'
;; mypl1.font_name = 'Helivetica'

;mylg.font_name = 'Times'
;; mylg.shadow = 0
mylg.font_size = '14'
;mylg.horizontal_spacing = 0.05

;; image file
;mypl1.save,'output.png'
;mypl1.save,'output.png',resolution=150

end
