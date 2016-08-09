;; dummy index
x=0 & y=1 & ro=2 & pr=3 & vx=4 & vy=5 & vz=6 & bx=7 & by=8 & bz=9

;; reading data ... (skipping the header)
data1 = (read_ascii('data/x-00000.dat',data_start=1)).FIELD01
data2 = (read_ascii('data/x-00006.dat',data_start=1)).FIELD01

mypl1 = plot(data1[x,*],data1[ro,*],xtitle='x',ytitle='$\rho$',name='t=000')
;; ,renderer=1) ;; use software rendering over a remote connection.
mypl2 = plot(data2[x,*],data2[ro,*],/overplot,name='t=600')
;; mypl.symbol = 'circle'
;; mypl.color  = 'red'
mypl1.font_size = 16
mypl1.xrange = [0,10000]
mypl1.linestyle = 'dash'
mypl2.linestyle = 'solid'
;; mypl1.font_name = 'Times'
;; mypl1.font_name = 'Helivetica'
mylg = legend(target=[mypl1,mypl2],position=[3000,6],/data,/auto_text_color)
;mylg.font_name = 'Times'
mylg.font_size = '14'
;; mylg.shadow = 0
mylg.horizontal_spacing = 0.05

;; image file
;mypl1.save,'output.png'
;mypl1.save,'output.png',resolution=150

end
