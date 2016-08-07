;; dummy index
x=0 & y=1 & ro=2 & pr=3 & vx=4 & vy=5 & vz=6 & bx=7 & by=8 & bz=9

;; reading data ... (skipping the header)
data = (read_ascii('data/x-00002.dat',data_start=1)).FIELD01

mypl= plot(data[x,*],data[ro,*],position=[-0.5,0.5],xtitle='x',ytitle='$\rho$')
;mypl = plot(data[x,*],data[by,*],position=[-0.5,0.5],xtitle='x',ytitle='$B_y$')

;mypl.xrange = [-0.5,0.5]
mypl.symbol = 'circle'
mypl.color  = 'red'
mypl.font_size = 16
;mypl.font_name = 'Times'
;mypl.font_name = 'Helivetica'
mylg = legend(target=mypl,position=[0.45,0.7],/data,/auto_text_color)
;mylg.font_name = 'Times'
mylg.font_size = '14'
mylg.label = '$\rho$'
;mylg.label = '$B_y$'
;mylg.shadow = 0
mylg.horizontal_spacing = 0.05

;; image file
;mypl1.save,'output.png'
;mypl1.save,'output.png',resolution=150

end
