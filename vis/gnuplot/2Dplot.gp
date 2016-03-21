# tested on gnuplot 5.0.3 x11 terminal
# gnuplot> load "2Dplot.gp"

unset key
set multiplot
set dgrid3d 30,100,1

set view map
set ticslevel 0
#splot "jz.dat" binary format='%float%float%float' us 1:2:3 w pm3d
splot "ro.dat" binary format='%float%float%float' us 1:2:3 w pm3d

unset colorbox
set contour
unset surface
lcval=0  # color of contour lines
nline=20 # number of contour lines
set cntrparam levels nline
set cntrparam cubicspline
set style increment user
do for [num=1 : nline+1: 1]{ set style line num lc lcval } # set all contour colors to "lcval"
splot "Az.dat" binary format='%float%float%float' us 1:2:3 w l

unset contour
set colorbox
set surface
unset multiplot
set key

# end
