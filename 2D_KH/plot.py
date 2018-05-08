import matplotlib.pyplot as plt
#import numpy as np
import openmhd
# dummy index
vx=0;vy=1;vz=2;pr=3;ro=4;bx=5;by=6;bz=7;ps=8

# reading the data ...
x,y,t,data = openmhd.data_read(15)

# clearing the current figure, if any
plt.clf()
# extent: [left, right, bottom, top]
extent=[x[0],x[-1],y[0],y[-1]]
# Note: ().T is necessary, because the imshow routine uses the image coordinates
# 2D plot (vmin: minimum value, vmax: max value)
myimg = plt.imshow(data[:,:,ro].T,vmin=0,origin='lower',cmap='bwr',extent=extent)

# image operations (e.g. colormaps)
# myimg.set_cmap('jet')
# myimg.set_cmap('RdBu_r')  # colortable(70,/reverse) in IDL
# myimg.set_cmap('seismic')
# myimg.set_cmap('bwr')
# myimg.set_cmap('gnuplot2')

# useful options
# plt.grid()
plt.xlabel("X",size=16)
plt.ylabel("Y",size=16)
plt.title('Density (t = %6.1f)' % t, size=20)

# color bar (next to myimg)
plt.colorbar()

# flow vectors
# Note: ().T is necessary (the image coordinates)
myxsub = 40; myysub = 20; myxsub0 = int(myxsub/2); myysub0 = int(myysub/2)
myvec = plt.quiver( x[myxsub0::myxsub],
                    y[myysub0::myysub],
                    data[myxsub0::myxsub,myysub0::myysub,vx].T,
                    data[myxsub0::myxsub,myysub0::myysub,vy].T,
                    alpha=0.7,angles='xy')

# plot
plt.show()

# image file
# plt.savefig('output.png', dpi=144)

# end
