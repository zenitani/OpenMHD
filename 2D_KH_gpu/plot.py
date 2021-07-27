import matplotlib.pyplot as plt
import numpy as np
import openmhd
# dummy index
vx=0;vy=1;vz=2;pr=3;ro=4;bx=5;by=6;bz=7;ps=8

# reading the data ...
x,y,t,data = openmhd.data_read("data/field-00015.dat")
# reading the data (partial domain: [ix1,ix2] x [jx1,jx2])
# x,y,t,data = openmhd.data_read("data/field-00015.dat",ix2=101,jx1=50,jx2=101)

# clearing the current figure, if any
plt.clf()
# extent: [left, right, bottom, top]
extent=[x[0],x[-1],y[0],y[-1]]
# 2D plot (vmin/mymin: minimum value, vmax/mymax: max value)
# Note: ().T is necessary for 2-D plot routines (imshow/pcolormesh...)
tmp = np.ndarray((x.size,y.size),np.double)
tmp[:,:] = data[:,:,ro]
mymax = max(tmp.max(), -tmp.min()) if( tmp.max() > 0.0 ) else 0.0
mymin = min(tmp.min(), -tmp.max()) if( tmp.min() < 0.0 ) else 0.0
myimg = plt.imshow(tmp.T,origin='lower',vmin=mymin,vmax=mymax,cmap='bwr',extent=extent)

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
# Note: ().T is necessary for a 2-D array
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
