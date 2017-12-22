import matplotlib.pyplot as plt
#import numpy as np
import openmhd
# dummy index
vx=0;vy=1;vz=2;pr=3;ro=4;bx=5;by=6;bz=7;ps=8

# reading the data ...
x,y,t,data = openmhd.data_read(15)

plt.clf()

# 2D image
# extent: [left, right, bottom, top]
# Note: ().T is necessary, because the imshow routine uses the image coordinates
extent=[x[0],x[-1],y[0],y[-1]]
myimg = plt.imshow((data[:,:,ro]).T,origin='lower',cmap='seismic',extent=extent)

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
#plt.savefig('output.png', dpi=144)

# end
