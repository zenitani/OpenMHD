import matplotlib.pyplot as plt
import numpy as np
import openmhd
from pylab import *
# dummy index
vx=0;vy=1;vz=2;pr=3;ro=4;bx=5;by=6;bz=7;ps=8

x,y,t,data = openmhd.data_read(20)
#x,y,t,data = openmhd.data_read(20,ix1=0,ix2=100,jx1=11)

plt.clf()
# extent: [left, right, bottom, top]
extent=[x[0],x[-1],y[0],y[-1]]
plt.imshow(data[:,:,pr].T,origin='lower',extent=extent)

# plt.grid()
plt.xlabel("X",size=16)
plt.ylabel("Y",size=16)
plt.title('Pressure (t = %6.1f)' % t, size=20)
plt.colorbar()

# preparing Vector potential
az = np.ndarray((x.size,y.size),np.double)
az[0,0] = 0.0
for j in range(1,y.size):
    az[0,j] = az[0,j-1] + 0.5*(data[0,j-1,bx]+data[0,j,bx])

for i in range(1,x.size):
    az[i,0] = az[i-1,0] - 0.5*(data[i-1,0,by]+data[i,0,by])
    for j in range(1,y.size):
        az[i,j] = az[i,j-1] + 0.5*(data[i,j-1,bx]+data[i,j,bx])

contour(az[:,:].T,origin='lower',extent=extent)

plt.show()
#plt.savefig('output.png', dpi=144)
#plt.clf()

# end
