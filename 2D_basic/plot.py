import matplotlib.pyplot as plt
import numpy as np
import openmhd
# dummy index
vx=0;vy=1;vz=2;pr=3;ro=4;bx=5;by=6;bz=7;ps=8

# reading the data ...
x,y,t,data = openmhd.data_read(20)
#x,y,t,data = openmhd.data_read(20,ix1=0,ix2=100,jx1=11)

# 2D image
plt.clf()
# extent: [left, right, bottom, top]
extent=[x[0],x[-1],y[0],y[-1]]
plt.imshow(data[:,:,pr].T,origin='lower',extent=extent)

# useful options
# plt.grid()
plt.xlabel("X",size=16)
plt.ylabel("Y",size=16)
plt.title('Pressure (t = %6.1f)' % t, size=20)
# colorbar
plt.colorbar()

# preparing Vector potential (Az)
az = np.ndarray((x.size,y.size),np.double)
# az[0,0] = 0.0
az[0,-1] = 0.5*(data[0,-1,bx] - data[0,-1,by])
for j in range(y.size-1,0,-1):
    az[0,j-1] = az[0,j] - 0.5*(data[0,j-1,bx]+data[0,j,bx])
for i in range(1,x.size):
    az[i,:] = az[i-1,:] - 0.5*(data[i-1,:,by]+data[i,:,by])

# contour of Az = magnetic field lines
plt.contour(az[:,:].T,origin='lower',extent=extent,colors='w',linestyles='solid')

# plot
plt.show()

# image file
#plt.savefig('output.png', dpi=144)

# end
