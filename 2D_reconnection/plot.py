import matplotlib.pyplot as plt
import numpy as np
import openmhd
# dummy index
vx=0;vy=1;vz=2;pr=3;ro=4;bx=5;by=6;bz=7;ps=8

# reading the data ...
x,y,t,data = openmhd.data_read(8)
#x,y,t,data = openmhd.data_read(10,ix1=0,ix2=1301,jx1=0,jx2=151)

# preparing the canvas
fig = plt.figure(figsize=(12, 6), dpi=80)
# fig.clear()
plt.clf()

# 2D image
# extent: [left, right, bottom, top]
extent=[x[0],x[-1],y[0],y[-1]]
plt.imshow(data[:,:,vx].T,origin='lower',extent=extent,aspect='auto')

# useful options
# plt.grid()
plt.xlabel("X",size=16)
plt.ylabel("Y",size=16)
plt.title('Outflow speed (t = %6.1f)' % t, size=20)
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
plt.contour(az.T,extent=extent,colors='w',linestyles='solid')

# plot
plt.show()

# image file
#plt.savefig('output.png', dpi=80)

# end
