import matplotlib.pyplot as plt
import numpy as np
import openmhd
# dummy index
vx=0;vy=1;vz=2;pr=3;ro=4;bx=5;by=6;bz=7;ps=8

# reading the data ...
#x,y,t,data = openmhd.data_read(8)
# reading the data (partial domain: [ix1,ix2] x [jx1,jx2])
x,y,t,data = openmhd.data_read(10,ix1=0,ix2=1301,jx1=0,jx2=151)

# 2D mirroring (This depends on the BC)
ix = x.size
jx = 2*y.size-2
jxh= y.size-1
tmp  = data
data = np.ndarray((ix,jx,9),np.double)
data[:,jxh:,:]   =  tmp[:,1:,:]
data[:,0:jxh, :] =  tmp[:,-1:-jxh-1:-1, :]
data[:,0:jxh,vy] = -tmp[:,-1:-jxh-1:-1,vy]
data[:,0:jxh,vz] = -tmp[:,-1:-jxh-1:-1,vz]
data[:,0:jxh,bx] = -tmp[:,-1:-jxh-1:-1,bx]
data[:,0:jxh,ps] = -tmp[:,-1:-jxh-1:-1,ps]
tmp = y
y = np.ndarray((jx),np.double)
y[jxh:]  =  tmp[1:]
y[0:jxh] = -tmp[-1:-jxh-1:-1]

# preparing the canvas
fig = plt.figure(figsize=(12, 6), dpi=80)
# fig.clear()
plt.clf()

# 2D image
# extent: [left, right, bottom, top]
extent=[x[0],x[-1],y[0],y[-1]]
myimg = plt.imshow(data[:,:,vx].T,origin='lower',cmap='seismic',extent=extent,aspect='auto')

# image operations (e.g. color map)
#myimg.set_cmap('jet')
#myimg.set_cmap('seismic')

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
