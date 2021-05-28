import matplotlib.pyplot as plt
import numpy as np
import openmhd
import gc #  garbage collection

# dummy index
vx=0;vy=1;vz=2;pr=3;ro=4;bx=5;by=6;bz=7;ps=8

# reading the data ...
# x,y,t,data = openmhd.data_read("data/field-00015.dat")
# reading the data (subdomain: [ix1,ix2] x [jx1,jx2] or xrange (x1,x2) x yrange (y1,y2))
# x,y,t,data = openmhd.data_read("data/field-00015.dat",ix1=0,ix2=301,jx1=0,jx2=51)
x,y,t,data = openmhd.data_read("data/field-00015.dat",xrange=(0.0,30.0),yrange=(0.0,5.0))

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
# releasing the memory, because this tmp could be large
del tmp
gc.collect()

tmp = y
y = np.ndarray((jx),np.double)
y[jxh:]  =  tmp[1:]
y[0:jxh] = -tmp[-1:-jxh-1:-1]

# preparing the canvas
fig = plt.figure(figsize=(10, 5), dpi=80)
# fig.clear()
plt.clf()

# extent: [left, right, bottom, top]
extent=[x[0],x[-1],y[0],y[-1]]
# 2D plot (vmin/mymin: minimum value, vmax/mymax: max value)
# Note: ().T is necessary for 2-D plot routines (imshow/pcolormesh...)
tmp = np.ndarray((x.size,y.size),np.double)
tmp[:,:] = data[:,:,vx]
mymax = max(tmp.max(), -tmp.min()) if( tmp.max() > 0.0 ) else 0.0
mymin = min(tmp.min(), -tmp.max()) if( tmp.min() < 0.0 ) else 0.0
myimg = plt.imshow(tmp.T,origin='lower',vmin=mymin,vmax=mymax,cmap='jet',extent=extent,aspect='auto')

# image operations (e.g. colormaps)
# myimg.set_cmap('jet')
# myimg.set_cmap('RdBu_r')  # colortable(70,/reverse) in IDL
# myimg.set_cmap('seismic')
# myimg.set_cmap('bwr')
# myimg.set_cmap('gist_ncar_r')
# myimg.set_cmap('Pastel1')

# useful options
# plt.grid()
plt.xlabel("X",size=16)
plt.ylabel("Y",size=16)
plt.title('Outflow speed (t = %6.1f)' % t, size=20)
# colorbar
plt.colorbar()

# preparing Vector potential (Az)
az = np.ndarray((x.size,y.size),np.double)
fx = 0.5*(x[1]-x[0])
fy = 0.5*(y[1]-y[0])
az[0,0] = (fy*data[0,0,bx] - fx*data[0,0,by])
for j in range(1,y.size):
    az[0,j] = az[0,j-1] + fy*(data[0,j-1,bx]+data[0,j,bx])
for i in range(1,x.size):
    az[i,:] = az[i-1,:] - fx*(data[i-1,:,by]+data[i,:,by])

# contour of Az = magnetic field lines
plt.contour(az.T,extent=extent,colors='w',linestyles='solid')

# plot
plt.show()

# adjusting the margins...
plt.tight_layout()

# image file
# plt.savefig('output.png', dpi=80)

# end
