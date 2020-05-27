import matplotlib.pyplot as plt
import numpy as np
import openmhd
import os # for movies

# dummy index
vx=0;vy=1;vz=2;pr=3;ro=4;bx=5;by=6;bz=7;ps=8

# ---- for movies --------------------
# NOTE: "import os" is necessary
if not os.path.exists('movie'):
    print('This routine outputs image files in "./movie/".')
    os.mkdir('./movie/')

# ---- loop for movies ---------------
for ii in range(0,41):

    # reading the data ...
    x,y,t,data = openmhd.data_read("data/field-%05d.dat" % ii)
    
    # clearing the current figure, if any
    plt.clf()
    # extent: [left, right, bottom, top]
    extent=[x[0],x[-1],y[0],y[-1]]
    # 2D plot (vmin/mymin: minimum value, vmax/mymax: max value)
    # Note: ().T is necessary for 2-D plot routines (imshow/pcolormesh...)
    tmp = np.ndarray((x.size,y.size),np.double)
    tmp[:,:] = data[:,:,pr]
    mymax = max(tmp.max(), -tmp.min()) if( tmp.max() > 0.0 ) else 0.0
    mymin = min(tmp.min(), -tmp.max()) if( tmp.min() < 0.0 ) else 0.0
    myimg = plt.imshow(tmp.T,origin='lower',vmin=mymin,vmax=mymax,cmap='jet',extent=extent)

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
    plt.contour(az.T,extent=extent,colors='w',linestyles='solid')

    # plot
    # plt.show()

    # image file
    plt.savefig('movie/output-%05d.png' % ii, dpi=144)

# ---- for movies --------------------
plt.show()
print()
print('The image files should be found in "./movie/".')
