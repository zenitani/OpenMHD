import matplotlib.pyplot as plt
import numpy as np
import openmhd3d as openmhd
# dummy index
vx=0;vy=1;vz=2;pr=3;ro=4;bx=5;by=6;bz=7;ps=8

# reading the data ...
x,y,z,t,data = openmhd.data_read("data/field-00010.dat", xrange=[0,130],yrange=[0,15])
# reading the data (partial domain: [ix1,ix2] x [jx1,jx2] x [kx1,kx2])
# x,y,z,t,data = openmhd.data_read("data/field-00010.dat",ix1=0,ix2=100,jx1=11)

def plot_volume(ax, array, x, y, z, xcut=0.0, ycut=0.0, zcut=0.0, cmap='jet'):
    """For a given 3d *array* plot three cuts, with help from their subplanes.
Generalized from https://matplotlib.org/stable/gallery/mplot3d/intersecting_planes.html."""

    ix0, jx0, kx0 = array.shape
    for i in range(2,ix0-1):
        if( x[i] > xcut ):
            ix1 = i; xrlist = [ slice(None,ix1), slice(ix1,None) ]
            break
    for j in range(2,jx0-1):
        if( y[j] > ycut ):
            jx1 = j; yrlist = [ slice(None,jx1), slice(jx1,None) ]
            break
    for k in range(2,kx0-1):
        if( z[k] > zcut ):
            kx1 = k; zrlist = [ slice(None,kx1), slice(kx1,None) ]
            break

    cmap = plt.get_cmap(cmap)
    mymax = max(array.max(), -array.min()) if( array.max() > 0.0 ) else 0.0
    mymin = min(array.min(), -array.max()) if( array.min() < 0.0 ) else 0.0
    facecolors = cmap((array[ix1,:,:] - mymin) / (mymax - mymin))
    Z, Y = np.meshgrid( z, y )
    X = x[ix1] * np.ones_like(Y)
    for yr in yrlist:
        for zr in zrlist:
            ax.plot_surface(X[yr,zr], Y[yr,zr], Z[yr,zr], rstride=3, cstride=3, facecolors=facecolors[yr,zr], shade=False, alpha=0.5)

    facecolors = cmap((array[:,jx1,:] - mymin) / (mymax - mymin))
    Z, X = np.meshgrid( z, x )
    Y = y[jx1] * np.ones_like(Z)
    for xr in xrlist:
        for zr in zrlist:
            ax.plot_surface(X[xr,zr], Y[xr,zr], Z[xr,zr], rstride=3, cstride=3, facecolors=facecolors[xr,zr], shade=False)

    facecolors = cmap((array[:,:,kx1] - mymin) / (mymax - mymin))
    Y, X = np.meshgrid( y, x )
    Z = z[kx1] * np.ones_like(X)
    for xr in xrlist:
        for yr in yrlist:
            ax.plot_surface(X[xr,yr], Y[xr,yr], Z[xr,yr], rstride=3, cstride=3, facecolors=facecolors[xr,yr], shade=False, alpha=0.5)



fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.view_init(elev=120, azim=-120, roll=0)

# tmp = 0.5*( data[:,:,:,bx]**2 + data[:,:,:,by]**2 + data[:,:,:,bz]**2 )
tmp = data[::3,::2,::2,vx]
plot_volume(ax, tmp, x[::3], y[::2], z[::2], xcut=75.0, cmap='jet')

ax.set_xlabel("X",size=16)
ax.set_ylabel("Y",size=16)
ax.set_zlabel("Z",size=16)
ax.set_xlim( x[0],x[-1] )
ax.set_ylim( y[0],y[-1] )
ax.set_zlim( z[0],z[-1] )
ax.set_aspect('equal')
# plt.title('Magnetic energy (t = %6.1f)' % t, size=20)

# plot
plt.show()

# image file
# plt.savefig('output.png', dpi=144)

# end
