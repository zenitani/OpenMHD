import matplotlib.pyplot as plt
import openmhd
# dummy index
vx=0;vy=1;vz=2;pr=3;ro=4;bx=5;by=6;bz=7;ps=8

# reading the data ...
x,y,t,data = openmhd.data_read(25)

plt.clf()

# 2D image
# extent: [left, right, bottom, top]
plt.imshow(data[:,:,ro].T,origin='lower',extent=[x[0],x[-1],y[0],y[-1]])

# useful options
# plt.grid()
plt.xlabel("X",size=16)
plt.ylabel("Y",size=16)
plt.title('Density (t = %6.1f)' % t, size=20)
plt.colorbar()

# plot
plt.show()

# image file
#plt.savefig('output.png', dpi=144)

# end
