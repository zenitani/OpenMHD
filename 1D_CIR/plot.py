import numpy as np
import matplotlib.pyplot as plt
# dummy index
x=0;y=1;ro=2;pr=3;vx=4;vy=5;vz=6;bx=7;by=8;bz=9

# reading data...
d1 = np.loadtxt('data/x-00000.dat')
d2 = np.loadtxt('data/x-00006.dat')

# Y axis (density)
#line1, = plt.plot(d1[:,x],d1[:,ro],linestyle=':',color='gray')
#line2, = plt.plot(d2[:,x],d2[:,ro],linestyle='-',color='black')
#plt.ylabel('Density', fontsize=16)
#plt.ylabel(r'$\rho$', fontsize=20)
# Y axis (|B|)
line1, = plt.plot(d1[:,x],np.sqrt(d1[:,bx]**2+d1[:,by]**2+d1[:,bz]**2),linestyle=':',color='gray')
line2, = plt.plot(d2[:,x],np.sqrt(d2[:,bx]**2+d2[:,by]**2+d2[:,bz]**2),linestyle='-',color='black')
plt.ylabel(r'$|B|$', fontsize=16)

# X axis
plt.xlim( 0, 10000)
plt.xlabel(r'$X$',fontsize=16)
plt.legend( (line1, line2), ('t=000', 't=600'), loc='upper left', shadow=True)
#plt.legend( (line1, line2), ('t=000', 't=600'), loc='best', shadow=True)

#plt.tight_layout() # if necessary
plt.show()
#plt.savefig("output.png")
