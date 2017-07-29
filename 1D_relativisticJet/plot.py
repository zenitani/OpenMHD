import numpy as np
import matplotlib.pyplot as plt
# dummy index
x=0;y=1;ro=2;pr=3;vx=4;vy=5;vz=6;ux=7;uy=8;uz=9;u0=10;bx=11;by=12;bz=13;ex=14;ey=15;ez=16;pt=17

# reading the data...
d1 = np.loadtxt('data/x-00000.dat')
d2 = np.loadtxt('data/x-00002.dat')

line1, = plt.plot(d1[:,x],d1[:,u0],linestyle=':',color='gray')
line2, = plt.plot(d2[:,x],d2[:,u0],linestyle='-',color='black')

plt.xlim( -0.1, 0.1)
plt.xlabel(r'$X$',fontsize=16)
plt.ylabel(r'$\gamma$', fontsize=16)

plt.legend( (line1, line2), ('t=0.0', 't=0.2'), loc='upper left', fontsize=14, shadow=True)

#plt.tight_layout() # if necessary
plt.show()
#plt.savefig("output.png")
