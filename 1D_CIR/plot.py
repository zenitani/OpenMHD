import numpy as np
import matplotlib.pyplot as plt
# dummy index
x=0;y=1;ro=2;pr=3;vx=4;vy=5;vz=6;bx=7;by=8;bz=9

data1 = np.loadtxt('data/x-00000.dat')
data2 = np.loadtxt('data/x-00006.dat')
line1, = plt.plot(data1[:,x],data1[:,ro],linestyle=':',color='gray')
line2, = plt.plot(data2[:,x],data2[:,ro],linestyle='-',color='black')
plt.xlabel('X',fontsize=16)
plt.ylabel('Density', fontsize=16)
#plt.ylabel(r'$\rho$', fontsize=20)
plt.legend( (line1, line2), ('t=000', 't=600'), loc='upper left', shadow=True)
#plt.legend( (line1, line2), ('t=000', 't=600'), loc='best', shadow=True)

plt.show()
#plt.savefig("output.png")
