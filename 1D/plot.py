import numpy as np
import matplotlib.pyplot as plt
# dummy index
x=0;y=1;ro=2;pr=3;vx=4;vy=5;vz=6;bx=7;by=8;bz=9

data = np.loadtxt('data/x-00002.dat')
plt.plot(data[:,x],data[:,pr])
plt.xlabel('X')
plt.ylabel('Pressure')
plt.show()
#plt.savefig("output.png")
