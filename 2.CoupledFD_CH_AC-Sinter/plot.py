import numpy as np
import matplotlib.pyplot as plt

# For plotting curves using TWO particles
#d1 = np.loadtxt("Results/2-particles/energy-time.txt")
#d2 = np.loadtxt("Results/2-particles/energy-step.txt")

# For plotting curves using NINE particles
d1 = np.loadtxt("energy-time.txt")
d2 = np.loadtxt("energy-step.txt")

#plt.ylable("Free Energy")
plt.figure(1)
plt.subplot(211)
plt.plot(d1[:,0],d1[:,1],'b-')

plt.subplot(212)
plt.plot(d2[:,0],d2[:,1],'r-')

plt.show()
