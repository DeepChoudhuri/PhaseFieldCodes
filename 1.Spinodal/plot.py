import numpy as np
import matplotlib.pyplot as plt

#d1 = np.loadtxt("energy-time.txt")
d1 = np.loadtxt("energy-time-FFT.txt")

plt.plot(d1[:,0],d1[:,1])

plt.show()
