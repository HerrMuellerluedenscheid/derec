import numpy as num
import matplotlib.pyplot as plt



freqs = num.fromfile("freqs.txt")
spec = num.fromfile('spectr.txt')

print freqs.shape
print spec.shape
print spec

plt.plot(spec)
plt.show()

