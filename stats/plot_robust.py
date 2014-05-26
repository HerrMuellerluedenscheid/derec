import sys 
import matplotlib.pyplot as plt


file_name = sys.argv[1]

f = open(file_name, 'r')
results = []

for l in f.readlines():
    results.append(map(float, l.split()))

fig = plt.figure()
ax = fig.add_subplot(111)

for a,b,c,d in results:
    if d==1.:
        c = 'bo'
    else:
        c = 'ro'
    ax.plot(a/1000.,b,c)
plt.xlabel('Mislocation [km]')
plt.ylabel('Angle [deg]')

plt.show()
