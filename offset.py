import numpy as num
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms

Trans = transforms.Transform()
Trans.transform((3,4))
x = num.arange(0,10, 0.1)

y1 = num.sin(x)
y2 = num.sin(x*0.5)
y3 = num.sin(x*1.5)

fig = plt.figure(facecolor='white')

ax = fig.add_subplot(111)
ax.axes.get_yaxis().set_visible(False)

dx = 0
dy = 40/72.
offset = transforms.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
trans = ax.transData + offset
print trans

line1 = plt.plot(x, y1)
line2 = plt.plot(x, y2, transform=trans)
line3 = plt.plot(x, y3)

plt.title('I\'d like to have this...')



plt.show()

