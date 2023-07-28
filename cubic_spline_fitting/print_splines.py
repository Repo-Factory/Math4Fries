import numpy as np
from matplotlib import pyplot as plt

input = np.loadtxt('./pogin.txt', delimiter=',')
input = np.transpose(input)

spline = np.loadtxt('./pogout.txt', delimiter=',')
spline = np.transpose(spline)


import matplotlib.colors as mcolors
colorItems = mcolors.TABLEAU_COLORS.items()
colorList = list(colorItems)

plt.rcParams["figure.figsize"] = [14.50, 7.50]
plt.rcParams["figure.autolayout"] = True

plt.scatter(input[0], input[1], color='blue',marker="o")

plt.plot(spline[0], spline[1], color=colorList[1%10][0])

plt.xlabel('x')
plt.ylabel('y')
plt.title('Cubic Spline Interpolation')
plt.grid()

plt.show()