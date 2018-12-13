import matplotlib.pyplot as plt
from statistics import median, variance, mean
from math import sqrt
import numpy as np

x = [i/100 for i in range(10000)]

y1 = [sqrt(0.5 * i**3) for i in x]
y2 = [sqrt(1 * i**3) for i in x]
y3 = [sqrt(1 * i) for i in x]

plt.figure(1)
plt.plot(x, y1, x, y2, x, y3)
plt.grid(True)

y4 = [sqrt(1 * i + 0) for i in x]
y5 = [sqrt(1 * i + 5) for i in x]
y6 = [sqrt(1 * i + 10) for i in x]

plt.figure(2)
plt.plot(x, y4, x, y5, x, y6)
plt.grid(True)
plt.show()