import numpy as np
import matplotlib.pyplot as plt

i = np.arange(1, 2000)
alpha = 1.5

a = np.abs(i+2)**(3-alpha) - 3*np.abs(i+1)**(3-alpha) + 3*np.abs(i)**(3-alpha) - np.abs(i-1)**(3-alpha)
b = -a * i**(alpha)
plt.plot(i, b)
plt.show()
