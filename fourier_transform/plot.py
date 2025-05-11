import matplotlib.pyplot as plt
import numpy as np

#plt.rcParams.update({'font.size': 10})

data10 = np.loadtxt("data10.csv", delimiter=",")
data30 = np.loadtxt("data30.csv", delimiter=",")
data100 = np.loadtxt("data100.csv", delimiter=",")
data500 = np.loadtxt("data500.csv", delimiter=",")

fig, axs = plt.subplots(2, 2, figsize=(12, 8), sharex="col")
axs[0, 0].plot(data10[:,0], data10[:,1])
axs[0, 0].set_title('M=10')
axs[0, 1].plot(data30[:,0], data30[:,1], 'tab:orange')
axs[0, 1].set_title('M=30')
axs[1, 0].plot(data100[:,0], data100[:,1], 'tab:green')
axs[1, 0].set_title('M=100')
axs[1, 1].plot(data500[:,0], data500[:,1], 'tab:red')
axs[1, 1].set_title('M=500')
#axs.label.set_size(10)
#axs.set_units

for ax in axs.flat:
    ax.set(xlabel='k', ylabel='Ff(k)')


# Hide x labels and tick labels for top plots and y ticks for right plots.
#for ax in axs.flat:
 #   ax.label_outer()
#DON'T HIDE

plt.savefig('plot.pdf')
plt.show()

quit()