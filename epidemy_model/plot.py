import matplotlib.pyplot as plt
import numpy as np

data10 = np.loadtxt("data.csv", delimiter=",")
#"Zeit, Susceptible, Exposed, Infectious, Removed" alle abhängig von beta -->13 spalten

fig, axs = plt.subplots(3)
axs[0].plot(data[:,0], data[:,1], 'tab:blue', label="Susceptible")
axs[0].plot(data[:,0], data[:,4], 'tab:orange', label="Exposed")
axs[0].plot(data[:,0], data[:,7], 'tab:green', label="Infectious")
axs[0].plot(data[:,0], data[:,10], 'tab:black', label="Removed")
axs[0].set_title('R_0=1.25')
axs[1].plot(data[:,0], data[:,2], 'tab:blue')
axs[1].plot(data[:,0], data[:,5], 'tab:orange')
axs[1].plot(data[:,0], data[:,8], 'tab:green')
axs[1].plot(data[:,0], data[:,11], 'tab:black')
axs[1].set_title('R_0=1.5')
axs[2].plot(data[:,0], data[:,3], 'tab:blue')
axs[2].plot(data[:,0], data[:,6], 'tab:orange')
axs[2].plot(data[:,0], data[:,9], 'tab:green')
axs[2].plot(data[:,0], data[:,12], 'tab:black')
axs[2].set_title('R_0=2')
ax.legend()

for ax in axs.flat:
    ax.set(xlabel='k', ylabel='Ff(k)')

# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()

plt.savefig('plot.pdf')
plt.show()

quit()