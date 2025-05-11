import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("temp.csv", delimiter=",", skiprows=1)

temp_data = np.delete(data, 0, axis=1)



print(temp_data)

plt.imshow(temp_data, aspect='auto')
plt.show()
plt.savefig("temp.png", dpi=200)