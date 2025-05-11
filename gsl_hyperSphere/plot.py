import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt # Bibliothek für Optimierungsprobleme (lineare Regression)


# 
data1_x = np.loadtxt("sphere_num.csv", delimiter=',', usecols=0) #Dimension
data1_y = np.loadtxt("sphere_num.csv", delimiter=',', usecols=1) #Integral
data2_x = np.loadtxt("sphere_ana.csv", delimiter=',', usecols=0)
data2_y = np.loadtxt("sphere_ana.csv", delimiter=',', usecols=1)

plt.figure(figsize=(10,6))

plt.scatter(data1_x, data1_y, c="g", label = "numerisch")
plt.scatter(data2_x, data2_y, c="r", marker="x", label = "analytisch")
plt.xlabel(r"D", fontsize="14")
plt.ylabel("Integral", fontsize="14")

plt.savefig("dim.png", dpi=200)

# 
data3_x = np.loadtxt("logError_montecarlo.csv", delimiter=',', usecols=0) #N
data3_y = np.loadtxt("logError_montecarlo.csv", delimiter=',', usecols=1) #Abweichung

def f(x, a):
    return a * 1/np.sqrt(x)

plt.figure(figsize=(10,6))

param, cov = opt.curve_fit(f, data3_x, data3_y)
print("a:", param[0])

y1 = f(data3_x, *param)

plt.scatter(data3_x, data3_y)
plt.plot(data3_x, y1,'r.', linewidth=1,label = "theoretischer Abfall")
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"N", fontsize="14")
plt.ylabel("Integrationsfehler", fontsize="14")

plt.savefig("error.png", dpi=200)


# 
data4_x = np.loadtxt("standardDeviation_montecarlo.csv", delimiter=',', usecols=0) #N
data4_y = np.loadtxt("standardDeviation_montecarlo.csv", delimiter=',', usecols=1) #Standardabweichung

plt.figure(figsize=(10,6))

param, cov = opt.curve_fit(f, data4_x, data4_y)
print("a:", param[0])

y2 = f(data4_x, *param)

plt.scatter(data4_x, data4_y)
plt.plot(data4_x, y2,'r.', linewidth=1,label = "theoretischer Abfall")

plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"N", fontsize="14")
plt.ylabel(r"Standardabweichung $\sigma$", fontsize="14")

plt.savefig("standDev.png", dpi=200)


plt.show()
