import matplotlib.pyplot as plt
import numpy as np


###################################### 
data1_x = np.loadtxt("analytisch.csv", delimiter=',', usecols=0) #n
data1_y = np.loadtxt("analytisch.csv", delimiter=',', usecols=1) #Eigenwert
data2_x = np.loadtxt("numerisch.csv", delimiter=',', usecols=0)
data2_y = np.loadtxt("numerisch.csv", delimiter=',', usecols=1)

plt.figure(figsize=(10,6))

plt.scatter(data1_x, data1_y, c="g", marker="x", label = "analytisch")
plt.scatter(data2_x, data2_y, c="r", marker="+", label = "numerisch")
plt.xlabel("n", fontsize="14")
plt.ylabel("Eigenwerte ohne Potential", fontsize="14")
plt.legend() #******************

plt.savefig("Vergleich noPot.png", dpi=200)


plt.figure(figsize=(10,6))

plt.scatter(data1_x, abs(data1_y-data2_y), c="g", marker="x", label = "Abweichung")
plt.xlabel("n", fontsize="14")
plt.ylabel("Abweichung Eigenwerte noPot", fontsize="14")
plt.legend() #******************

plt.savefig("Abweichung noPot.png", dpi=200)



data3_x = np.loadtxt("EigenVektoren.csv", delimiter=',', usecols=0) #x
data3_y1 = np.loadtxt("EigenVektoren.csv", delimiter=',', usecols=1) #psi(x, n=1)
data3_y2 = np.loadtxt("EigenVektoren.csv", delimiter=',', usecols=2) #psi(x, n=2)
data3_y3 = np.loadtxt("EigenVektoren.csv", delimiter=',', usecols=3) #psi(x, n=3)

plt.figure(figsize=(10,6))

plt.scatter(data3_x, data3_y1)
plt.scatter(data3_x, data3_y2)
plt.scatter(data3_x, data3_y3)
plt.xlabel("x", fontsize="14")
plt.ylabel(r"Eigenvektoren $\psi(x, n)$ ohne Potential", fontsize="14")

plt.savefig("EigenVektoren.png", dpi=200)

###################################### 
data4_x = np.loadtxt("analytisch_HARM.csv", delimiter=',', usecols=0) #n
data4_y = np.loadtxt("analytisch_HARM.csv", delimiter=',', usecols=1) #Eigenwert
data5_x = np.loadtxt("numerisch_HARM.csv", delimiter=',', usecols=0)
data5_y = np.loadtxt("numerisch_HARM.csv", delimiter=',', usecols=1)

plt.figure(figsize=(10,6))

plt.scatter(data4_x, data4_y, c="g", marker="x", label = "analytisch")
plt.scatter(data5_x, data5_y, c="r", marker="+", label = "numerisch")
plt.xlabel("n", fontsize="14")
plt.ylabel(r"Eigenwerte für $V(x) = \frac{k}{2} x^2$", fontsize="14")
plt.legend() #******************

plt.savefig("Vergleich harmPot.png", dpi=200)


plt.figure(figsize=(10,6))

plt.scatter(data4_x, abs(data4_y-data5_y), c="g", marker="x", label = "Abweichung")
plt.xlabel("n", fontsize="14")
plt.ylabel(r"Abweichung Eigenwerte für $V(x) = \frac{k}{2} x^2$", fontsize="14")
plt.legend() #******************

plt.savefig("Abweichung harmPot.png", dpi=200)

plt.show()
