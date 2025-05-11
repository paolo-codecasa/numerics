import numpy as np
import matplotlib.pyplot as plt


data1 = np.loadtxt("schwingereuler.csv",delimiter=",")
data2 = np.loadtxt("schwingerrk4.csv",delimiter=",")
data3 = np.loadtxt("schwingereulerlog.csv",delimiter=",")
data4 = np.loadtxt("schwingerrk4log.csv",delimiter=",")
data5 = np.loadtxt("Trajektorie1.csv",delimiter=",")
data6 = np.loadtxt("Trajektorie2.csv",delimiter=",")
data7 = np.loadtxt("Trajektorie3.csv",delimiter=",")


##############################################################

plt.figure(figsize=(10,6))

# gridlines erzeugen (wie man sieht, gibt es zig Möglichkeiten diese zu Zeichnen)
plt.minorticks_on()
plt.grid(visible=True, which='major', color='#666666', linestyle='-', alpha=0.6)
plt.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.3)

# Achsenbeschriftungen
plt.xlabel(r"t", fontsize="14")
plt.ylabel("x", fontsize="14")

# Titel für den Plot
k = 10;
plt.plot(data1[:,0], data1[:,1],color ='green',linewidth=1,label = "Schwinger numerisch euler")
plt.plot(data1[:,0],-0.5*np.cos(np.sqrt(2*k)*data1[:,0]), color='red', linewidth=1,label = "Schwinger analytisch")


plt.legend(loc = 1,fontsize ="10")

plt.savefig("schwingereuler.pdf")



######################################################
plt.figure(figsize=(10,6))

# gridlines erzeugen (wie man sieht, gibt es zig Möglichkeiten diese zu Zeichnen)
plt.minorticks_on()
plt.grid(visible=True, which='major', color='#666666', linestyle='-', alpha=0.6)
plt.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.3)

# Achsenbeschriftungen
plt.xlabel(r"t", fontsize="14")
plt.ylabel("x", fontsize="14")


plt.plot(data2[:,0], data2[:,1], color='green', linewidth=1,label = "Schwinger numerisch rk4")
plt.plot(data2[:,0],-0.5*np.cos(np.sqrt(2*k)*data2[:,0]), color='red', linewidth=1,label = "Schwinger analytisch")



plt.legend(loc = 1,fontsize ="10")

plt.savefig("schwingerrk4.pdf")




###########################################################
plt.figure(figsize=(10,6))

# gridlines erzeugen (wie man sieht, gibt es zig Möglichkeiten diese zu Zeichnen)
plt.minorticks_on()
plt.grid(visible=True, which='major', color='#666666', linestyle='-', alpha=0.6)
plt.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.3)

# Achsenbeschriftungen
plt.ylabel(r"$\sigma_x$", fontsize="14")
plt.xlabel(r"$\Delta t$", fontsize="14")


plt.plot(data3[:,0], data3[:,1],'g.',label = "Schwinger Fehler numerisch eulerlogarithmisch")
plt.plot(data4[:,0], data4[:,1],'r.',label = "Schwinger Fehler numerisch rk4logarithmisch")


plt.xscale('log')
plt.yscale('log')
plt.legend(loc = 1,fontsize ="10")

plt.savefig("Fehler der Integrationsverfahren.pdf")
###################################################################################################
plt.figure(figsize=(10,6))

# gridlines erzeugen (wie man sieht, gibt es zig Möglichkeiten diese zu Zeichnen)
plt.minorticks_on()
plt.grid(visible=True, which='major', color='#666666', linestyle='-', alpha=0.6)
plt.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.3)


# Achsenbeschriftungen
plt.xlabel(r"x", fontsize="14")
plt.ylabel("y", fontsize="14")


plt.plot(data5[:,0], data5[:,1], color='green', linewidth=1,label = "Trajektorie mit x =-0.5, y= 0.5,v_x = 0, v_y = 0")




plt.legend(loc = 1,fontsize ="10")

plt.savefig("Trajektorie1.pdf")
###################################################################################################
plt.figure(figsize=(10,6))

# gridlines erzeugen (wie man sieht, gibt es zig Möglichkeiten diese zu Zeichnen)
plt.minorticks_on()
plt.grid(visible=True, which='major', color='#666666', linestyle='-', alpha=0.6)
plt.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.3)


# Achsenbeschriftungen
plt.xlabel(r"x", fontsize="14")
plt.ylabel("y", fontsize="14")


plt.plot(data6[:,0], data6[:,1], color='green', linewidth=1,label = "Trajektorie mit x =-0.1, y= 0.5,v_x = 0, v_y = 0")




plt.legend(loc = 1,fontsize ="10")

plt.savefig("Trajektorie2.pdf")
###################################################################################################
plt.figure(figsize=(10,6))

# gridlines erzeugen (wie man sieht, gibt es zig Möglichkeiten diese zu Zeichnen)
plt.minorticks_on()
plt.grid(visible=True, which='major', color='#666666', linestyle='-', alpha=0.6)
plt.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.3)


# Achsenbeschriftungen
plt.xlabel(r"x", fontsize="14")
plt.ylabel("y", fontsize="14")


plt.plot(data7[:,0], data7[:,1], color='green', linewidth=1,label = "Trajektorie mit x =-0.1, y= 0.1,v_x = 0, v_y = 3")




plt.legend(loc = 1,fontsize ="10")
plt.show()
plt.savefig("Trajektorie3.pdf")
###################################################################################################