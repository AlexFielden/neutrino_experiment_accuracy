# -*- coding: utf-8 -*-
from matplotlib import pyplot as plt #imports the matplotlib library
import numpy as np #imports the numpy library
L = 53 # sets up variables for the baseline distance and the mixing parameters theta 12 and 13
t12 = 0.583638
t13 = 0.1496

y = np.linspace(0.000,0.490,1000) # this section of the code produces the energy spectrum of the neutrinos aswell as the expected number of neutrinos for that energy
y = np.append(y, np.linspace(0.490,0.742,1000))
y = np.append(y, np.linspace(0.742,0.699,1000))
y = np.append(y, np.linspace(0.699,0.500,1000))
y = np.append(y, np.linspace(0.500,0.321,1000))
y = np.append(y, np.linspace(0.321,0.167,1000))
y = np.append(y, np.linspace(0.167,0.066,1000))
y = np.append(y, np.linspace(0.066,0.007,1000))
x = np.linspace(0.0018,0.0026,1000)
x = np.append(x, np.linspace(0.0026,0.0034,1000))
x = np.append(x, np.linspace(0.0034,0.0042,1000))
x = np.append(x, np.linspace(0.0042,0.005,1000))
x = np.append(x, np.linspace(0.005,0.0058,1000))
x = np.append(x, np.linspace(0.0058,0.0066,1000))
x = np.append(x, np.linspace(0.0066,0.0074,1000))
x = np.append(x, np.linspace(0.0074,0.0082,1000))
numneutrinos = np.array([x,y])
plt.figure(1) #produces plot for the energy spectrum
plt.plot(numneutrinos[0],numneutrinos[1])
plt.xlabel("Neutrino energy / Gev")
plt.ylabel("number of neutrinos / A.U")
En = numneutrinos[0]

def energyvals(E,res): #function to add the noise to the energy values to a set amount
    y = np.empty_like(E)
    x = (res/100)
    change = np.random.uniform(low=(1-x),high=(1+x),size=8000)
    for i in range(0,len(E)):
        y[i] = ( E[i] * change[i])
    return y


d21 = 7.42E-5 * 1.27 * (L/(En)) # declares the variables for mass difference that do not change
d32 = 2.514E-3 * 1.27 * (L/(En))
d31 = d32 + d21 # declares the variables for normal mass hierarchy

pN = 1 - ((np.sin(2*t13))**2) * ( ((np.cos(t12))**2) * ((np.sin(d31))**2) + ((np.sin(t12))**2) * ((np.sin(d32))**2) ) - ((np.cos(t13))**4) * (((np.sin(2*t12))**2)) * ((np.sin(d21))**2)
#calculates the survival probability for normal hierarchy


d31 = d32 - d21 # declares the variables for inverse mass hierarchy

pI = 1 - ((np.sin(2*t13))**2) * ( ((np.cos(t12))**2) * ((np.sin(d31))**2) + ((np.sin(t12))**2) * ((np.sin(d31))**2) ) - ((np.cos(t13))**4) * (((np.sin(2*t12))**2)) * ((np.sin(d21))**2)
##calculates the survival probability for inverse hierarchy
yN = numneutrinos[1] * pN # calcuclates the expected number of neutrinos to be detected for each hierarchy ordering
yI = numneutrinos[1] * pI
plt.figure(2) #produces a plot of the survival probabilties
plt.plot(En,pI, En, pN)
plt.xlabel("Energy / GeV")
plt.ylabel("probability of survival / A.U")
plt.legend(["Inverse hierarchy","Normal hierarchy"])
plt.title("Probability of survival for electron antineutrino vs energy of the neutrino")
plt.figure(3) #produces a graph of the expected number of neutrinos for each hierarchy
plt.xlabel("Energy / GeV")
plt.ylabel("Number of electron antineutrinos / A.U")

plt.plot(x,yN,x,y,x,yI)
plt.legend(["Neutrino expectation for normal hierarchy","Neutrino expectation for no oscillation","Neutrino expectation for inverse hierarchy"],loc=1)
plt.show()



En = energyvals(numneutrinos[0],9) #calls function to add noise to energy values to a resolution of 9%
d21 = 7.42E-5 * 1.27 * (L/(En))
d32 = 2.514E-3 * 1.27 * (L/(En))
d31 = d32 + d21
pNN = 1 - ((np.sin(2*t13))**2) * ( ((np.cos(t12))**2) * ((np.sin(d31))**2) + ((np.sin(t12))**2) * ((np.sin(d32))**2) ) - ((np.cos(t13))**4) * (((np.sin(2*t12))**2)) * ((np.sin(d21))**2)
#this section of the code recalculates the expected number of neutrinos for each energy after noise has been added for both hierarchies
d31 = d32 - d21
pIN = 1 - ((np.sin(2*t13))**2) * ( ((np.cos(t12))**2) * ((np.sin(d31))**2) + ((np.sin(t12))**2) * ((np.sin(d32))**2) ) - ((np.cos(t13))**4) * (((np.sin(2*t12))**2)) * ((np.sin(d21))**2)
yNN = numneutrinos[1] * pNN
yIN = numneutrinos[1] * pIN

plt.figure(4) #this section of the code produces the graphs to compare the noise added plots with the expected values
plt.plot(En,yNN,En,yIN,zorder=1)
plt.plot(numneutrinos[0],yN,zorder=2)
plt.xlabel("Energy / GeV")
plt.ylabel("Number of neutrinos / A.U")
plt.legend(["Normal hierarchy with added noise","Inverse hierarchy with added noise","Expected number of neutrinos for normal hierarchy"])
plt.title("Noise added simulations with normal hierarchy, Number of neutrinos vs neutrino energy")
plt.show()
plt.figure(7)
plt.plot(En,yIN,En,yNN,zorder=1)
plt.plot(numneutrinos[0],yI,zorder=2)
plt.xlabel("Energy / GeV")
plt.ylabel("Number of neutrinos / A.U")
plt.legend(["Inverse hierarchy with added noise","Normal hierarchy with added noise","Expected number of neutrinos for inverse hierarchy"])
plt.title("Noise added simulations with Inverse hierarchy, Number of neutrinos vs neutrino energy")
plt.show()

NormaldiffnormalN = yN - yNN # calculates the difference as shown in the method
normaldiffinverseN = yN - yIN

inversediffnormalN = yI - yNN
inversediffinverseN = yI - yIN

plt.figure(5) #this section of the code produces graphs for the modulus of the difference.
plt.plot(En,abs(NormaldiffnormalN),En, abs(normaldiffinverseN))
plt.legend(["Modulus difference of normal and normal with noise","Modulus difference of normal and inverse with noise"])
plt.xlabel("Energy / GeV")
plt.ylabel("Modulus of differnece / A.U")
plt.figure(6)
plt.plot(En, abs(inversediffnormalN),zorder=2)
plt.plot(En,abs(inversediffinverseN),zorder=1)
plt.xlabel("Energy / GeV")
plt.ylabel("Modulus of differnece / A.U")
plt.legend(["Modulus difference of inverse and normal with noise","Modulus difference of inverse and inverse with noise"])
plt.show

