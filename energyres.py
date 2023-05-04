from matplotlib import pyplot as plt
import numpy as np




def energyvals(E,res):
    y = np.empty_like(E)
    x = (res/100)
    change = np.random.uniform(low=(1-x),high=(1+x),size=100000)
   
    print(change)
    #print(y)
    for i in range(0,len(E)):
        y[i] = ( E[i] * change[i])
    return y

E= 0.004
L=55
s13 = 2.34E-2
s12 = 3.08E-1
S23 = 4.37E-1
m21 = 7.54E-5
m31 = 2.47E-3
c13 = 1 - s13
c12 = 1 - s12
c413 = 1-(2*s13) +(s13)**2
d21 = 7.54E-5 * 1.27 * (L/(E)) """add hbar and c product swap 4 to 1.27 in numerator"""
d32 = 2.47E-3 * 1.27 * (L/(E))
mee = (np.cos(t12)**2)*(m31)+((np.sin(t12))**2)*m31
sphi = (c12*np.sin(2*s12*d21)-s12*np.sin(2*c12*d21))/((1-((np.sin(2*0.58833))**2)*((np.sin(d21))**2))**(1/2))
phi = np.arcsin(sphi)



#P = 1 - 2* s13 * c13 * - 4 * c413 * s12 * c12 * (np.sin(d21)**2) + 2 * s13* c13 * (1 - (4* s12 * c12 * (np.sin(d21))**2))**(1/2) * np.cos(2*d32*phi)



P = 1 - (1/2)*((np.sin(2*0.15357))**2) * (1-((1-(np.sin(2*0.58833)**2)*(np.sin(d21)**2))**(1/2))*np.cos(2*mee*phi)) 


c12*np.sin(2*s12*d21)-s12*np.sin(2*c12*d21)



"""

t12 = 33.44
t13 = 8.57
d31 = 2.47E-3  * (L/(4*E))
m31 = 2.47E-3
d21 = 7.54E-5  * (L/(4*E))
m21 = 7.54E-5 


sphi = (c12*np.sin(2*s12*d21)-s12*np.sin(2*c12*d21))/((1-2*s12*c12*((np.sin(d21))**2))**(1/2))
cphi = (c12*np.cos(2*s12*d21)-s12*np.cos(2*c12*d21))/((1-2*s12*c12*(np.sin(d21))**2)**(1/2))
mee = (np.cos(t12)**2)*(m31)+((np.sin(t12))**2)*m31
s1 = (-1/2)*((2*s13*c13**2)(1-(1-(np.sin(t12*2)**2)*(np.sin(d21))**2)**(1/2)))
s2 = -((np.cos(t13))**4)*(np.cos(2*t12))**2(np.sin(d21))**2


plt.figure(1)
plt.plot(P,E)
plt





"""