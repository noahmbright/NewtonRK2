import numpy as np
from matplotlib import pyplot as plt

#Define parameters
G= 66425.33725 #m^3/kg/yr^2
m1 = 5.972e24 #kg
m2 = 1.989e30 #kg
mu = m1*m2/(m1+m2) #kg
r0 = 151.6e9 #m
phidot0 = np.pi*2  #rad/yr
l = phidot0*mu*(r0**2)
gamma = G*m1*m2
params = gamma*mu/(l**2)

t0 = 0
tf= 1
n=5000
deltat=(tf-t0)/(n-1)

t=np.linspace(t0,tf,n)

k1 = np.zeros([n])
k2 = np.zeros([n])
l1 = np.zeros([n])
l2 =np.zeros([n])
j1=np.zeros([n])
j2=np.zeros([n])

r = np.zeros([n])
dr = np.zeros([n])
d2r = np.zeros([n])

phi=np.zeros([n])

r[0] = r0
dr[0] = 0
d2r = (-gamma/r**2 + l**2/mu/r**3)/mu

phi[0] = 0



for i in range(1,n):
    l1[i-1] = deltat * d2r[i-1]
    l2[i-1] = deltat * (d2r[i-1] + l1[i-1])
    dr[i] = dr[i - 1] + (l1[i - 1] + l2[i - 1])/2

    k1[i - 1] = deltat * dr[i - 1]
    k2[i - 1] = deltat * (dr[i - 1] + k1[i - 1])
    r[i] = r[i-1] +(k1[i-1]+k2[i-1])/2
    d2r[i] = (-gamma/r[i]**2 + l**2/mu/r[i]**3)/mu

    j1[i-1] = deltat * l/mu/r[i-1]**2
    j2[i-1] = deltat * (l/mu/(r[i-1]+j1[i-1])**2)
    phi[i] = phi[i-1] + (j1[i-1]+j2[i-1])/2





graph = plt.subplot(111, projection = "polar")
graph.plot(phi,r)
plt.show()

print(phi[n-1])