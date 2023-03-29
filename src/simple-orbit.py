#Simple Orbit

import numpy as np
import math as math
import matplotlib.pyplot as plt
import utility as util

#Calculate and plot simple orbit for 2-body w/ body 1 at origin w/ mass m1 and body 2 orbiting w/ mass m2


#Body Parameters
mu = 1                              #DU^3 / TU^2

#Initial positions/velo
r = [4, 0, 0]                       #DU
v = [0, 0.4, 0]                       #DU / TU


h = np.cross(r,v)
SME = util.spec_mech_energy(r,v,mu)
a = - mu / (2*SME)
p = util.magnitude(h)**2 / mu
e = math.sqrt(1 - (p / a))



print("Angular Momentum Vector: ", h)
print("Specific Mechanical Energy: ", SME)
print("Semi-Major Axis: ", a)
print("Semi-Latus Rectum: ", p)
print("Eccentricity: ", e)


nu_vals = np.linspace(0, 2*math.pi, 100)
rad = []
for nu in nu_vals:
    rad.append(util.radius(p,e,nu))



fig,ax = plt.subplots(subplot_kw={'projection':'polar'})
ax.plot(nu_vals,rad)
ax.set_rmax(np.max(rad))
ax.set_rticks(np.linspace(0,round(np.max(rad)+1)),5)
ax.grid(True)

ax.set_title("test")
plt.show()


