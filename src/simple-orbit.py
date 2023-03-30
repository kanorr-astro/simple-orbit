#Simple Orbit

import numpy as np
import math as math
import matplotlib.pyplot as plt
import utility as util

#Body Parameters
mu = 1                              #DU^3 / TU^2

#Initial positions/velo
r = np.array([0, 1, 0])                       #DU
v = np.array([0, 0, 0.5])                     #DU / TU

#calculate angular momentum (h)
h = np.cross(r,v)

#calculate line of nodes(n):
n = np.cross([0,0,1],h)

#calculate e vector:
term1 = (util.magnitude(v)**2 - (mu / util.magnitude(r))) * r
term2 = np.dot(r,v) * v
e_vect = (1/mu) * (term1 - term2)


#calculate specific mechanical energy
SME = util.spec_mech_energy(r,v,mu)

#calculate semi-major axis (a)
if SME == 0:
    a = "inf"
else:    
    a = - mu / (2*SME)

#calculate semi-latus rectum (p):    
p = util.magnitude(h)**2 / mu

#calculate eccentricity (e):
if a == "inf":
    e = 1
else:
    e = math.sqrt(1 - (p / a))

#calculate orbit type:
if util.magnitude(h) > 0:
    if e == 0:
        conic = "circular"
    elif (e > 0) and (e < 1):
        conic = "ellipse"
    elif e == 1:
        conic = "parabolic"
    elif e > 1:
        conic = "hyperbolic"
    else:
        conic = "unknown"
else:
    conic = "degenerate"

#calculate inclination (i):
i = math.acos(h[2]/util.magnitude(h))

#calculate longitude of the ascending node:
long_asc = math.acos(n[0]/util.magnitude(n))

#calculate argument of periapsis:
arg_peri = math.acos(np.dot(n,e_vect) / (util.magnitude(n) * util.magnitude(e_vect)))

print("Angular Momentum Vector(h): ", h)
print("Eccentricity Vector: ", e_vect)
print("Line of Nodes Vector: ", n)
print("Specific Mechanical Energy: ", SME)
print("Semi-Major Axis(a): ", a)
print("Semi-Latus Rectum(p): ", p)
print("Eccentricity(e): ", e)
print("Orbit Type: ", conic)
print("Inclination(i): ", math.degrees(i))
print("Long. of Ascending Node: ", long_asc)
print("Argument of Periapsis: ", math.degrees(arg_peri))

nu_raw = np.linspace(0, 2*math.pi, 100)
nu_mod = []
if e >= 1:
    for nu in nu_raw:
        if (nu > (7*math.pi / 8)) and (nu < (9*math.pi / 8)):
            None
        else:
            nu_mod.append(nu)
else:
    nu_mod = nu_raw

rad = []
for nu in nu_mod:
    rad.append(util.radius(p,e,nu))


fig,ax = plt.subplots(subplot_kw={'projection':'polar'})
ax.plot(nu_mod, rad)
ax.set_rmax(np.max(rad))
ax.set_rticks(np.linspace(0,round(np.max(rad)+1)),5)
ax.grid(True)

ax.set_title("Orbital Diagram in PQW Frame")
plt.show()


