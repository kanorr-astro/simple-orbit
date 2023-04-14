#Simple Orbit 3d

#Module Imports
import numpy as np
import math as m
import matplotlib.pyplot as plt
import utility as util

#Body Parameters - Canonical Units
mu = 1
ae = 1

#Initial Position/Velocity Vectors - Canonical Units
r = np.array([0,1,1])
v = np.array([-1,0,-0.05])

#calculate angular momentum (h):
h = np.cross(r,v)

#calculate the line of nodes (n):
n = np.cross([0,0,1],h)

#calculate eccentricity vector (e):
term1 = (util.magnitude(v)**2 - (mu / util.magnitude(r))) * r
term2 = np.dot(r,v) * v
e_vect = (1/mu)*(term1 - term2)

#calculate specific mechanical energy (SME):
SME = util.spec_mech_energy(r,v,mu)

#calculate semi-major axis (a):
if SME >= 0:
    a = "inf"
else:
    a = -mu / (2*SME)

#calculate semi-latus rectum (p):
p = util.magnitude(h)**2 / mu

#calculate eccentricity (e):
e = m.sqrt(1 + ((2*SME*util.magnitude(h)**2)/(mu**2)))

#calculate orbit conic type:
if util.magnitude(h) > 0:
    if e == 0:
        conic = "circular"
    elif (e > 0) and (e < 1):
        conic = 'elliptic'
    elif e == 1:
        conic = 'parabolic'
    elif e > 1:
        conic = 'hyperbolic'
else:
    conic = 'degenerate'

#calculate inclination (i):
i = m.acos(h[2]/util.magnitude(h))

#calculate longitude of the ascending node(Omega):
ascend = m.acos(n[0]/util.magnitude(n))

#calculate argument of periapsis (omege):
peri = m.acos(np.dot(n,e_vect) / (util.magnitude(n) * util.magnitude(e_vect)))

#Generate linspace for nu angles
nu_raw = np.linspace(0, 2*m.pi, 100)
nu_mod = []
if e >= 1:
    for nu in nu_raw:
        if (nu > (7*m.pi / 8)) and (nu < (9*m.pi / 8)):
            None
        else:
            nu_mod.append(nu)
else:
    nu_mod = nu_raw

#Using nu angles generate list of points in PQW orientation:
PQW = []
for nu in nu_mod:
    PQW.append(np.array([util.radius(p,e,nu)*m.cos(nu), util.radius(p,e,nu)*m.sin(nu), 0]))



#Break PQW coords into individual coord lists for plotting
pcoord = []
qcoord = []
wcoord = []
for point in PQW:
    pcoord.append(point[0])
    qcoord.append(point[1])
    wcoord.append(point[2])

#define PQW to IJK rotation matrix (R):
R_raw = [m.cos(ascend)*m.cos(peri)-m.sin(ascend)*m.sin(peri)*m.cos(i),
         -m.cos(ascend)*m.sin(peri)-m.sin(ascend)*m.cos(peri)*m.cos(i),
         m.sin(ascend)*m.sin(i),
         m.sin(ascend)*m.cos(peri)+m.cos(ascend)*m.sin(peri)*m.cos(i),
         -m.sin(ascend)*m.sin(peri)+m.cos(ascend)*m.cos(peri)*m.cos(i),
         -m.cos(ascend)*m.sin(i),
         m.sin(peri)*m.sin(i),
         m.cos(peri)*m.sin(i),
         m.cos(i)]
R = np.array(R_raw).reshape(3,3)

#Transform from PQW to IJK coords:
IJK = []
for point in PQW:
    IJK.append(np.dot(R,point))

#Break IJK coords into individual coord lists for plotting
icoord = []
jcoord = []
kcoord = []
for point in IJK:
    icoord.append(point[0])
    jcoord.append(point[1])
    kcoord.append(point[2])

#generate date for wireframe Earth:
u,v =  np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    
px = ae*np.cos(u)*np.sin(v)
py = ae*np.sin(u)*np.sin(v)
pz = ae*np.cos(v)

#Output Orbital Parameters:
print("Orbit Conic Type: ", conic)
print("Semi-Major Axis(a): ", a)
print("Semi-Latus Rectum(p): ", p)
print("Eccentricity(e): ", e)
print("Inclination(i): ", m.degrees(i))
print("Longitude of Ascending Node(Omega): ", m.degrees(ascend))
print("Argument of Periapsis: ", m.degrees(peri))


fig = plt.figure(figsize=(12,6))

ax1 = fig.add_subplot(121, projection = '3d')
ax1.set_xlim(-2, 2)
ax1.set_ylim(-2, 2)
ax1.set_zlim(-2, 2)
ax1.scatter(pcoord, qcoord, wcoord)
ax1.title.set_text('PQW Plot of Orbit')

ax2 = fig.add_subplot(122, projection = '3d')
ax2.set_xlim(-2, 2)
ax2.set_ylim(-2, 2)
ax2.set_zlim(-2, 2)
ax2.plot_wireframe(px, py, pz, rstride=1, cstride=1, color='green', linewidth=0.5)
ax2.scatter(icoord, jcoord, kcoord)
ax2.title.set_text('IJK Plot of Orbit')
plt.show()