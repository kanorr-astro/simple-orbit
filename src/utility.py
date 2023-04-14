#Utility Functions for simple-orbit
import numpy as np
import math as math


def magnitude(vector):
    return math.sqrt(sum(pow(element, 2) for element in vector))

def spec_mech_energy(r,v,mu):
    SME = (magnitude(v)**2 / 2) - (mu / magnitude(r))
    return SME   

def radius(p, e, nu):
    return p / (1 + e * math.cos(nu))

