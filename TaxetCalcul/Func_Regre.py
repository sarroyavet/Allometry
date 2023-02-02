import numpy as np
from cmath import e

def PopHirt():
    return([25.5, 0.26, 22, -0.6])

# Allometric functions
def regfunct(X,a,b):
    return a * np.power(X, b)

def SpeedHirt(X,a,b,h,i):
    return a * np.power(X, b)*(1-np.power(e,-h*np.power(X,i)))