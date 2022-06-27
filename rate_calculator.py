import sys
import Hirogen_Functions

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson, quad
import scipy.constants as constants
from astropy import units as u

very_shortnames = ["'SDSS J0748'", "'SDSS J0938'", "'SDSS J0952'", "'SDSS J1055'", "'SDSS J1241'", "'SDSS J1342'", "'SDSS J1350'"]

very_shortname = very_shortnames[0]


def sigmoid(x, A, K, B, V, Q):
    return A+((K-A)/((1+Q*np.exp(-B*x))**(1/V)))

def power_law(t, M, N):
    return -M*t**(-N)

def curve(time, A, K, B, V, Q, M, N):
    evo = power_law(time, M, N)
    return evo * sigmoid(evo, A, K, B, V, Q)


if very_shortname in "'SDSS J0952', 'SDSS J1241'":
    a, k, b, v, q = np.genfromtxt('FeVII_parameters.csv')
elif very_shortname in "'SDSS J0748', 'SDSS J1342', 'SDSS J1350'":
    a, k, b, v, q = np.genfromtxt('Non_FeVII_parameters.csv')

m, n = np.genfromtxt('power_params.csv')

time = np.arange(0,6000,1)

no_ecles = 5

no_objects = 681730

def curve(time, a, k, b, v, q, m, n):
    evo = power_law(time, m, n)
    return evo * sigmoid(evo, a, k, b, v, q)

visibility_time = -quad(curve, time[0], time[-1], args = (a, k, b, v, q, m, n))[0]

day_rate = no_ecles/(visibility_time*no_objects)
year_rate = day_rate*365