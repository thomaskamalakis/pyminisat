import numpy as np
from scipy.special import erfc, erfinv
"""
Some basic utility functions to be used in the simulations
"""

# Q function
def Q(x):
    return 0.5 * erfc(x / np.sqrt(2))

# Inverse Q function
def invQ(y):
    return np.sqrt(2)*erfinv(1-2*y)

# From dB to linear scale
def to_linear(x):
    return 10 ** (x/10)

# From linear to dB scale
def to_dB(x):
    return 10 * np.log10(x)
    
# From mW to dBm 
def to_dBm(x):
    return 10 * np.log10(x/1e-3)