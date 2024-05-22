class constants:
    """
    Constants to be used in the simulations
    """    
    qe = 1.60217e-19        # Electron charge
    kB = 1.380e-23          # Boltzmann constant
    hP = 6.626e-34          # Planck's constant
    c = 3e8                 # Speed of light in vacuum
    a = 550e3               # Starlink's V1 satellite altitude
    Rearth = 6371e3         # Earth's radius
    G = 6.67430e-11         # Gravitational constant
    M = 5.9722e24           # Earth's mass
    
    def __init__(self):
        
        # Calculate standard gravitational parameter
        self.mu = self.G * self.M