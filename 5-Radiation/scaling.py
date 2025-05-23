import numpy as np

class CGS():
    name = 'CGS'
    m_Earth = 5.972e27      # g
    m_Sun = 1.989e33        # g
    r_Earth = 6.371e8       # cm
    grav = 6.67430e-8         # cm^3 g^-1 s^-2
    G = grav
    yr = 3.156e+7           # s
    au = 1.498e+13          # cm
    AU = au
    mu = 2.4                # molecular weight per atomic mass unit
    m_u = 1.6726e-24        # g
    m_p = m_u
    k_B = 1.3807e-16        # erg per K
    h_P = 6.62606e-27       # erg per Hertz
    e = 4.8032e-10          # stat coulomb
    c = 2.9979e10           # cm/s
    Stefan = 5.670374419e-5 # erg cm^-2 s^-1 K^-4
    pc = 180./np.pi*3600.*au
    kyr = 1e3*yr
    Myr = 1e6*yr
    kms = 1e5
    Angstrom = 1e-8         # Å / cm
    micron = 1e-4           # mu / cm

class MKS():
    name = 'MKS'
    m_Earth = 5.972e24   # kg
    m_Sun = 1.989e30     # kg
    r_Earth = 6.371e6    # m
    grav = 6.67430e-11   # m^3 kg^-1 s^-2
    G = grav
    yr = 3.156e+7        # s
    au = 1.498e+11       # m
    AU = au
    mu = 2.4             # molecular weight per atomic mass unit
    m_u = 1.6726e-27     # kg
    m_p = m_u
    k_B = 1.3807e-23     # J per K
    h_P = 6.62606e-34    # J per Hertz
    e = 1.6022e-19       # coulomb
    c = 2.9979e8         # m/s
    Stefan = 5.670374419e-8 # J m^-2 s^-1 K^-4
    pc = 180./np.pi*3600.*au
    kyr = 1e3*yr
    Myr = 1e6*yr
    kms = 1e3
    Angstrom = 1e-10     # Å / m
    micron = 1e-6        # mu / m

class SI(MKS):
    name = 'SI'

class scaling():
    '''
    Return a structure with scaling constants.  Use e.g.

        scgs=scaling(cgs)
        sSI=scaling(SI)
        print (cgs.k_b, SI.k_b)
    '''
    def __init__(self,system=CGS,l=0,m=0,t=0,D=0,v=0,mu=2.4,verbose=0):
        """ Object holding scaling information for code units """

        if type(system)==type('str'):
            if   system=='cgs' or system=='CGS':
                system = CGS
            elif system=='mks' or system=='MKS':
                system = MKS
            elif system=='si'  or system=='SI':
                system = SI
        self.system = system

        if verbose>0:
            print("using "+system.name+" units")

        # Interpret string values
        if l == 'pc'    : l = system.pc
        if m == 'Solar' : m = system.m_Sun
        if m == 'm_Sun' : m = system.m_Sun
        if t == 'yr'    : t = system.yr
        if t == 'kyr'   : t = system.kyr
        if t == 'Myr'   : t = system.Myr
        if v == 'kms'   : v = system.kms

        # Count input scaling constants
        inputs = [l,m,t,D,v]
        nonzero = [i>0 for i in inputs].count(True)
        if verbose > 2:
            print('there','is' if nonzero==1 else 'are',nonzero,'input value'+('s' if nonzero>1 else ''))
        if nonzero != 3:
            print('exactly 3 input values must be specified (there {} {})'.\
                  format('was' if nonzero==1 else 'were',nonzero))
            return None
        
        self.l = l
        self.m = m
        self.t = t
        self.D = D
        self.v = v

        # If velocity units were given, either length or time units were not
        if v > 0:
            if l == 0:
                self.l = v*t
            else:
                self.t = l/v
        else:
            self.v = l/t

        # If density units were given, either length or mass unis were not
        if D > 0:
            if m == 0:
                self.m = D*l**3
            else:
                self.l = (m/D)**(1/3)
        else:
            self.D = m/l**3

        # Compute derived scaling constants for other quantities
        self.m = self.D*self.l**3                           # mass
        self.P = self.D*self.v**2                           # pressure
        self.e = self.m*self.v**2                           # energy
        self.E = self.D*self.v**2                           # energy density
        self.mu = mu                                        # mean molecular weight
        self.T = mu*(system.m_u)/(system.k_B) * self.v**2   # temperature given as P/rho

        # Scaling of constants of nature into code units
        self.G = system.G*self.D*self.t**2                  # constant of gravity
        self.Stefan = system.Stefan/(self.E*self.v)         # Stefan's constant, rescaled
        self.h_P = system.h_P/(self.t*self.e)               # Planck's constant, rescaled
        self.k_B = system.k_B/self.e                        # Boltzmann's constant, rescaled
        self.c = system.c/self.v                            # speed of light, rescaled
        self.m_u = system.m_u/self.m                        # atomic mass unit, rescaled
        if verbose>1:
            print(vars(self))