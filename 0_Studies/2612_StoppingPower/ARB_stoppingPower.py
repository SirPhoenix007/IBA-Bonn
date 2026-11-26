# Ashley, Ritchie & Brandt Stopping Power
#by H. Schumacher

lazy import argparse
lazy import numpy as np
lazy import matplotlib.pyplot as plt
lazy import time
lazy import tqdm
lazy import math
lazy from functools import wraps
lazy from scipy import integrate
lazy import scipy.special as scsp

'''
SOURCES:
Ashley_1972.pdf -> (1)
Fernandez_Varea_2026.pdf -> (2)
Sannikov_1958.pdf -> (3)
Jensen_1932.pdf -> (4)
Lenz_1932.pdf -> (5)
'''

#F-Function (1.14) / (2.13/14)
#Xi-Function
#x & w
#eta, chi
#L(x) function (1.13)

V_ZERO = 2.19*1e6
BOHR_RADIUS = 0.0529e-9
MYSTIC_CONSTANT = 0.8853

def x_parameter(Z1:int, Z2:int, E1:float):
    '''
    Index 1: Projectile\n
    Z1, E1 in MeV\n
    \n
    Index 2: Target\n
    Z2\n
    \n
    Conversion from (1)\n
    v_1^2/v_0^2 = 40.2*E1 in MeV/amu
    '''
    x = float((40.2*E1/Z1)/Z2)
    v1 = float(np.sqrt(40.2*E1/Z1 * V_ZERO**2))
    
    return x, v1

#########
# Rho Lenz-Jensen
#########
A_0 = 0.52918 #Angström per a_0
kappa = 0.6 * (3*np.pi / (8*np.sqrt(2))) ** (2/3) #Jensen (1a)

def P0(c1):
    return 2.0 * sum(math.comb(3,k) * math.factorial(k + 2) * c1**k for k in range(4))

def lambda_LJ(Z,T):
    return 1.0 / (Z**(1/3) * T)

def Jensen_Energy_Integrals(c):
    P_0 = P0(c)
    f = lambda s: np.exp(-s)*(1+c*s)**3 / (4*np.pi*P_0*s**3) # density function
    P_1, _ = integrate.quad(lambda s: 4*np.pi*s**2 * f(s) * 2*s, 0, 80) #(7a)
    
    Q = lambda s: integrate.quad(lambda t: 4*np.pi*t**4 * f(t) * 2*t, 0, s)[0]
    out = lambda s: integrate.quad(lambda t: 4*np.pi*t**2 * f(t) * 2*t, s, 80)[0]
    
    P_2 = 0.5*integrate.quad(lambda s: 4*np.pi*s**4 * f(s) * (Q(s)/s**2 + out(s))*2*s, 1e-9, 60, limit=200)[0] #(7b)
    P_3 = kappa*(4*np.pi)**(2/3) * integrate.quad(lambda s: f(s)**(5/3)*4*np.pi*s**4*2*s, 1e-9, 80, limit=200)[0] #(7c)
    
    return P_1, P_2, P_3

def F12(c, Z, N):
    P_1, P_2, P_3 = Jensen_Energy_Integrals(c)
    return (N/Z)*P_1 - (N/Z)**2*P_2, (N/Z)**(5/3)*P_3 #(7d)
   
def lambda_calculation(Z, N):
    cgrid = np.linspace(0.10, 0.40, 151)
    best_c1 = max(cgrid, key=lambda c: (lambda F: F[0]**2 / F[1])(F12(c, Z, N)))
    F1, F2 = F12(best_c1, Z, N)
    return best_c1, 1.0/(Z**(1/3) * 0.5 * F1/F2) #c1, lambda in a_0
    
def rho_LJ(r, N, c1, lambda_value):
    x = np.sqrt(r/lambda_value)
    return N/(4 * np.pi * lambda_value**3 * P0(c1)) * np.exp(-x) * (1 + c1*x)**3 / x**3

# def rho_Lenz_Jensen(r, mu, Z1):
#     '''
#     (3)
#     '''
#     xi = r / mu
#     lambda_LJ = BOHR_RADIUS/(Z1**(1/3)*3.458) #(4)
#     x = lambda_LJ * xi
#     C = 0.265 #(4)
    
#     return Z1/(4*np.pi*mu**3) * lambda_LJ**3 * x**(-3/2) * np.exp(-x**(1/2)) * (1 + C*x**0.5)**3

# def plot_Lenz_Jensen(Z1):
#     mu = MYSTIC_CONSTANT * BOHR_RADIUS / Z1**(1/3)
#     plt.figure(figsize=(4,4), dpi=200)
#     r = np.linspace(0.01,1e2,20000)
#     Rho = lambda r: r**2 * 4 * np.pi * rho_Lenz_Jensen(r, mu, Z1)
#     I_LJ, _ = integrate.quad(Rho, 0.01*BOHR_RADIUS, 2000*BOHR_RADIUS, limit=100, epsabs=0.0, epsrel=1e-5)
#     print(f'P = {I_LJ/(4*Z1)}')
#     # plt.plot(r, rho_Lenz_Jensen(r, mu, Z1))
#     plt.plot(r,r**0.5
#             * np.exp(-np.sqrt(r))
#             * (1 + 0.265 * np.sqrt(r))**3)
#     # plt.xlim(0,1000*BOHR_RADIUS)
#     plt.show()



# import numpy as np
# from scipy.integrate import quad
# from scipy.constants import physical_constants

# a0 = physical_constants["Bohr radius"][0]  # meters

# def rho_LJ(r, Z, P, lam, C=0.265):
#     mu = 0.8853 * a0 * Z**(-1/3)
#     x = lam * r / mu

#     if np.isscalar(r) and r == 0:
#         return np.inf

#     return (
#         Z / (4 * np.pi * mu**3)
#         * lam**3 / (4 * P)
#         * x**(-1.5)
#         * np.exp(-np.sqrt(x))
#         * (1 + C * np.sqrt(x))**3
#     )

# def electron_number_x(Z, P, C=0.333):

#     integrand = lambda x: (
#         x**0.5
#         * np.exp(-np.sqrt(x))
#         * (1 + C * np.sqrt(x))**3
#     )

#     integral, error = quad(
#         integrand, 0, np.inf,
#         epsabs=1e-10,
#         epsrel=1e-12
#     )

#     Ne = Z / (4 * P) * integral

#     return Ne, error


if __name__ == "__main__":
    #Test K+ Jensen Values
    Z = 37
    N = 37
    c1, L = lambda_calculation(Z,N)
    s = np.linspace(1e-6, 40, 400001) # s = x = sqrt(r/lambda)
    r = L * s**2
    total_Q = integrate.quad(lambda r: 4 * np.pi * r**2 * rho_LJ(r, N, c1, L), s[0], s[-1])
    print(total_Q)
    plt.plot(r, 4*np.pi*r**2*rho_LJ(r, N, c1, L)) #x-axis: r in a_0, y-axis in electrons/a_0
    plt.xlim(-0.1,6)
    plt.show()
    
    # Z = 50
    # mu = 0.8853 * a0 * Z**(-1/3)
    # lam = 6.119     #UMEDA / Replace with the corresponding lambda
    # P = 10.2     #JENSEN / Replace with the value from the original parametrization
    
    # Ne, error = electron_number_x(Z, P)

    # print(f"Integrated electron number: {Ne:.10f}")
    # print(f"Expected electron number:   {Z}")
    # print(f"Integration error:          {error:.3e}")
    
    # plot_Lenz_Jensen(Z)