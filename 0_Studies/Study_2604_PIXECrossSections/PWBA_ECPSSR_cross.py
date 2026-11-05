# PWBA cross section calculation
# by Henry Schumacher

lazy import argparse
lazy import numpy as np
lazy import matplotlib.pyplot as plt
lazy import xraydb as xy

lazy from scipy import integrate
lazy from scipy.special import expn
# lazy import medeleev as md #TODO to be added in future, when medeleev is functional under 3.15

# ---------------------------------------------------------------- Konstanten
RYDBERG_EV  = 13.605693          # eV
E2_EV_ANG   = 14.399645          # e^2 in eV*Angstrom
A0_SQ_BARN  = 2.800285e7         # a0^2 in barn (a0^2 = 2.800285e-17 cm^2)
U_OVER_ME   = 1822.888486        # atomare Masseneinheit / Elektronenmasse
U_IN_EV     = 931.4941e6
ME_EV       = 0.51099e6          # eV

# print(xy.xray_edge('Ti','K',energy_only=True))

def slater_screening(Z2:int, edge:str):
    
    slater_screening_values = {'K':0.3,
                         'L':4.15,
                         'M1':11.25,
                         'M2':11.25,
                         'M3':21.15,
                         'M4':21.15,
                         'M5':21.15}
    
    if edge not in slater_screening_values.keys():
        return 'ERROR: Incorrect edge given.'
    else:
        return float(Z2 - slater_screening_values[edge])

def reduced_mass(M1:float, M2:float):
    return (M1 * M2)/(M1 + M2)

def eta_s(E1:float,Z2:int, M1:float, M2:float, edge:str):
    Z2s = slater_screening(Z2, edge)
    M = reduced_mass(M1, M2)
    return (ME_EV * E1)/(M1 * Z2s**2 * RYDBERG_EV) #TODO This equation should not have M1 but M. Check ECPSSR results once done.

def theta_s(n:int, Z2:int, edge:str):
    E_s = xy.xray_edge(Z2, edge, energy_only=True)          #eV
    Z2s = slater_screening(Z2, edge)
    return (n**2 * E_s)/(Z2s**2 * RYDBERG_EV)

def sigma_0s(j2:float, Z1:int, Z2:int, edge:str):
    Z2s = slater_screening(Z2, edge)
    
    
    prefactor   = 2*j2 + 1
    numerator   = 4 * np.pi * Z1**2 * A0_SQ_BARN
    denominator = Z2s**4
    
    return prefactor * numerator / denominator
    
def form_factor_F_WK(W:float, Q:float, n:int):
    '''
    |F_WK(Q)|^2
    '''
    
    if W > (1.0 + 1e-12):
        k = np.sqrt(W - 1.0/n)
        
        prefactor = 128 / (1 - np.exp(-2 * np.pi / k))
        first_term = (Q * (Q + (k**2 + 1) / 3)) / ((Q - k**2 + 1)**2 + 4 * k**2)**3
        second_term = np.exp(-(2 / k) * np.arctan2(2 * k, Q - k**2 + 1))
        
        return prefactor * first_term * second_term
    else:
        t = np.sqrt(max(1.0 - W, 0.0))
        first_term = 128 * (Q * (Q - (t**2 - 1) / 3)) / ((Q + t**2 + 1)**2 - 4 * t**2)**3
        # second_term = ((Q + (t - 1)**2)/(Q + (t + 1)**2))**(1/t)
        if t < 1e-7:
            exponent = -4.0 / (Q + 1.0)
        else:
            exponent = np.log1p(-4.0 * t / (Q + (t + 1.0) ** 2)) / t
        
        return first_term * np.exp(exponent)
    
def W_limits(eta:float, theta:float, M:float, n:int):
    return [theta/(n**2), M * eta / ME_EV] # Wmin, Wmax

def Q_limits(W:float, eta:float, M:float):
    prefactor = (M / ME_EV)**2 * eta
    
    if ME_EV * W / (M * eta) >= 1.0:
        root_term = 0.0
    else:
        root_term = np.sqrt(1 - ME_EV * W / (M * eta))
    
    Qmin = prefactor * (1 - root_term)**2
    Qmax = prefactor * (1 + root_term)**2
    return [Qmin, Qmax]
    
def dQ_integration_form_factor(W:float, eta:float, M:float, n:int, edge:str):
      
    Q_lim = np.array(Q_limits(W, eta, M))
    
    if edge == 'K':
        f = lambda x: form_factor_F_WK(W, np.exp(x), n) / np.exp(x)
        log_lim = np.log(Q_lim)
        pts = [p for p in (0.0,) if log_lim[0] < p < log_lim[1]]
        result, _ = integrate.quad(f, log_lim[0], log_lim[1], points=pts or None, limit=300, epsabs=0.0, epsrel=1e-8)
        return result
    
def dW_integration_form_factor(eta:float, theta:float, M:float, n:int, edge:str):
    
    Wmin = theta / n**2
    Wmax = M * eta / ME_EV
    
    total = 0.0
    if Wmax <= Wmin:
        return 0.0  
    # Bereich W > 1 :  W = 1 + k^2, dW = 2 k dk
    if Wmax > 1.0:
        kmax = np.sqrt(Wmax - 1.0)
        edges = [0.0]
        e = 0.25
        while e < kmax:
            edges.append(e)
            e *= 2.0
        edges.append(kmax)
        h = lambda k: 2.0 * k * dQ_integration_form_factor(1.0 + k**2, eta, M, n, edge)
        for lo, hi in zip(edges[:-1], edges[1:]):
            int_value, _ = integrate.quad(h, lo, hi, limit=400, epsabs=0.0, epsrel=1e-7)
            total += int_value
    
    # Bereich W < 1 :  W = 1 - t^2, dW = 2 t dt
    if Wmin < 1.0:
        w_hi = min(1.0, Wmax)
        t_lo, t_hi = np.sqrt(1.0 - w_hi), np.sqrt(1.0 - Wmin)
        g = lambda t: 2.0 * t * dQ_integration_form_factor(1.0 - t**2, eta, M, n, edge)
        int_value, _ = integrate.quad(g, t_lo, t_hi, limit=400, epsabs=0.0, epsrel=1e-7)
        total += int_value
    
    return total

def sigma_PWBA(eta:float, theta:float, M:float, j2:float, Z1:int, Z2:int, n:int, edge:str):
    
    sigma_0 = sigma_0s(j2, Z1, Z2, edge)
    f_PWBA = dW_integration_form_factor(eta, theta, M, n, edge)
    
    return sigma_0 * f_PWBA / eta



def test_PWBA(E1_MeV:float, M1_u:float, M2_u:float, Z1:int, Z2:int, n:int, edge:str):
    
    E1 = E1_MeV * 1e6
    M1 = M1_u * U_IN_EV
    M2 = M2_u * U_IN_EV
    
    j2 = 0.5
    
    eta = eta_s(E1, Z2, M1, M2, edge)
    theta = theta_s(n, Z2, edge)
    M = reduced_mass(M1, M2)
    
    return sigma_PWBA(eta, theta, M, j2, Z1, Z2, n, edge)

'''
sigma_ECPSSR = correction_term * sigma_PWBA(msR*eta, zeta_s*theta,...)
'''


if __name__ == "__main__":
    result = test_PWBA(E1_MeV=1.0, M1_u=1.007276, M2_u=63.546, Z1=1, Z2=29, n=1, edge='K')
    print(f'{result:12.4f} b')
    print(result)