# PWBA cross section calculation
# by Henry Schumacher

lazy import argparse
lazy import numpy as np
lazy import matplotlib.pyplot as plt
lazy import xraydb as xy
lazy import time
lazy import tqdm
lazy from functools import wraps
lazy from scipy import integrate
lazy from scipy.special import expn
# lazy import medeleev as md #TODO to be added in future, when medeleev is functional under 3.15

lazy from M_coeff_matrix import load_coefficients

#-----------------------------------------------------------------
def timeit(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        start = time.perf_counter()
        result = func(*args, **kwargs)
        elapsed = time.perf_counter() - start

        print(f"{func.__name__} took {elapsed:.6f} seconds")
        return result

    return wrapper


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
                         'L1':4.15,
                         'L2':4.15,
                         'L3':4.15,
                         'M1':11.25,
                         'M2':11.25,
                         'M3':11.25,
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
    # M = reduced_mass(M1, M2)
    return (ME_EV * E1)/(M1 * Z2s**2 * RYDBERG_EV)

def theta_s(n:int, Z2:int, edge:str):
    E_s = xy.xray_edge(Z2, edge, energy_only=True)
    Z2s = slater_screening(Z2, edge)
    return (n**2 * E_s)/(Z2s**2 * RYDBERG_EV)

def sigma_0s(j2:float, Z1:int, Z2:int, edge:str):
    Z2s = slater_screening(Z2, edge)
    
    
    prefactor   = 2*j2 + 1
    numerator   = 4 * np.pi * Z1**2 * A0_SQ_BARN
    denominator = Z2s**4
    
    return prefactor * numerator / denominator
    
#===============FORM FACTORS===============#    

#===============K-Shell===============#

def form_factor_F_WK(W:float, Q:float, n:int):
    '''
    |F_WK(Q)|^2
    '''
    
    if W > ((1/(n**2)) + 1e-12):
        k = np.sqrt(W - (1/(n**2)))
        
        prefactor = 128 / (1 - np.exp(-2 * np.pi / k))
        first_term = (Q * (Q + (k**2 + 1) / 3)) / ((Q - k**2 + 1)**2 + 4 * k**2)**3
        second_term = np.exp(-(2 / k) * np.arctan2(2 * k, Q - k**2 + 1))
        
        return prefactor * first_term * second_term
    else:
        t = np.sqrt(max((1/(n**2)) - W, 0.0))
        first_term = 128 * (Q * (Q - (t**2 - 1) / 3)) / ((Q + t**2 + 1)**2 - 4 * t**2)**3
        # second_term = ((Q + (t - 1)**2)/(Q + (t + 1)**2))**(1/t)
        if t < 1e-7:
            exponent = -4.0 / (Q + 1.0)
        else:
            exponent = np.log1p(-4.0 * t / (Q + (t + 1.0) ** 2)) / t
        
        return first_term * np.exp(exponent)

#===============L-Shell===============#

def form_factor_F_WL1(W:float, Q:float, n:float):
    '''
    |F_W,2s(Q)|^2 ------- Choi 1973
    '''
    
    if W > ((1/(n**2)) + 1e-12):
        k = np.sqrt(W - 1.0/(n**2))
        
        prefactor = 16 / (1 - np.exp(-2 * np.pi / k))
        first_term = Q / ((Q - k**2 + 1/4)**2 + k**2)**5
        second_term = np.exp(-(2 / k) * np.arctan2(k, Q - k**2 + 1/4))
        third_term = Q**5 - (8/3 + 11/3*k**2)*Q**4 + (41/24 + 6*k**2 + 14/3*k**4)*Q**3
        fourth_term = (5/48 - 31/24*k**2 - 10/3*k**4 - 2*k**6)*Q**2
        fifth_term = (47/3840 - 41/120*k**4 - 2/3*k**6 - 1/3*k**8)*Q
        sixth_term = (1/768 + 17/768*k**2 + 7/48*k**4 + 11/24*k**6 + 2/3*k**8 + 1/3*k**10)
        
        
        return prefactor * first_term * second_term * (third_term + fourth_term + fifth_term + sixth_term)
    else:
        t = np.sqrt(max(1/(n**2) - W, 0.0))
        first_term = 16 * Q / ((Q + t**2 + 1/4)**2 - t**2)**5
        # second_term = ((Q + (t - 1)**2)/(Q + (t + 1)**2))**(1/t)
        third_term = Q**5 - (8/3 - 11/3*t**2)*Q**4 + (41/24 - 6*t**2 + 14/3*t**4)*Q**3
        forth_term = (5/48 + 31/24*t**2 - 10/3*t**4 + 2*t**6)*Q**2
        fifth_term = (47/3840 - 41/120*t**4 + 2/3*t**6 - 1/3*t**8)*Q
        sixth_term = ((1/768) - 17/768*t**2 + 7/48*t**4 - 11/24*t**6 + 2/3*t**8 - 1/3*t**10)
        
        if t < 1e-7:
            exponent = -2.0 / (Q + 1/4)
        else:
            exponent = np.log1p(-2.0 * t / (Q + (t + 1/2) ** 2)) / t
        
        return first_term * np.exp(exponent) * (third_term + forth_term + fifth_term + sixth_term)

def form_factor_F_WL2(W:float, Q:float, n:float):
    '''
    1/3 x |F_W,2p(Q)|^2 ------- Choi 1973
    '''
    
    if W > ((1/(n**2)) + 1e-12):
        k = np.sqrt(W - 1.0/(n**2))
        
        prefactor = 16 / (1 - np.exp(-2 * np.pi / k))
        first_term = Q / ((Q - k**2 + 1/4)**2 + k**2)**5
        second_term = np.exp(-(2 / k) * np.arctan2(k, Q - k**2 + 1/4))
        third_term = 9/4*Q**4 - (3/4 + k**2)*Q**3 + (19/32 - 3/4*k**2 - 1/2*k**4)*Q**2
        forth_term = (107/960 + 41/48*k**2 + 113/60*k**4 + k**6)*Q
        fifth_term = (11/3072 + 3/64*k**2 + 7/32*k**4 + 5/12*k**6 + 1/4*k**8)
        
        return 1/3 * prefactor * first_term * second_term * (third_term + forth_term + fifth_term)
    else:
        t = np.sqrt(max(1/(n**2) - W, 0.0))
        first_term = 16 * Q / ((Q + t**2 + 1/4)**2 - t**2)**5
        third_term = 9/4*Q**4 - (3/4 - t**2)*Q**3 + (19/32 + 3/4*t**2 - 1/2*t**4)*Q**2
        forth_term = (107/960 - 41/48*t**2 + 113/60*t**4 - t**6)*Q
        fifth_term = (11/3072 - 3/64*t**2 + 7/32*t**4 - 5/12*t**6 + 1/4*t**8)
        
        if t < 1e-7:
            exponent = -2.0 / (Q + 1/4)
        else:
            exponent = np.log1p(-2.0 * t / (Q + (t + 1/2) ** 2)) / t
        
        return 1/3 * first_term * np.exp(exponent) * (third_term + forth_term + fifth_term)

def form_factor_F_WL3(W:float, Q:float, n:float):
    '''
    2/3 x |F_W,2p(Q)|^2 ------- Choi 1973
    '''
    
    if W > ((1/(n**2)) + 1e-12):
        k = np.sqrt(W - 1.0/(n**2))
        
        prefactor = 16 / (1 - np.exp(-2 * np.pi / k))
        first_term = Q / ((Q - k**2 + 1/4)**2 + k**2)**5
        second_term = np.exp(-(2 / k) * np.arctan2(k, Q - k**2 + 1/4))
        third_term = 9/4*Q**4 - (3/4 + k**2)*Q**3 + (19/32 - 3/4*k**2 - 1/2*k**4)*Q**2
        forth_term = (107/960 + 41/48*k**2 + 113/60*k**4 + k**6)*Q
        fifth_term = (11/3072 + 3/64*k**2 + 7/32*k**4 + 5/12*k**6 + 1/4*k**8)
        
        return 2/3 * prefactor * first_term * second_term * (third_term + forth_term + fifth_term)
    else:
        t = np.sqrt(max(1/(n**2) - W, 0.0))
        first_term = 16 * Q / ((Q + t**2 + 1/4)**2 - t**2)**5
        third_term = 9/4*Q**4 - (3/4 - t**2)*Q**3 + (19/32 + 3/4*t**2 - 1/2*t**4)*Q**2
        forth_term = (107/960 - 41/48*t**2 + 113/60*t**4 - t**6)*Q
        fifth_term = (11/3072 - 3/64*t**2 + 7/32*t**4 - 5/12*t**6 + 1/4*t**8)
        
        if t < 1e-7:
            exponent = -2.0 / (Q + 1/4)
        else:
            exponent = np.log1p(-2.0 * t / (Q + (t + 1/2) ** 2)) / t
        
        return 2/3 * first_term * np.exp(exponent) * (third_term + forth_term + fifth_term)

#===============M-Shell===============#

def form_factor_F_WMs(W:float, Q:float, n:float, **kwargs):
    '''
    |F_W,M'_3s(Q)|^2 ------- Choi 1973a
    '''
    
    C_3s, pow_3s = kwargs['C_3s'], kwargs['pow_3s']
    
    if W > ((1/(n**2)) + 1e-12):
        k = np.sqrt(W - 1.0/(n**2))
        
        prefactor = 128*Q / (27 * (1 - np.exp(-2 * np.pi / k)))
        second_term = np.exp(-(2 / k) * np.arctan2(2/3 * k, Q - k**2 + 1/9))
        third_term = ((Q - k**2 + 1/9)**2 + 4/9*k**2)**7
        # coefficient_term = 0
        # for i in range(0,10):
        #     for j in range(0,10):
        #         coefficient_term += C_3s[j,i] * k**(2*i) * Q**j
        C_3s_k = C_3s @ ((k**2) ** pow_3s)
        coefficient_term = np.polyval(C_3s_k[::-1], Q)
        return (prefactor * second_term / third_term) * coefficient_term
    else:
        t = np.sqrt(max(1/(n**2) - W, 0.0))
        third_term = 128*Q / (27*((Q + t**2 + 1/9)**2 - 4/9*t**2)**7)
        # coefficient_term = 0
        # for i in range(0,10):
        #     for j in range(0,10):
        #         coefficient_term += C_3s[j,i] * ((-1)**i*t**(2*i)) * Q**j
        C_3s_t = C_3s @ ((-t**2) ** pow_3s)
        coefficient_term = np.polyval(C_3s_t[::-1], Q)
            
        if t < 1e-7:
            exponent = -4/3 / (Q + 1/(n**2)) # TODO: check against Cohen_1985
        else:
            exponent = np.log1p(-4/3 * t / (Q + (t + 1/3) ** 2)) / t
        
        return third_term * np.exp(exponent) * coefficient_term

def form_factor_F_WM2(W:float, Q:float, n:float, **kwargs):
    '''
    |F_W,M'_3p(Q)|^2 ------- Choi 1973a
    '''
    
    C_3p, pow_3p = kwargs['C_3p'], kwargs['pow_3p']
    
    if W > ((1/(n**2)) + 1e-12):
        k = np.sqrt(W - 1.0/(n**2))
        
        prefactor = 128*Q / (27 * (1 - np.exp(-2 * np.pi / k)))
        second_term = np.exp(-(2 / k) * np.arctan2(2/3 * k, Q - k**2 + 1/9))
        third_term = ((Q - k**2 + 1/9)**2 + 4/9*k**2)**7
        # coefficient_term = 0
        # for i in range(0,10):
        #     for j in range(0,10):
        #         coefficient_term += C_3s[j,i] * k**(2*i) * Q**j
        C_3p_k = C_3p @ ((k**2) ** pow_3p)
        coefficient_term = np.polyval(C_3p_k[::-1], Q)
        return (prefactor * second_term / third_term) * coefficient_term
    else:
        t = np.sqrt(max(1/(n**2) - W, 0.0))
        third_term = 128*Q / (27*((Q + t**2 + 1/9)**2 - 4/9*t**2)**7)
        # coefficient_term = 0
        # for i in range(0,10):
        #     for j in range(0,10):
        #         coefficient_term += C_3s[j,i] * ((-1)**i*t**(2*i)) * Q**j
        C_3p_t = C_3p @ ((-t**2) ** pow_3p)
        coefficient_term = np.polyval(C_3p_t[::-1], Q)
            
        if t < 1e-7:
            exponent = -4/3 / (Q + 1/(n**2)) # TODO: check against Cohen_1985
        else:
            exponent = np.log1p(-4/3 * t / (Q + (t + 1/3) ** 2)) / t
        
        return 1/3* third_term * np.exp(exponent) * coefficient_term

def form_factor_F_WM3(W:float, Q:float, n:float, **kwargs):
    '''
    |F_W,M'_3p(Q)|^2 ------- Choi 1973a
    '''
    
    C_3p, pow_3p = kwargs['C_3p'], kwargs['pow_3p']
    
    if W > ((1/(n**2)) + 1e-12):
        k = np.sqrt(W - 1.0/(n**2))
        
        prefactor = 128*Q / (27 * (1 - np.exp(-2 * np.pi / k)))
        second_term = np.exp(-(2 / k) * np.arctan2(2/3 * k, Q - k**2 + 1/9))
        third_term = ((Q - k**2 + 1/9)**2 + 4/9*k**2)**7
        # coefficient_term = 0
        # for i in range(0,10):
        #     for j in range(0,10):
        #         coefficient_term += C_3s[j,i] * k**(2*i) * Q**j
        C_3p_k = C_3p @ ((k**2) ** pow_3p)
        coefficient_term = np.polyval(C_3p_k[::-1], Q)
        return (prefactor * second_term / third_term) * coefficient_term
    else:
        t = np.sqrt(max(1/(n**2) - W, 0.0))
        third_term = 128*Q / (27*((Q + t**2 + 1/9)**2 - 4/9*t**2)**7)
        # coefficient_term = 0
        # for i in range(0,10):
        #     for j in range(0,10):
        #         coefficient_term += C_3s[j,i] * ((-1)**i*t**(2*i)) * Q**j
        C_3p_t = C_3p @ ((-t**2) ** pow_3p)
        coefficient_term = np.polyval(C_3p_t[::-1], Q)
            
        if t < 1e-7:
            exponent = -4/3 / (Q + 1/(n**2)) # TODO: check against Cohen_1985
        else:
            exponent = np.log1p(-4/3 * t / (Q + (t + 1/3) ** 2)) / t
        
        return 2/3 * third_term * np.exp(exponent) * coefficient_term
    
#========^=======FORM FACTORS=======^========#

#========v=======ECPSSR CORRECTIONS======v=========#
def xi_s(eta, theta, n): #
    return 2 * n * np.sqrt(eta) / theta

def g_correction_function(xi, edge):
    denominator = (1 + xi)
    if edge == 'K':
        return (1 + 9*xi + 31*xi**2 + 98*xi**3 + 12*xi**4 + 25*xi**5 + 4.2*xi**6 + 0.515*xi**7)/(denominator**9)
    elif edge == 'L1':
        return (1 + 9*xi + 31*xi**2 + 49*xi**3 + 162*xi**4 + 63*xi**5 + 18*xi**6 + 1.97*xi**7)/(denominator**9)
    elif edge == 'L2':
        return (1 + 10*xi + 45*xi**2 + 102*xi**3 + 331*xi**4 + 6.7*xi**5 + 58*xi**6 + 7.8*xi**7 + 0.888*xi**8)/(denominator**10)

def I_correction_function(x):
    if x <= 0.0:
        return float("inf")
    elif 0 < x <= 0.035:
        return 0.75*np.pi*(np.log(x**(-2)) - 1)
    elif 0.035 < x <= 3.1:
        return np.exp(-2*x)*(0.031 + 0.210*np.sqrt(x) + 0.005*x - 0.069*x**(3/2) + 0.324*x**2)**(-1)
    elif 3.1 < x <= 11.0:
        return np.exp(-2*x)/(x**1.6)
    else:
        return 0.0
    
def h_correction_function(xi, theta, n, edge):
    C_S = {'K':1.5, 'L1':1.5, 'L2':1.25, 'L3':1.25}
    return (2*n/(theta*xi**3) * I_correction_function(C_S[edge] * n / xi))

def zeta_s(theta, xi, n, Z1, Z2, edge):
    '''
    BINDING-POLARIZATION FACTOR zeta(xi)
    '''
    Z2s = slater_screening(Z2, edge)
    factor1 = (2*Z1/(Z2s * theta))
    g_value = g_correction_function(xi, edge)
    h_value = h_correction_function(xi, theta, n, edge)
    return 1 + factor1*(g_value - h_value)

def Coulomb_defl(arg0, edge):
    
    l2 = {'K': 0, 'L1': 0, 'L2': 1, 'L3': 1}
    nu = 10 + 2*l2[edge]
    return (9 + 2*l2[edge]) * expn(nu, arg0)

def rel_correction(xi, n, Z2, edge):
    Z2s = slater_screening(Z2, edge)
    
    if edge in ['K', 'L1']:
        ys = 0.4*(Z2s/137)**2 / (n*xi)
    if edge in ['L2', 'L3']:
        ys = 0.15*(Z2s/137)**2 / xi
        
    return np.sqrt(1 + 1.1*ys**2) + ys

def f_corr_energy_loss(zeta, M1, E1, M, E_edge, edge):
    
    l2 = {'K': 0, 'L1': 0, 'L2': 1, 'L3': 1}
    z = np.sqrt(1.0 - zeta * E_edge * M1 / (M * E1))
    nu = 9 + 2*l2[edge]
    
    return (1/(2**nu * (nu-1)) * ((nu*z -1)*(1+z)**nu + (nu*z + 1)*(1-z)**nu),z)
    

#========^=======ECPSSR CORRECTIONS======^=========#

#========v=======INTEGRATION (Q,W)=======v========#
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
        
def dQ_integration_form_factor(W:float, eta:float, M:float, n:int, ff, GL_points, GL_weights, **kwargs):
    
    Q_lim = np.array(Q_limits(W, eta, M))
    log_lim = np.log(Q_lim)
    if log_lim[1] <= log_lim[0]:
        return 0.0
    Q = np.exp(0.5*(log_lim[1] - log_lim[0])*GL_points + 0.5*(log_lim[1] + log_lim[0]))
    
    return 0.5*(log_lim[1] - log_lim[0]) * np.sum(GL_weights * ff(W, Q, n, **kwargs) / Q)

def dW_integration_form_factor(eta:float, theta:float, M:float, n:int, ff, GL_points, GL_weights, **kwargs):
    
    Wmin = theta / n**2
    Wmax = M * eta / ME_EV
    
    total = 0.0
    if Wmax <= Wmin:
        return 0.0  
    # Bereich W > 1 :  W = 1 + k^2, dW = 2 k dk
    if Wmax > (1/(n**2)):
        kmax = np.sqrt(Wmax - (1/(n**2)))
        steps = [0.0]
        step_value = 0.25
        while step_value < kmax:
            steps.append(step_value)
            step_value *= 2.0
        steps.append(kmax)
        h = lambda k: 2.0 * k * dQ_integration_form_factor(1.0/(n**2) + k**2, eta, M, n, ff, GL_points, GL_weights, **kwargs)
    
        for lo, hi in zip(steps[:-1], steps[1:]):
            int_value, _ = integrate.quad(h, lo, hi, limit=80, epsabs=0.0, epsrel=1e-4)
            total += int_value
    
    # Bereich W < 1 :  W = 1 - t^2, dW = 2 t dt
    if Wmin < (1/(n**2)):
        w_hi = min((1/(n**2)), Wmax)
        t_lo, t_hi = np.sqrt((1/(n**2)) - w_hi), np.sqrt((1/(n**2)) - Wmin)
        g = lambda t: 2.0 * t * dQ_integration_form_factor((1/(n**2)) - t**2, eta, M, n, ff, GL_points, GL_weights,**kwargs)
        int_value, _ = integrate.quad(g, t_lo, t_hi, limit=80, epsabs=0.0, epsrel=1e-4)
        total += int_value
    
    return total
#========^=======INTEGRATION (Q,W)=======^========#

#========v=======CROSS SECTIONS=======v========#
def sigma_PWBA(eta:float, theta:float, M:float, j2:float, Z1:int, Z2:int, n:int, edge:str):
    
    GL_points, GL_weights = np.polynomial.legendre.leggauss(48)
    FORM_FACTORS = {'K': form_factor_F_WK,
                        'L1': form_factor_F_WL1,
                        'L2': form_factor_F_WL2,
                        'L3': form_factor_F_WL3,
                        'M1': form_factor_F_WMs,
                        'M2': form_factor_F_WM2,
                        'M3': form_factor_F_WM3,}
    
    ff = FORM_FACTORS[edge]
    
    sigma_0 = sigma_0s(j2, Z1, Z2, edge)
    
    if edge == 'M1':
        C_3s = load_coefficients(M_='3s')
        pow_3s = np.arange(10)
        f_PWBA = dW_integration_form_factor(eta, theta, M, n, ff, GL_points, GL_weights, **{'C_3s':C_3s, 'pow_3s':pow_3s})
    elif edge in ['M2', 'M3']:
        C_3p = load_coefficients(M_='3p')
        pow_3p = np.arange(9)
        f_PWBA = dW_integration_form_factor(eta, theta, M, n, ff, GL_points, GL_weights, **{'C_3p':C_3p, 'pow_3p':pow_3p})
    else:
        f_PWBA = dW_integration_form_factor(eta, theta, M, n, ff, GL_points, GL_weights)
    
    return sigma_0 * f_PWBA / eta

def sigma_ECPSSR(j2:float, E1, M1:float, M2:float, Z1:int, Z2:int, n:int, edge:str, verbose=False):
    
    sigma_corr = {}
    
    eta = eta_s(E1, Z2, M1, M2, edge)
    theta = theta_s(n, Z2, edge)
    M = reduced_mass(M1, M2)
    
    sigma_corr['PWBA'] = sigma_PWBA(eta, theta, M, j2, Z1, Z2, n, edge)
    if verbose: print('PWBA CROSS SECTION CALCULATED.')
    
    xi = xi_s(eta, theta, n)
    zeta = zeta_s(theta, xi, n, Z1, Z2, edge)
    theta_PSS = zeta*theta
    sigma_corr['PSS'] = sigma_PWBA(eta, theta_PSS, M, j2, Z1, Z2, n, edge)
    if verbose: print('PWBA-PSS CROSS SECTION CALCULATED.')
    
    E_edge = xy.xray_edge(element=Z2, edge=edge).energy
    d = Z1 * Z2 * E2_EV_ANG / (2.0 * E1)
    v1_over_v0 = np.sqrt(E1*ME_EV / M1 / RYDBERG_EV)
    q0 = E_edge / (E2_EV_ANG * v1_over_v0)
    arg0 = np.pi * d * q0 * zeta
    C_correction_0 = Coulomb_defl(arg0, edge)
    sigma_corr['CPSS'] = sigma_corr['PSS'] * C_correction_0
    if verbose: print('PWBA-CPSS CROSS SECTION CALCULATED.')
    
    m_R = rel_correction(xi, n, Z2, edge)
    eta_mR = m_R * eta
    sigma_corr['CPSSR'] = sigma_PWBA(eta_mR, theta_PSS, M, j2, Z1, Z2, n, edge) * C_correction_0
    if verbose: print('PWBA-CPSSR CROSS SECTION CALCULATED.')
    
    _ , z = f_corr_energy_loss(zeta, M1, E1, M, E_edge, edge)
    arg_el = (2 * arg0) / (z * (1.0 + z))
    C_correction_el = Coulomb_defl(arg_el, edge)
    sigma_corr['ECPSSR'] = sigma_PWBA(eta_mR, theta_PSS, M, j2, Z1, Z2, n, edge) * C_correction_el
    if verbose: print('ECPSSR CROSS SECTION CALCULATED.')

    return sigma_corr

#========^=======CROSS SECTIONS=======^========#

#========v=======PRODUCTION OF SIM. DATA=======v========#
def PWBA_ECPSSR_calculate(E1_MeV, M1_u:float, M2_u:float, Z1:int, Z2:int, edge:str, verbose:bool=False):
    E1 = E1_MeV * 1e6    
    M1 = M1_u * U_IN_EV
    M2 = M2_u * U_IN_EV
    
    EDGE_N = {'K': 1, 'L1': 2, 'L2': 2, 'L3': 2, 'M1': 3, 'M2': 3, 'M3': 3, 'M4': 3, 'M5': 3}
    EDGE_J2 = {'K': 0.5, 'L1': 0.5, 'L2': 0.5, 'L3': 1.5, 'M1': 0.5, 'M2': 0.5, 'M3': 1.5, 'M4': 1.5, 'M5': 2.5}

    sigma_corr = sigma_ECPSSR(EDGE_J2[edge], E1, M1, M2, Z1, Z2, EDGE_N[edge], edge, verbose)
    
    return sigma_corr


if __name__ == "__main__":
    print(PWBA_ECPSSR_calculate(E1_MeV=15, M1_u=xy.atomic_mass('H'), M2_u=xy.atomic_mass('Pu'), Z1=1, Z2=xy.atomic_number('Pu'), edge='L1', verbose=False))