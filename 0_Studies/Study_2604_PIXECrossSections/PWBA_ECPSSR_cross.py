# PWBA cross section calculation
# by Henry Schumacher

lazy import argparse
lazy import numpy as np
lazy import matplotlib.pyplot as plt
lazy import xraydb as xy

lazy from scipy import integrate
lazy from scipy.special import expn
# lazy import medeleev as md #TODO to be added in future, when medeleev is functional under 3.15

lazy from M_coeff_matrix import load_coefficients


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
        second_term = np.exp(-(2 / k) * np.arctan2(2 * k, Q - k**2 + 1/4))
        third_term = Q**5 - (8/3 + 11/3*k**2)*Q**4 + (41/24 + 6*k**2 + 14/3*k**4)*Q**3
        fourth_term = (5/48 - 31/24*k**2 - 10/3*k**4 - 2*k**6)*Q**2
        fifth_term = (47/3840 - 41/120*k**4 - 2/3*k**6 - 1/3*k**8)*Q
        sixth_term = (1/768 + 17/768*k**2 + 7/48*k**4 + 11/24*k**6 + 2/3*k**8 + 1/3*k**10)
        
        
        return prefactor * first_term * second_term * (third_term + fourth_term + fifth_term + sixth_term)
    else:
        t = np.sqrt(max(1/(n**2) - W, 0.0))
        first_term = 16 * Q / ((Q + t**2 + 1/4)**2 - t**2)**5
        # second_term = ((Q + (t - 1)**2)/(Q + (t + 1)**2))**(1/t)
        third_term = Q**5 - (8/3 + 11/3*t**2)*Q**4 + (41/24 + 6*t**2 + 14/3*t**4)*Q**3
        forth_term = (5/48 + 31/24*t**2 + 10/3*t**4 + 2*t**6)*Q**2
        fifth_term = (47/3840 + 41/120*t**4 + 2/3*t**6 + 1/3*t**8)*Q
        sixth_term = -(-(1/768) + 17/768*t**2 + 7/48*t**4 + 11/24*t**6 + 2/3*t**8 + 1/3*t**10)
        
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
        second_term = np.exp(-(2 / k) * np.arctan2(2 * k, Q - k**2 + 1/4))
        third_term = 9/4*Q**4 - (3/4 + k**2)*Q**3 + (19/32 - 3/4*k**2 - 1/2*k**4)*Q**2
        forth_term = (107/960 + 41/48*k**2 + 113/60*k**4 + k**6)*Q
        fifth_term = (11/3072 + 3/64*k**2 + 7/32*k**4 + 5/12*k**6 + 1/4*k**8)
        
        return 1/3 * prefactor * first_term * second_term * (third_term + forth_term + fifth_term)
    else:
        t = np.sqrt(max(1/(n**2) - W, 0.0))
        first_term = 16 * Q / ((Q + t**2 + 1/4)**2 - t**2)**5
        third_term = 9/4*Q**4 - (3/4 - t**2)*Q**3 + (19/32 + 3/4*t**2 + 1/2*t**4)*Q**2
        forth_term = (107/960 - 41/48*t**2 - 113/60*t**4 - t**6)*Q
        fifth_term = (11/3072 - 3/64*t**2 - 7/32*t**4 - 5/12*t**6 - 1/4*t**8)
        
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
        second_term = np.exp(-(2 / k) * np.arctan2(2 * k, Q - k**2 + 1/4))
        third_term = 9/4*Q**4 - (3/4 + k**2)*Q**3 + (19/32 - 3/4*k**2 - 1/2*k**4)*Q**2
        forth_term = (107/960 + 41/48*k**2 + 113/60*k**4 + k**6)*Q
        fifth_term = (11/3072 + 3/64*k**2 + 7/32*k**4 + 5/12*k**6 + 1/4*k**8)
        
        return 1/3 * prefactor * first_term * second_term * (third_term + forth_term + fifth_term)
    else:
        t = np.sqrt(max(1/(n**2) - W, 0.0))
        first_term = 16 * Q / ((Q + t**2 + 1/4)**2 - t**2)**5
        third_term = 9/4*Q**4 - (3/4 - t**2)*Q**3 + (19/32 + 3/4*t**2 + 1/2*t**4)*Q**2
        forth_term = (107/960 - 41/48*t**2 - 113/60*t**4 - t**6)*Q
        fifth_term = (11/3072 - 3/64*t**2 - 7/32*t**4 - 5/12*t**6 - 1/4*t**8)
        
        if t < 1e-7:
            exponent = -2.0 / (Q + 1/4)
        else:
            exponent = np.log1p(-2.0 * t / (Q + (t + 1/2) ** 2)) / t
        
        return 2/3 * first_term * np.exp(exponent) * (third_term + forth_term + fifth_term)

#===============M-Shell===============#

def form_factor_F_WM1(W:float, Q:float, n:float):
    '''
    |F_W,M'(Q)|^2 ------- Choi 1973a
    '''
    
    C_3s = load_coefficients(M_='3s')
    
    if W > ((1/(n**2)) + 1e-12):
        k = np.sqrt(W - 1.0/(n**2))
        
        prefactor = 128 / (27 * (1 - np.exp(-2 * np.pi / k)))
        second_term = np.exp(-(2 / k) * np.arctan2(2/3 * k, Q - k**2 + 1/9))
        third_term = ((Q - k**2 + 1/9)**2 + 4/9*k**2)**7
        coefficient_term = 0
        for i in range(0,10):
            for j in range(0,10):
                coefficient_term += C_3s[j,i] * k**(2*i) * Q**j
        return (prefactor * second_term / third_term) * coefficient_term
    else:
        t = np.sqrt(max(1/(n**2) - W, 0.0))
        third_term = 128 / (27*((Q + t**2 + 1/9)**2 - 4/9*t**2)**7)
        coefficient_term = 0
        for i in range(0,10):
            for j in range(0,10):
                coefficient_term += C_3s[j,i] * (t**(2*i)) * Q**j
            
        if t < 1e-7:
            exponent = -4/3 / (Q + 1/(n**2)) # TODO: CHEK WIT PAPS
        else:
            exponent = np.log1p(-4/3 * t / (Q + (t + 1/3) ** 2)) / t
        
        return third_term * np.exp(exponent) * coefficient_term
#===============FORM FACTORS===============#

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
        result, _ = integrate.quad(f, log_lim[0], log_lim[1], points=pts or None, limit=200, epsabs=0.0, epsrel=1e-8)
        return result
    elif edge == 'L1':
        f = lambda x: form_factor_F_WL1(W, np.exp(x), n) / np.exp(x)
        log_lim = np.log(Q_lim)
        pts = [p for p in (0.0,) if log_lim[0] < p < log_lim[1]]
        result, _ = integrate.quad(f, log_lim[0], log_lim[1], points=pts or None, limit=200, epsabs=0.0, epsrel=1e-8)
        return result
    elif edge == 'L2':
        f = lambda x: form_factor_F_WL2(W, np.exp(x), n) / np.exp(x)
        log_lim = np.log(Q_lim)
        pts = [p for p in (0.0,) if log_lim[0] < p < log_lim[1]]
        result, _ = integrate.quad(f, log_lim[0], log_lim[1], points=pts or None, limit=200, epsabs=0.0, epsrel=1e-8)
        return result
    elif edge == 'L3':
        f = lambda x: form_factor_F_WL3(W, np.exp(x), n) / np.exp(x)
        log_lim = np.log(Q_lim)
        pts = [p for p in (0.0,) if log_lim[0] < p < log_lim[1]]
        result, _ = integrate.quad(f, log_lim[0], log_lim[1], points=pts or None, limit=200, epsabs=0.0, epsrel=1e-8)
        return result
    elif edge == 'M1':
        f = lambda x: form_factor_F_WM1(W, np.exp(x), n) / np.exp(x)
        log_lim = np.log(Q_lim)
        pts = [p for p in (0.0,) if log_lim[0] < p < log_lim[1]]
        result, _ = integrate.quad(f, log_lim[0], log_lim[1], points=pts or None, limit=200, epsabs=0.0, epsrel=1e-8)
        return result

def dW_integration_form_factor(eta:float, theta:float, M:float, n:int, edge:str):
    
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
        h = lambda k: 2.0 * k * dQ_integration_form_factor(1.0/(n**2) + k**2, eta, M, n, edge)
    
        for lo, hi in zip(steps[:-1], steps[1:]):
            int_value, _ = integrate.quad(h, lo, hi, limit=400, epsabs=0.0, epsrel=1e-7)
            total += int_value
    
    # Bereich W < 1 :  W = 1 - t^2, dW = 2 t dt
    if Wmin < (1/(n**2)):
        w_hi = min((1/(n**2)), Wmax)
        t_lo, t_hi = np.sqrt((1/(n**2)) - w_hi), np.sqrt((1/(n**2)) - Wmin)
        g = lambda t: 2.0 * t * dQ_integration_form_factor((1/(n**2)) - t**2, eta, M, n, edge)
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
    
    #TODO check j2 values
    if edge in ['K', 'L1', 'M1', 'M2']:
        j2 = 0.5
    elif edge in ['L2', 'L3', 'M3', 'M4']:
        j2 = 1.5
    elif edge in ['M5']:
        j2 = 2.5

    eta = eta_s(E1, Z2, M1, M2, edge)
    theta = theta_s(n, Z2, edge)
    M = reduced_mass(M1, M2)
    
    return sigma_PWBA(eta, theta, M, j2, Z1, Z2, n, edge)

'''
sigma_ECPSSR = correction_term * sigma_PWBA(msR*eta, zeta_s*theta,...)
'''


if __name__ == "__main__":
    # print('PROTON ON AMERICIUM')
    # result = test_PWBA(E1_MeV=1.0, M1_u=1.007276, M2_u=241.057, Z1=1, Z2=95, n=1, edge='K')
    # print(f'K-shell: {result:12.3f} b')
    # result = test_PWBA(E1_MeV=1.0, M1_u=1.007276, M2_u=241.057, Z1=1, Z2=95, n=2, edge='L1')
    # print(f'L1-shell: {result:12.3f} b')
    # result = test_PWBA(E1_MeV=1.0, M1_u=1.007276, M2_u=241.057, Z1=1, Z2=95, n=2, edge='L2')
    # print(f'L2-shell: {result:12.3f} b')
    # result = test_PWBA(E1_MeV=1.0, M1_u=1.007276, M2_u=241.057, Z1=1, Z2=95, n=2, edge='L3')
    # print(f'L3-shell: {result:12.3f} b')
    # result = test_PWBA(E1_MeV=1.0, M1_u=1.007276, M2_u=241.057, Z1=1, Z2=95, n=3, edge='M1')
    # print(f'M1-shell: {result:12.3f} b')
    for E in np.linspace(0.2,2,40):
        result = test_PWBA(E1_MeV=E, M1_u=1.007276, M2_u=xy.atomic_mass('K'), Z1=1, Z2=xy.atomic_number('K'), n=3, edge='M1')
        print(f' {E} MeV: {result:12.3f} b')
    
    # print('ALPHA:')
    
    # result = test_PWBA(E1_MeV=1.0, M1_u=2*1.007276, M2_u=241.057, Z1=2, Z2=29, n=1, edge='K')
    # print(f'K-shell: {result:12.3f} b')
    # result = test_PWBA(E1_MeV=1.0, M1_u=2*1.007276, M2_u=241.057, Z1=2, Z2=29, n=2, edge='L1')
    # print(f'L1-shell: {result:12.3f} b')
    # result = test_PWBA(E1_MeV=1.0, M1_u=2*1.007276, M2_u=241.057, Z1=2, Z2=29, n=2, edge='L2')
    # print(f'L2-shell: {result:12.3f} b')
    # result = test_PWBA(E1_MeV=1.0, M1_u=2*1.007276, M2_u=241.057, Z1=2, Z2=29, n=2, edge='L3')
    # print(f'L3-shell: {result:12.3f} b')