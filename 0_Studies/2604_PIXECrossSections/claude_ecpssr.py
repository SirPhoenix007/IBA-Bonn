#!/usr/bin/env python3
"""
ECPSSR-Wirkungsquerschnitte für die K-Schalen-Ionisation durch leichte Ionen.

Implementiert nach T. Mukoyama, "Inner-shell ionization by light-ion impact",
Int. J. PIXE 1, 209 (1991), Abschnitte 2 und 3 (Brandt & Lapicki, PRA 23, 1717 (1981)).

Enthalten:
  * PWBA-Wirkungsquerschnitt (Gl. 2.4-2.15) mit wasserstoffähnlichen K-Formfaktoren
    (Gl. 2.14 / 2.15), numerisch integriert, wahlweise mit
      - Näherungsgrenzen (Qmin = W^2/(4 eta), Qmax = inf)  -> Basis für Gl. 3.22
      - exakten Grenzen (Gl. 2.9 / 2.10)                   -> Basis für Gl. 3.23
  * Bindungs-/Polarisationskorrektur  zeta   (Gl. 3.1-3.9)
  * Coulomb-Ablenkung                 C      (Gl. 3.14)
  * Relativistische Korrektur         m^R    (Gl. 3.16)
  * Energieverlust                    f(z)   (Gl. 3.20, 3.21)

Nur die K-Schale ist implementiert (das Paper liefert nur dafür die
Formfaktoren in geschlossener Form).

Aufruf:
  python ecpssr.py                       # Test: 1 MeV p auf Mn
  python ecpssr.py --Z2 29 --M2 63.546 --EK 8979 --E 2.0
  python ecpssr.py --selftest
"""
import argparse
from math import sqrt, exp, log, log1p, expm1, atan2, pi

from scipy import integrate
from scipy.special import expn

# ---------------------------------------------------------------- Konstanten
RYDBERG_EV = 13.605693          # eV
E2_EV_ANG = 14.399645           # e^2 in eV*Angstrom
A0_SQ_BARN = 2.800285e7         # a0^2 in barn (a0^2 = 2.800285e-17 cm^2)
U_OVER_ME = 1822.888486         # atomare Masseneinheit / Elektronenmasse

# Einige K-Bindungsenergien (eV, Bearden & Burr) und Atommassen (u)
ELEMENTS = {
    22: ("Ti", 4966.0, 47.867), 23: ("V", 5465.0, 50.942),
    24: ("Cr", 5989.0, 51.996), 25: ("Mn", 6539.0, 54.938),
    26: ("Fe", 7112.0, 55.845), 27: ("Co", 7709.0, 58.933),
    28: ("Ni", 8333.0, 58.693), 29: ("Cu", 8979.0, 63.546),
    30: ("Zn", 9659.0, 65.380),
}


# ------------------------------------------------- K-Schalen-Formfaktor (2.14/2.15)
def form_factor_K(W, Q):
    """|F_WK(Q)|^2 pro dW (dimensionslose Einheiten)."""
    if W > 1.0 + 1e-12:                                   # Gl. (2.14), k reell
        k2 = W - 1.0
        k = sqrt(k2)
        den = ((Q - k2 + 1.0) ** 2 + 4.0 * k2) ** 3
        pref = 128.0 * Q * (Q + (k2 + 1.0) / 3.0) / den
        ang = atan2(2.0 * k, Q - k2 + 1.0)                # Zweig in (0, pi)
        return pref * exp(-2.0 * ang / k) / (-expm1(-2.0 * pi / k))
    t2 = max(1.0 - W, 0.0)                                # Gl. (2.15), k = i t
    t = sqrt(t2)
    B = Q + (t + 1.0) ** 2
    den = ((Q + t2 + 1.0) ** 2 - 4.0 * t2) ** 3
    pref = 128.0 * Q * (Q - t2 / 3.0 + 1.0 / 3.0) / den
    if t < 1e-7:
        expo = -4.0 / (Q + 1.0)
    else:
        expo = log1p(-4.0 * t / B) / t                    # ln[(Q+(t-1)^2)/(Q+(t+1)^2)]/t
    return pref * exp(expo)


def _q_limits(W, eta, mu_m):
    """Exakte Grenzen Qmin, Qmax nach Gl. (2.9), (2.10); mu_m = M(reduziert)/m."""
    x = W / (mu_m * eta)
    if x >= 1.0:
        s = 0.0
    else:
        s = sqrt(1.0 - x)
    qmin = W * W / (eta * (1.0 + s) ** 2)                 # numerisch stabile Form
    qmax = (mu_m ** 2) * eta * (1.0 + s) ** 2
    return qmin, qmax


def _inner(W, eta, mu_m):
    """int dQ/Q^2 |F|^2  (in ln Q integriert)."""
    if mu_m is None:                                       # Näherungsgrenzen
        qmin, qmax = W * W / (4.0 * eta), 1e12
    else:
        qmin, qmax = _q_limits(W, eta, mu_m)
    if qmax <= qmin:
        return 0.0
    f = lambda u: form_factor_K(W, exp(u)) / exp(u)
    a, b = log(qmin), log(qmax)
    # Maximum des Integranden liegt grob bei Q ~ 1; Stützpunkt dort setzen
    pts = [p for p in (0.0,) if a < p < b]
    val, _ = integrate.quad(f, a, b, points=pts or None, limit=300,
                            epsabs=0.0, epsrel=1e-8)
    return val


def f_pwba(eta, theta, mu_m=None):
    """f_K(eta, theta) aus Gl. (2.5); mu_m=None -> Näherungsgrenzen (Wmax = inf)."""
    wmin = theta
    wmax = 3000.0 if mu_m is None else mu_m * eta
    total = 0.0
    # Bereich W < 1 :  W = 1 - t^2, dW = 2 t dt
    if wmin < 1.0:
        w_hi = min(1.0, wmax)
        t_lo, t_hi = sqrt(1.0 - w_hi), sqrt(1.0 - wmin)
        g = lambda t: 2.0 * t * _inner(1.0 - t * t, eta, mu_m)
        v, _ = integrate.quad(g, t_lo, t_hi, limit=200, epsabs=0.0, epsrel=1e-7)
        total += v
    # Bereich W > 1 :  W = 1 + k^2, dW = 2 k dk
    if wmax > 1.0:
        kmax = sqrt(wmax - 1.0)
        edges = [0.0]
        e = 0.25
        while e < kmax:
            edges.append(e)
            e *= 2.0
        edges.append(kmax)
        h = lambda k: 2.0 * k * _inner(1.0 + k * k, eta, mu_m)
        for lo, hi in zip(edges[:-1], edges[1:]):
            v, _ = integrate.quad(h, lo, hi, limit=200, epsabs=0.0, epsrel=1e-7)
            total += v
    return total


def sigma_pwba(Z1, Z2s, eta, theta, mu_m=None):
    """PWBA-K-Schalenquerschnitt in barn, Gl. (2.11): sigma = sigma0 * f / eta."""
    sigma0 = 4.0 * pi * 2.0 * Z1 ** 2 * A0_SQ_BARN / Z2s ** 4    # 2j2+1 = 2
    return sigma0 * f_pwba(eta, theta, mu_m) / eta


# --------------------------------------------------------- Korrekturfunktionen
def g_K(xi):                                            # Gl. (3.5)
    num = (1 + 9 * xi + 31 * xi**2 + 98 * xi**3 + 12 * xi**4
           + 25 * xi**5 + 4.2 * xi**6 + 0.515 * xi**7)
    return num / (1 + xi) ** 9


def I_func(x):                                          # Funktion I(x) unter Gl. (3.8)
    if x <= 0.0:
        return float("inf")
    if x <= 0.035:
        return 0.75 * pi * (log(1.0 / x**2) - 1.0)
    if x <= 3.1:
        return exp(-2 * x) / (0.031 + 0.210 * x**0.5 + 0.005 * x
                              - 0.069 * x**1.5 + 0.324 * x**2)
    if x <= 11.0:
        return exp(-2 * x) / x**1.6
    return 0.0


def zeta_K(Z1, Z2s, theta, xi, n=1, c=1.5):            # Gl. (3.8), (3.9)
    h = 2.0 * n / (theta * xi**3) * I_func(c * n / xi)
    return 1.0 + 2.0 * Z1 / (Z2s * theta) * (g_K(xi) - h)


def coulomb_factor(arg, l2=0):                          # Gl. (3.14)
    nu = 10 + 2 * l2
    return (9 + 2 * l2) * expn(nu, arg)


def f_energy_loss(z, nu=9):                             # Gl. (3.20)
    return ((nu * z - 1.0) * (1 + z) ** nu
            + (nu * z + 1.0) * (1 - z) ** nu) / (2.0**nu * (nu - 1.0))


def m_relativistic(Z2s, xi, n=1):                       # Gl. (3.16), K-Schale
    y = 0.40 * (Z2s / 137.0) ** 2 / (n * xi)
    return sqrt(1.0 + 1.1 * y * y) + y


# ------------------------------------------------------------------ ECPSSR
def ecpssr_K(Z1, M1_u, Z2, M2_u, E1_MeV, EK_eV, verbose=True):
    """
    K-Schalen-Ionisationsquerschnitt (barn) in der ECPSSR-Theorie.
    Z1, M1_u : Projektil-Ladung und -Masse (u)
    Z2, M2_u : Target-Ordnungszahl und -Masse (u)
    E1_MeV   : Projektilenergie (Labor, MeV)
    EK_eV    : gemessene K-Bindungsenergie (eV)
    Rückgabe : dict mit verschiedenen Varianten
    """
    E1 = E1_MeV * 1e6                                    # eV
    M1 = M1_u * U_OVER_ME                                # in Elektronenmassen
    M2 = M2_u * U_OVER_ME
    mu_m = M1 * M2 / (M1 + M2)                           # reduzierte Masse / m

    Z2s = Z2 - 0.3                                       # Slater-Abschirmung K (S.212)
    theta = EK_eV / (Z2s**2 * RYDBERG_EV)                # Gl. (2.8), n = 1
    eta = (E1 / M1) / (Z2s**2 * RYDBERG_EV)              # Gl. (2.3), m/M1 * E1/(Z^2 Ry)
    xi = 2.0 * sqrt(eta) / theta                         # xi_K = 2 n sqrt(eta)/theta

    zeta = zeta_K(Z1, Z2s, theta, xi)                    # Gl. (3.9)
    mR = m_relativistic(Z2s, xi)                         # Gl. (3.16)

    # Coulomb-Ablenkung (3.14): d = Z1 Z2 e^2 / (2 E1),  hbar q0 = E_K / v1
    d = Z1 * Z2 * E2_EV_ANG / (2.0 * E1)                 # Angstrom
    v1_over_v0 = sqrt(E1 / M1 / RYDBERG_EV)              # v1/v0, v0 = e^2/hbar
    q0 = EK_eV / (E2_EV_ANG * v1_over_v0)                # 1/Angstrom
    arg0 = pi * d * q0 * zeta

    # Energieverlust (3.21)
    z = sqrt(1.0 - zeta * EK_eV * M1 / (mu_m * E1))
    arg_el = 2.0 * arg0 / (z * (1.0 + z))                # Argument in Gl. (3.22)

    out = {}
    # (a) reine PWBA, Näherungsgrenzen, Gl. (2.11)
    out["PWBA (Näherungsgrenzen)"] = sigma_pwba(Z1, Z2s, eta, theta)
    # (b) PWBA mit exakten Grenzen
    out["PWBA (exakte Grenzen)"] = sigma_pwba(Z1, Z2s, eta, theta, mu_m)
    # (c) Standard-ECPSSR, Gl. (3.22): f(z) * sigma_PWBA(mR*eta, zeta*theta)
    s_pwba_mod = sigma_pwba(Z1, Z2s, mR * eta, zeta * theta)
    out["ECPSSR Gl. (3.22)"] = coulomb_factor(arg_el) * f_energy_loss(z) * s_pwba_mod
    # (d) ECPSSR mit exakten Grenzen, Gl. (3.23)
    s_epwba_mod = sigma_pwba(Z1, Z2s, mR * eta, zeta * theta, mu_m)
    out["ECPSSR Gl. (3.23)"] = coulomb_factor(arg_el) * s_epwba_mod

    if verbose:
        print(f"Projektil Z1={Z1}, M1={M1_u} u | Target Z2={Z2}, M2={M2_u} u")
        print(f"E1 = {E1_MeV} MeV, E_K = {EK_eV} eV")
        print("--- Parameter ---")
        print(f"  Z2K = {Z2s:.3f}")
        print(f"  theta_K = {theta:.5f}")
        print(f"  eta_K   = {eta:.5f}")
        print(f"  xi_K    = {xi:.5f}")
        print(f"  zeta_K  = {zeta:.5f}   (Bindung/Polarisation)")
        print(f"  m^R     = {mR:.5f}   (relativistisch)")
        print(f"  d q0    = {d*q0:.5f}   -> C(pi d q0 zeta) = {coulomb_factor(arg0):.5f}")
        print(f"  z       = {z:.5f}   -> f(z) = {f_energy_loss(z):.5f}")
        print("--- K-Schalen-Querschnitte (barn) ---")
        for k, v in out.items():
            print(f"  {k:28s} {v:12.4f} b")
    return out


# --------------------------------------------------------------- Selbsttest
def selftest():
    """Oszillatorstärke des Kontinuums von H(1s) bei Q -> 0 muss 0.4350 ergeben."""
    Qs = 1e-8
    fW = lambda W: W * form_factor_K(W, Qs) / Qs
    # W = 1 + k^2
    v, _ = integrate.quad(lambda k: 2 * k * fW(1 + k * k), 0, 1, limit=200)
    v2, _ = integrate.quad(lambda k: 2 * k * fW(1 + k * k), 1, 1e4, limit=500)
    print(f"Kontinuums-Oszillatorstärke H(1s): {v + v2:.4f}  (Soll: 0.4350)")


def main():
    ap = argparse.ArgumentParser(description="K-Schalen ECPSSR (Mukoyama 1991)")
    ap.add_argument("--Z1", type=int, default=1)
    ap.add_argument("--M1", type=float, default=1.007276, help="Projektilmasse in u")
    ap.add_argument("--Z2", type=int, default=25)
    ap.add_argument("--M2", type=float, default=None, help="Targetmasse in u")
    ap.add_argument("--EK", type=float, default=None, help="K-Bindungsenergie in eV")
    ap.add_argument("--E", type=float, default=1.0, help="Projektilenergie in MeV")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        selftest()
        return
    name, ek, m2 = ELEMENTS.get(a.Z2, (None, None, None))
    EK = a.EK if a.EK is not None else ek
    M2 = a.M2 if a.M2 is not None else m2
    if EK is None or M2 is None:
        raise SystemExit("Für dieses Z2 bitte --M2 und --EK angeben.")
    ecpssr_K(a.Z1, a.M1, a.Z2, M2, a.E, EK)


if __name__ == "__main__":
    main()
