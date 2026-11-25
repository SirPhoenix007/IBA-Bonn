"""
Lenz-Jensen electron density rho_LJ(r) and normalisation check.

Sources
-------
Jensen, Z. Phys. 77, 722 (1932): eqs. (4), (5a), (9), (11), Table 1
Sannikov, JETP 34, 1650 (1958): eq. (6)

Units: Hartree atomic units. a_H = 1 (length), e = 1 (charge).
    length  : a_H  (Bohr radius, 0.529177 Angstrom)
    density : electrons / a_H^3
Everything below is therefore dimensionless numbers in these units.
"""
import numpy as np
from math import comb, factorial

A_H_ANGSTROM = 0.529177210903  # a_H in Angstrom (only used for printing)

# ---------------------------------------------------------------------------
# Jensen Table 1:  Z, N, c1, 1/lambda   (1/lambda is given in units Z^(1/3)/a_H)
# ---------------------------------------------------------------------------
TABLE1 = {
    "Na+":     dict(Z=11, N=10, c1=0.298, inv_lam=12.87),
    "K+":      dict(Z=19, N=18, c1=0.285, inv_lam=12.04),
    "Rb+":     dict(Z=37, N=36, c1=0.275, inv_lam=11.47),
    "neutral": dict(Z=37, N=37, c1=0.265, inv_lam=10.91),  # Z arbitrary, scales out
    "Br-":     dict(Z=35, N=36, c1=0.254, inv_lam=10.33),
    "Cl-":     dict(Z=17, N=18, c1=0.243, inv_lam=9.76),
}


def P0(c1):
    """Jensen eq. (11): P0 = 2 * sum_n C(3,n) (n+2)! c1^n   (dimensionless).
    Sannikov's 4P equals P0 with C = c1."""
    return 2.0 * sum(comb(3, n) * factorial(n + 2) * c1**n for n in range(4))


def lambda_jensen(Z, inv_lam):
    """Jensen's length lambda [a_H].
    Table 1 gives 1/lambda in units Z^(1/3)/a_H  ->  lambda = 1/(Z^(1/3) * inv_lam)."""
    return 1.0 / (Z ** (1.0 / 3.0) * inv_lam)


def rho_LJ(r, N, c1, lam):
    """Electron number density [electrons / a_H^3] at radius r [a_H].

    Jensen:   n(r) = N / (4 pi lam^3 P0) * exp(-x_J) * x_J^-3 * (1 + c1 x_J)^3
              with x_J = sqrt(r/lam)                       (dimensionless)
    Sannikov: identical, written with x_S = r/lam = x_J^2,
              x_S^(-3/2) exp(-x_S^(1/2)) (1 + C x_S^(1/2))^3.
    lam [a_H] plays the role of mu/lambda_S in Sannikov's notation.
    """
    xJ = np.sqrt(r / lam)                                   # dimensionless
    pref = N / (4.0 * np.pi * lam**3 * P0(c1))              # [1/a_H^3]
    return pref * np.exp(-xJ) * xJ**-3 * (1.0 + c1 * xJ) ** 3


def check_normalisation(name, smax=40.0, npts=200001):
    p = TABLE1[name]
    Z, N, c1 = p["Z"], p["N"], p["c1"]
    lam = lambda_jensen(Z, p["inv_lam"])                    # [a_H]

    # Uniform grid in s = sqrt(x_S) = x_J (dimensionless).
    # r = lam * s^2 [a_H],  dr = 2 lam s ds [a_H].
    # This removes the r^(-3/2) cusp at the nucleus. s = 0 is r = 0.
    s = np.linspace(0.0, smax, npts)
    r = lam * s**2                                          # [a_H]
    dr_ds = 2.0 * lam * s                                   # [a_H]

    integrand = np.zeros_like(s)                            # d(charge)/ds [electrons]
    rho = rho_LJ(r[1:], N, c1, lam)                         # [e/a_H^3], skip r=0 (rho -> inf)
    integrand[1:] = 4.0 * np.pi * r[1:] ** 2 * rho * dr_ds[1:]
    # (integrand -> 0 at s = 0 since r^2 rho dr/ds ~ s)

    # Simpson-free trapezoid is fine at this resolution
    Q_num = np.trapezoid(integrand, s)                      # [electrons]
    # Analytic: integral = (2N/P0) * sum_n C(3,n) c1^n (n+2)! = N exactly
    Q_cum = np.cumsum(integrand) * (s[1] - s[0])            # N(<r) vs r

    print(f"--- {name}:  Z={Z}, N={N}, c1={c1}, 1/lambda={p['inv_lam']} [Z^(1/3)/a_H]")
    print(f"lambda_J        = {lam:.5f} a_H = {lam*A_H_ANGSTROM:.5f} Angstrom")
    print(f"P0(c1)          = {P0(c1):.5f}")
    print(f"x_S = r/lambda_J: x_S = 100 -> r = {100*lam:.3f} a_H = {100*lam*A_H_ANGSTROM:.3f} Angstrom")
    print(f"grid: s in [0,{smax}] -> x_S up to {smax**2:.0f}, r up to {lam*smax**2:.1f} a_H")
    print(f"integral rho 4 pi r^2 dr = {Q_num:.6f}   (target N = {N})")
    # Charge missing if x_S were cut at 100 (s = 10):
    i100 = np.searchsorted(s, 10.0)
    print(f"charge inside x_S<=100   = {Q_cum[i100]:.4f}  ({100*Q_cum[i100]/N:.2f} % of N)")
    return r, rho


if __name__ == "__main__":
    for name in ["Rb+", "Br-", "neutral"]:
        check_normalisation(name)
        print()
