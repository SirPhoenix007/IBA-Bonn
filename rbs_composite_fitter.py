#!/usr/bin/env python3
"""
falling_composite_fit.py
════════════════════════
Fit a continuous composite function to 1-D data:

         ⎧ A_L · exp(-(x-x_L)² / (2σ²))        x ≤ x_L
  f(x) = ⎨ falling polynomial of degree 1 or 2  x_L < x < x_R
         ⎩ A_R · exp(-(x-x_R)² / (2σ²))        x ≥ x_R

Model constraints
─────────────────
- The two Gaussian halves have identical width σ.
- The left and right Gaussian heights A_L and A_R may differ.
- The connecting polynomial is degree 1 or degree 2 only.
- The connecting polynomial is strictly falling.

Parameter order
───────────────
For poly_degree = 1:

    [x_L, x_R, A_L, A_R, sigma]

For poly_degree = 2:

    [x_L, x_R, A_L, A_R, sigma, q]

The degree-1 connector is

    P(t) = A_L + (A_R - A_L)t

The degree-2 connector is

    P(t) = A_L + (A_R - A_L)t + q t(1 - t)

with

    t = (x - x_L)/(x_R - x_L).

Strict falling requires

    x_R > x_L,
    A_L > A_R,

and for the quadratic connector additionally

    |q| < A_L - A_R.

Important note
──────────────
With a degree-1 or degree-2 middle polynomial, the model is continuous but
generally not C¹-continuous at x_L and x_R.

Dependencies
────────────
numpy, scipy, matplotlib
"""

from __future__ import annotations

import warnings
from typing import List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


# ══════════════════════════════════════════════════════════════════════════════
# Composite model factory
# ══════════════════════════════════════════════════════════════════════════════

def _make_composite(poly_degree: int):
    """
    Build composite function:

        left Gaussian half       x <= x_L
        falling linear/quadratic x_L < x < x_R
        right Gaussian half      x >= x_R

    Shared Gaussian width:

        sigma_L = sigma_R = sigma

    Parameter order:

        x_L, x_R, A_L, A_R, sigma[, q]

    For poly_degree = 1:
        middle is the unique line from (x_L, A_L) to (x_R, A_R).

    For poly_degree = 2:
        middle is

            P(t) = A_L + (A_R - A_L)t + q t(1 - t)

        where

            t = (x - x_L)/(x_R - x_L).

    To keep the quadratic strictly falling on t in [0, 1], enforce

        A_L > A_R
        -d < q < d

    where

        d = A_L - A_R > 0.
    """
    if poly_degree not in (1, 2):
        raise ValueError(
            f"poly_degree must be either 1 or 2, got {poly_degree}."
        )

    n_extra = 0 if poly_degree == 1 else 1

    def _f(x, *params):
        x = np.asarray(x, dtype=float)

        x_L, x_R, A_L, A_R, sigma = params[:5]
        q = params[5] if poly_degree == 2 else 0.0

        delta_x = x_R - x_L
        delta_A = A_L - A_R

        # Required for a meaningful strictly falling connector.
        if delta_x <= 0.0 or delta_A <= 0.0 or sigma <= 0.0:
            return np.full_like(x, 1e300)

        # Strict monotonicity condition for quadratic connector:
        #
        # P(t) = A_L - delta_A t + q t(1 - t)
        # P'(t) = -delta_A + q(1 - 2t)
        #
        # For P'(t) < 0 for all t in [0, 1], need:
        #
        #     -delta_A + |q| < 0
        #
        # therefore:
        #
        #     |q| < delta_A
        if poly_degree == 2 and abs(q) >= delta_A:
            return np.full_like(x, 1e300)

        t = (x - x_L) / delta_x
        out = np.empty_like(x)

        # Left half-Gaussian
        mL = t <= 0.0
        if mL.any():
            out[mL] = A_L * np.exp(
                -((x[mL] - x_L) ** 2) / (2.0 * sigma ** 2)
            )

        # Right half-Gaussian
        mR = t >= 1.0
        if mR.any():
            out[mR] = A_R * np.exp(
                -((x[mR] - x_R) ** 2) / (2.0 * sigma ** 2)
            )

        # Middle polynomial
        mM = ~mL & ~mR
        if mM.any():
            tm = t[mM]

            if poly_degree == 1:
                out[mM] = A_L + (A_R - A_L) * tm

            else:
                out[mM] = (
                    A_L
                    + (A_R - A_L) * tm
                    + q * tm * (1.0 - tm)
                )

        return out

    return _f, n_extra


# ══════════════════════════════════════════════════════════════════════════════
# Public fitter class
# ══════════════════════════════════════════════════════════════════════════════

class CompositeGaussPolyFitter:
    """
    Fit a continuous half-Gauss | falling poly | half-Gauss composite to data.

    New model constraints
    ---------------------
    - The two Gaussian halves share one common width sigma.
    - The heights A_L and A_R may differ.
    - The connecting polynomial is degree 1 or 2 only.
    - The connecting polynomial is strictly falling.

    Parameter order
    ---------------
    poly_degree = 1:

        [x_L, x_R, A_L, A_R, sigma]

    poly_degree = 2:

        [x_L, x_R, A_L, A_R, sigma, q]

    The quadratic connector is

        P(t) = A_L + (A_R - A_L)t + q t(1 - t)

    with

        t = (x - x_L)/(x_R - x_L).

    Strict falling requires

        x_R > x_L,
        A_L > A_R,
        |q| < A_L - A_R.

    Note
    ----
    With degree 1 or 2, the model is generally continuous but not C¹-continuous
    at the join points.
    """

    def __init__(self, poly_degree: int = 1) -> None:
        self._func, self._n_extra = _make_composite(poly_degree)
        self.poly_degree = poly_degree
        self.n_params = 5 + self._n_extra

        self.param_names: List[str] = [
            "x_L",
            "x_R",
            "A_L",
            "A_R",
            "sigma",
        ]

        if poly_degree == 2:
            self.param_names.append("q")

        self.popt: Optional[np.ndarray] = None
        self.perr: Optional[np.ndarray] = None
        self.pcov: Optional[np.ndarray] = None

    # ── fit ──────────────────────────────────────────────────────────────────

    def fit(
        self,
        x_data: np.ndarray,
        y_data: np.ndarray,
        p0: Union[List[float], np.ndarray],
        bounds: Tuple = (-np.inf, np.inf),
        sigma_y: Optional[np.ndarray] = None,
        maxfev: int = 20_000,
        **kwargs,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Fit the composite function to x_data and y_data.

        Parameters
        ----------
        x_data, y_data:
            Data arrays of equal length.

        p0:
            Initial guesses in the correct parameter order.

            For poly_degree = 1:

                [x_L, x_R, A_L, A_R, sigma]

            For poly_degree = 2:

                [x_L, x_R, A_L, A_R, sigma, q]

        bounds:
            Lower and upper bounds for the parameters.

        sigma_y:
            Optional y-data uncertainties. If supplied, the fit is weighted.

        maxfev:
            Maximum number of function evaluations.

        Returns
        -------
        popt:
            Best-fit parameter values.

        perr:
            One-standard-deviation parameter uncertainties.
        """
        p0 = list(p0)

        if len(p0) != self.n_params:
            raise ValueError(
                f"p0 has {len(p0)} element(s), but {self.n_params} are required.\n"
                f"Parameter order: {self.param_names}"
            )

        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            popt, pcov = curve_fit(
                self._func,
                np.asarray(x_data, dtype=float),
                np.asarray(y_data, dtype=float),
                p0=p0,
                bounds=bounds,
                sigma=sigma_y,
                absolute_sigma=(sigma_y is not None),
                maxfev=maxfev,
                **kwargs,
            )

        diag = np.diag(pcov)

        if np.any(~np.isfinite(diag)):
            warnings.warn(
                "Covariance matrix contains inf/nan; fit may be poorly constrained.",
                RuntimeWarning,
                stacklevel=2,
            )

        if np.any(diag < 0):
            warnings.warn(
                "Negative covariance diagonal entries; uncertainties are unreliable.",
                RuntimeWarning,
                stacklevel=2,
            )

        perr = np.sqrt(np.abs(diag))

        self.popt = popt
        self.perr = perr
        self.pcov = pcov

        return popt, perr

    # ── evaluate ─────────────────────────────────────────────────────────────

    def evaluate(
        self,
        x: np.ndarray,
        params: Optional[np.ndarray] = None,
    ) -> np.ndarray:
        """
        Evaluate the composite model at positions x.
        """
        if params is None:
            if self.popt is None:
                raise RuntimeError("Run .fit() first, or supply params explicitly.")
            params = self.popt

        return self._func(np.asarray(x, dtype=float), *params)

    # ── uncertainty band ─────────────────────────────────────────────────────

    def uncertainty_band(
        self,
        x: np.ndarray,
        n_samples: int = 500,
        ci: float = 0.68,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Estimate a pointwise uncertainty band via Monte Carlo propagation.
        """
        if self.popt is None or self.pcov is None:
            raise RuntimeError("Run .fit() first.")

        rng = np.random.default_rng()

        samples = rng.multivariate_normal(
            self.popt,
            self.pcov,
            size=n_samples,
            check_valid="warn",
        )

        x = np.asarray(x, dtype=float)

        curves = np.array([
            self._func(x, *s) for s in samples
        ])

        lo_pct = 100.0 * (1.0 - ci) / 2.0
        hi_pct = 100.0 - lo_pct

        return (
            np.percentile(curves, lo_pct, axis=0),
            np.percentile(curves, hi_pct, axis=0),
        )

    # ── print results ────────────────────────────────────────────────────────

    def print_results(self) -> None:
        """
        Print a formatted summary table of all fit results.
        """
        if self.popt is None:
            print("No fit performed yet.")
            return

        w = max(len(n) for n in self.param_names)
        bar = "─" * (w + 42)

        print(f"\n{bar}")
        print(
            f"  Composite Gauss–Poly Fit "
            f"(poly_degree = {self.poly_degree})"
        )
        print(bar)
        print(f"  {'Parameter':<{w}}   {'Value':>16}   {'± 1σ':>16}")
        print(f"  {'─' * (w + 36)}")

        for name, val, err in zip(self.param_names, self.popt, self.perr):
            print(f"  {name:<{w}}   {val:>16.6g}   ±{err:>15.6g}")

        print(f"{bar}\n")

    # ── plot ─────────────────────────────────────────────────────────────────

    def plot(
        self,
        x_data: np.ndarray,
        y_data: np.ndarray,
        *,
        n_plot: int = 2000,
        show_band: bool = False,
        show_residuals: bool = True,
        show_join_pts: bool = True,
        title: str = "Composite Gaussian–Polynomial Fit",
    ) -> plt.Figure:
        """
        Plot data, fit, optional uncertainty band, and optional residuals.
        """
        if self.popt is None:
            raise RuntimeError("Run .fit() first.")

        x_data = np.asarray(x_data, dtype=float)
        y_data = np.asarray(y_data, dtype=float)

        x_plot = np.linspace(x_data.min(), x_data.max(), n_plot)
        y_fit = self.evaluate(x_plot)

        if show_residuals:
            fig, (ax, ax_res) = plt.subplots(
                2,
                1,
                figsize=(9, 6.5),
                sharex=True,
                gridspec_kw={"height_ratios": [3, 1]},
            )
        else:
            fig, ax = plt.subplots(figsize=(9, 4.5))

        ax.scatter(
            x_data,
            y_data,
            s=14,
            color="steelblue",
            alpha=0.6,
            label="Data",
            zorder=3,
        )

        ax.plot(
            x_plot,
            y_fit,
            color="crimson",
            lw=2,
            label="Fit",
            zorder=4,
        )

        if show_band:
            y_lo, y_hi = self.uncertainty_band(x_plot)
            ax.fill_between(
                x_plot,
                y_lo,
                y_hi,
                color="crimson",
                alpha=0.15,
                label="1σ band",
            )

        if show_join_pts:
            x_L, x_R, A_L, A_R = self.popt[:4]

            for xv, av, lbl in [
                (x_L, A_L, f"x_L = {x_L:.4g}"),
                (x_R, A_R, f"x_R = {x_R:.4g}"),
            ]:
                ax.axvline(xv, color="gray", ls="--", lw=0.9, alpha=0.7)
                ax.scatter(
                    [xv],
                    [av],
                    marker="D",
                    color="orange",
                    s=70,
                    zorder=5,
                    label=lbl,
                )

        ax.set_ylabel("y")
        ax.set_title(title)
        ax.legend(framealpha=0.88, fontsize=8)
        ax.grid(alpha=0.2)

        if show_residuals:
            residuals = y_data - self.evaluate(x_data)
            ax_res.scatter(
                x_data,
                residuals,
                s=8,
                color="steelblue",
                alpha=0.5,
            )
            ax_res.axhline(0, color="crimson", lw=0.9)
            ax_res.set_xlabel("x")
            ax_res.set_ylabel("Residuals")
            ax_res.grid(alpha=0.2)
        else:
            ax.set_xlabel("x")

        fig.tight_layout()
        return fig


# ══════════════════════════════════════════════════════════════════════════════
# Example 1: linear falling connector
# ══════════════════════════════════════════════════════════════════════════════

def demo_linear() -> None:
    """
    Demonstration using a degree-1 strictly falling connector.
    """
    rng = np.random.default_rng(42)

    fitter = CompositeGaussPolyFitter(poly_degree=1)

    # True parameter order:
    # x_L, x_R, A_L, A_R, sigma
    true_params = np.array([
        25.0,
        75.0,
        800.0,
        600.0,
        8.0,
    ])

    x = np.linspace(0.0, 100.0, 500)
    y_clean = fitter.evaluate(x, true_params)
    y = y_clean + rng.normal(0.0, 20.0, size=x.size)

    # Initial guess:
    # x_L, x_R, A_L, A_R, sigma
    p0 = [
        22.0,
        78.0,
        750.0,
        550.0,
        10.0,
    ]

    lo = [
        5.0,
        50.0,
        0.0,
        0.0,
        1.0,
    ]

    hi = [
        45.0,
        95.0,
        5000.0,
        5000.0,
        40.0,
    ]

    popt, perr = fitter.fit(x, y, p0=p0, bounds=(lo, hi))
    fitter.print_results()

    print("Linear connector parameter comparison:")
    w = max(len(n) for n in fitter.param_names)

    for name, tv, fv, fe in zip(fitter.param_names, true_params, popt, perr):
        pull = (fv - tv) / fe if fe > 0 else float("nan")
        print(
            f"  {name:<{w}}  true={tv:>10.4g}  "
            f"fit={fv:>10.4g} ±{fe:.4g}   pull={pull:+.2f}σ"
        )

    fig = fitter.plot(
        x,
        y,
        show_band=True,
        show_residuals=True,
        title="Demo — Shared-Width Gaussians with Linear Falling Connector",
    )

    plt.savefig(
        "falling_composite_linear_demo.png",
        dpi=150,
        bbox_inches="tight",
    )

    plt.show()


# ══════════════════════════════════════════════════════════════════════════════
# Example 2: quadratic falling connector
# ══════════════════════════════════════════════════════════════════════════════

def demo_quadratic() -> None:
    """
    Demonstration using a degree-2 strictly falling connector.
    """
    rng = np.random.default_rng(123)

    fitter = CompositeGaussPolyFitter(poly_degree=2)

    # True parameter order:
    # x_L, x_R, A_L, A_R, sigma, q
    #
    # Strict falling requires:
    #
    #     |q| < A_L - A_R
    #
    # Here A_L - A_R = 250, so q = 100 is valid.
    true_params = np.array([
        25.0,
        75.0,
        850.0,
        600.0,
        8.0,
        100.0,
    ])

    x = np.linspace(0.0, 100.0, 500)
    y_clean = fitter.evaluate(x, true_params)
    y = y_clean + rng.normal(0.0, 20.0, size=x.size)

    # Initial guess:
    # x_L, x_R, A_L, A_R, sigma, q
    p0 = [
        22.0,
        78.0,
        800.0,
        550.0,
        10.0,
        0.0,
    ]

    lo = [
        5.0,
        50.0,
        0.0,
        0.0,
        1.0,
        -5000.0,
    ]

    hi = [
        45.0,
        95.0,
        5000.0,
        5000.0,
        40.0,
        5000.0,
    ]

    popt, perr = fitter.fit(x, y, p0=p0, bounds=(lo, hi))
    fitter.print_results()

    print("Quadratic connector parameter comparison:")
    w = max(len(n) for n in fitter.param_names)

    for name, tv, fv, fe in zip(fitter.param_names, true_params, popt, perr):
        pull = (fv - tv) / fe if fe > 0 else float("nan")
        print(
            f"  {name:<{w}}  true={tv:>10.4g}  "
            f"fit={fv:>10.4g} ±{fe:.4g}   pull={pull:+.2f}σ"
        )

    fig = fitter.plot(
        x,
        y,
        show_band=True,
        show_residuals=True,
        title="Demo — Shared-Width Gaussians with Quadratic Falling Connector",
    )

    plt.savefig(
        "falling_composite_quadratic_demo.png",
        dpi=150,
        bbox_inches="tight",
    )

    plt.show()


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("\nRunning linear demo...\n")
    demo_linear()

    print("\nRunning quadratic demo...\n")
    demo_quadratic()