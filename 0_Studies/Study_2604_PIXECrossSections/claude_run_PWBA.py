"""
Run PWBA ionisation cross sections for all populated shells of an element.

Uses PWBA_ECPSSR_cross.py (and M_coeff_matrix.py) from the same folder.
Only shells with a non-zero fluorescence yield in xraydb are calculated,
and only those for which a form factor is implemented (K, L1, L2, L3, M1).

Example:
    python run_PWBA.py Cu --projectile p --emin 0.5 --emax 5 --npoints 10
"""
import argparse
import csv
import os
import time
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed

import numpy as np
import xraydb as xy

from PWBA_ECPSSR_cross import test_PWBA

# projectile: (charge Z1, mass in u)
PROJECTILES = {
    'p': (1, 1.007276),
    'd': (1, 2.013553),
    'a': (2, 4.001506),
}

# principal quantum number per edge, and the edges that have form factors implemented
EDGE_N = {'K': 1, 'L1': 2, 'L2': 2, 'L3': 2, 'M1': 3}


def PWBA_point(E_MeV: float, M1_u: float, M2_u: float, Z1: int, Z2: int, edge: str):
    """One (edge, energy) cross section in barn."""
    n = EDGE_N[edge]
    return edge, E_MeV, test_PWBA(E1_MeV=E_MeV, M1_u=M1_u, M2_u=M2_u,
                                  Z1=Z1, Z2=Z2, n=n, edge=edge)


def populated_edges(elem: str):
    """Edges with non-zero fluorescence yield that this code can handle."""
    element_edges = xy.xray_edges(element=elem)
    edges = [e for e in element_edges.keys()
             if element_edges[e].fyield != 0 and e in EDGE_N]
    skipped = [e for e in element_edges.keys()
               if element_edges[e].fyield != 0 and e not in EDGE_N]
    return edges, [element_edges[e].fyield for e in edges], skipped


def PWBA_single_elem(elem: str, energy_range, projectile='p', workers=None,
                     use_processes=False):
    Z2 = xy.atomic_number(elem)
    M2_u = xy.atomic_mass(elem)
    Z1, M1_u = PROJECTILES[projectile]

    edges, fyields, skipped = populated_edges(elem)
    print(f'Element {elem} (Z={Z2}, M={M2_u:.4f} u)')
    print('Populated shells :', edges)
    print('Fluorescence yld :', [round(f, 5) for f in fyields])
    if skipped:
        print('Skipped (no form factor implemented):', skipped)

    results = {e: {} for e in edges}
    tasks = [(e, E) for e in edges for E in energy_range]
    Executor = ProcessPoolExecutor if use_processes else ThreadPoolExecutor
    done, t0 = 0, time.time()

    with Executor(max_workers=workers) as pool:
        futures = {pool.submit(PWBA_point, E, M1_u, M2_u, Z1, Z2, e): (e, E)
                   for e, E in tasks}
        for fut in as_completed(futures):
            e, E = futures[fut]
            try:
                _, _, sigma = fut.result()
            except Exception as exc:  # keep going if one point fails
                print(f'  {e} @ {E:.3f} MeV failed: {exc}')
                sigma = np.nan
            results[e][E] = sigma
            done += 1
            print(f'  {e:>2} @ {E:7.3f} MeV: {sigma:14.6g} b '
                  f'[{done}/{len(tasks)}, {100 * done / len(tasks):.1f} %]')

    print(f'Finished in {time.time() - t0:.1f} s')
    # return sorted arrays per edge
    return {e: np.array([results[e][E] for E in energy_range]) for e in edges}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('element', help='target element symbol, e.g. Cu')
    ap.add_argument('--projectile', choices=PROJECTILES, default='p')
    ap.add_argument('--emin', type=float, default=1.0, help='MeV')
    ap.add_argument('--emax', type=float, default=1.0, help='MeV')
    ap.add_argument('--npoints', type=int, default=1)
    ap.add_argument('--log', action='store_true', help='log-spaced energies')
    ap.add_argument('--workers', type=int, default=os.cpu_count())
    ap.add_argument('--processes', action='store_true',
                    help='use processes instead of threads (true parallelism with the GIL)')
    ap.add_argument('--out', default=None, help='CSV output file')
    ap.add_argument('--plot', action='store_true')
    args = ap.parse_args()

    if args.log and args.emin > 0:
        energy_range = np.geomspace(args.emin, args.emax, args.npoints)
    else:
        energy_range = np.linspace(args.emin, args.emax, args.npoints)

    res = PWBA_single_elem(args.element, energy_range, args.projectile,
                           workers=args.workers, use_processes=args.processes)

    out = args.out or f'PWBA_{args.element}_{args.projectile}.csv'
    with open(out, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['E_MeV'] + [f'sigma_{e}_barn' for e in res])
        for i, E in enumerate(energy_range):
            w.writerow([E] + [res[e][i] for e in res])
    print(f'Wrote {out}')

    if args.plot and len(energy_range) > 1:
        import matplotlib.pyplot as plt
        for e, s in res.items():
            plt.loglog(energy_range, s, marker='o', label=e)
        plt.xlabel('Energy (MeV)')
        plt.ylabel('PWBA cross section (b)')
        plt.title(f'{args.projectile} on {args.element}')
        plt.legend()
        plt.grid(True, which='both', alpha=0.3)
        plt.show()


if __name__ == '__main__':
    main()
