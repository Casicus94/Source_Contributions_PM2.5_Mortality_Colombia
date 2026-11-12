#!/usr/bin/env python3
"""
Table 1 (mortality percent contributions by sector and department), GeoHealth revision.

Inputs (Brendan's files, folder for_alejandro_Sep29_2026/, all on the same 429x429 lon-lat grid):
  mort_total.nc                         mort_dist, mort_dist_min, mort_dist_max  (deaths per grid cell, FUSION, ACAG exposure)
  PM2_5_<RUN>_ymean_lonlat.nc           annual mean PM2.5 of each WRF-Chem run
and Departments_gridded.nc (shape_id per grid cell) + departments_list.csv (names).

Two attribution methods:
  --method brendan    reproduces calculate_fractions.py (Jan 2025):
                      sectors (incl. fires) scaled to  NoFIRES - NoANTHRO_FIRES  (anthropogenic only),
                      natural = reference - that total  -> the PM2.5 from fires ends up in 'Oth'.
  --method corrected  (used in the revised manuscript):
                      sectors (8 anthropogenic + fires) scaled to  reference - NoANTHRO_FIRES  (anthropogenic + fires),
                      Oth = NoANTHRO_FIRES / reference.
Sector contribution: NoFIRES - <sector run> (sector runs have fires off); fires: reference - NoFIRES.

Deaths per department and sector = sum over grid cells of deaths x fraction.
  --floor   rounds each department-sector death count DOWN to an integer before summing (Brendan's convention);
            by default no rounding is applied until display (revised manuscript).
Also prints the closure of the zero-out experiments (sum of individual contributions / combined contribution).

Usage (from the MPIC_project folder):
  python Scripts_submit/make_table1.py --method brendan --floor      # -> reproduces the submitted Table 1 (12,817 deaths)
  python Scripts_submit/make_table1.py --method corrected            # -> revised Table 1 (12,921 deaths)
"""
import argparse, os
import numpy as np, pandas as pd, netCDF4 as nc

SEC = {'Agr': 'Agriculture', 'Avi': 'Aviation', 'Ene': 'Energy_Use', 'Ind': 'Industry',
       'ORT': 'Off_Road_Transport', 'RT': 'Road_Transport', 'Shp': 'Shipping', 'Wst': 'Waste'}
KEYS = ['Agr', 'Avi', 'Ene', 'Ind', 'BB', 'ORT', 'RT', 'Shp', 'Wst', 'Oth']


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', default='.')
    ap.add_argument('--data', default='for_alejandro_Sep29_2026')
    ap.add_argument('--method', choices=['brendan', 'corrected'], default='corrected')
    ap.add_argument('--floor', action='store_true')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    R = lambda p: os.path.join(a.root, p)
    B = lambda p: os.path.join(a.root, a.data, p)
    L = lambda run: np.ma.filled(nc.Dataset(B(f'PM2_5_{run}_ymean_lonlat.nc'))['PM2_5_DRY'][:], np.nan).astype(float).squeeze()

    m = nc.Dataset(B('mort_total.nc'))
    D, Dmin, Dmax = (np.ma.filled(m[v][:], 0.).astype(float) for v in ('mort_dist', 'mort_dist_min', 'mort_dist_max'))
    sid = np.array(nc.Dataset(R('Departments_gridded.nc'))['shape_id'][:])
    lst = pd.read_csv(R('departments_list.csv')); names = dict(zip(lst.iloc[:, 1].astype(int), lst.iloc[:, 0]))

    ref, nf, naf = L('TOTALS'), L('TOTALS_No_FIRES'), L('TOTALS_No_ANTHRO_FIRES')
    S = {k: np.clip(nf - L(v), 0, None) for k, v in SEC.items()}
    S['BB'] = np.clip(ref - nf, 0, None)
    tot = sum(S.values())
    if a.method == 'brendan':
        target = np.clip(nf - naf, 0, None)                 # anthropogenic only (the error)
        natural = np.clip(ref - target, 0, None)            # includes the fire PM2.5
    else:
        target = np.clip(ref - naf, 0, None)                # anthropogenic + fires
        natural = np.clip(naf, 0, None)
    norm = np.where(tot == 0, 1.0, target / np.where(tot == 0, 1, tot))
    C = {k: S[k] * norm for k in S}; C['Oth'] = natural
    total = sum(C.values())
    F = {k: np.nan_to_num(C[k] / total) for k in KEYS}

    # closure of the zero-out experiments (deaths-weighted), independent of the method
    k = (sid > 0) & (D > 0) & np.isfinite(ref)
    comb = np.clip(ref - naf, 0, None)
    r = np.where(comb > 0, tot / np.where(comb > 0, comb, 1), np.nan)[k]; w = D[k]; ok = np.isfinite(r)
    o = np.argsort(r[ok]); cw = np.cumsum(w[ok][o]) / w[ok].sum()
    print('closure (sum of sectors+fires / combined): deaths-weighted mean %.3f, p5 %.3f, p95 %.3f'
          % ((w[ok] * r[ok]).sum() / w[ok].sum(), r[ok][o][np.searchsorted(cw, .05)], r[ok][o][np.searchsorted(cw, .95)]))

    rows = []
    for i, n in sorted(names.items()):
        msk = sid == i
        if D[msk].sum() == 0:
            continue                                        # San Andres: outside the inner domain
        dd = {x: (D * F[x])[msk].sum() for x in KEYS}
        if a.floor:
            dd = {x: np.floor(v) for x, v in dd.items()}
        tot_d = sum(dd.values())
        row = {'dept': n, **{x: 100 * dd[x] / tot_d if tot_d > 0 else np.nan for x in KEYS},
               'Tot': tot_d, 'min': Dmin[msk].sum(), 'max': Dmax[msk].sum()}
        row['_deaths'] = dd; rows.append(row)
    nat = {x: sum(r_['_deaths'][x] for r_ in rows) for x in KEYS}; nt = sum(nat.values())
    rows.append({'dept': 'COLOMBIA', **{x: 100 * nat[x] / nt for x in KEYS}, 'Tot': nt,
                 'min': Dmin[sid > 0].sum(), 'max': Dmax[sid > 0].sum()})
    T = pd.DataFrame([{kk: vv for kk, vv in r_.items() if kk != '_deaths'} for r_ in rows]).set_index('dept')
    out = a.out or f"Table1_{a.method}{'_floor' if a.floor else ''}.csv"
    T.round(2).to_csv(out)
    pd.set_option('display.width', 200)
    print(T.round(2).to_string()); print('wrote', out)


if __name__ == '__main__':
    main()
