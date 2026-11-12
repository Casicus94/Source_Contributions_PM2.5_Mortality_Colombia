#!/usr/bin/env python3
"""
Reviewer 3, major comment 3: sensitivity of the national sectoral mortality shares to the (non-sector-neutral) model bias.
Annual version. q = model/ACAG annual mean PM2.5 (Jan-Nov 2019, model x3.6) in each grid cell, mapped to the mortality grid.
Baseline sector fractions f_k (corrected method, Table 1). Scenarios (applied where the model is biased):
  S1  q<1: missing mass is natural (e.g. biogenic SOA)       -> all sector contributions keep their absolute value, rest -> Oth
  S2  q<1: missing mass is anthropogenic (inventory deficit) -> missing mass split among anthropogenic sectors in proportion to f_k
  S3  q>1: excess is fire smoke                              -> excess removed from BB first (then proportionally from the rest)
  S1+S3, S2+S3 combined.
National share = deaths-weighted mean of the fractions (deaths from the ACAG exposure, unchanged).
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd, netCDF4 as nc
from common import D, sid, L, B, P
sec = {'Agr':'Agriculture','Avi':'Aviation','Ene':'Energy_Use','Ind':'Industry','ORT':'Off_Road_Transport','RT':'Road_Transport','Shp':'Shipping','Wst':'Waste'}
ANTH = list(sec)
ref = L('TOTALS').squeeze(); nf = L('TOTALS_No_FIRES').squeeze(); naf = L('TOTALS_No_ANTHRO_FIRES').squeeze()
S = {k: np.clip(nf - L(v).squeeze(), 0, None) for k, v in sec.items()}; S['BB'] = np.clip(ref - nf, 0, None)
AF = np.clip(ref - naf, 0, None); tot = sum(S.values()); norm = np.where(tot > 0, AF / np.where(tot > 0, tot, 1), 0)
F = {k: np.nan_to_num(S[k] * norm / ref) for k in S}; F['Oth'] = np.nan_to_num(np.clip(naf, 0, None) / ref)
KEYS = ANTH + ['BB', 'Oth']

# bias ratio on the WRF lon-lat grid, then nearest-lon mapping to the mortality grid (same latitudes)
w = nc.Dataset(P + 'WRF-CHEM/PM2_5_TOTALS_ymonmean_lonlat_grid.nc'); M = np.ma.filled(w['PM2_5_DRY'][:], np.nan).astype(float) * 3.6
A = np.ma.filled(nc.Dataset(P + 'ACAG/PM25_NWSA_2019_regrid_all_V6.nc')['PM25'][:], np.nan).astype(float)[:len(M)]
A = np.where(A > 0, A, np.nan)
q_w = np.nanmean(M, 0) / np.nanmean(A, 0)
lonw = w['lon'][:]; lonm = nc.Dataset(B + 'mort_total.nc')['lon'][:]
j = np.abs(lonm[:, None] - lonw[None, :]).argmin(1); inside = (lonm >= lonw.min()) & (lonm <= lonw.max())
q = np.where(inside[None, :], q_w[:, j], np.nan); q = np.where(np.isfinite(q), q, 1.0)

def scen(low=None, high=False):
    G = {k: F[k].copy() for k in KEYS}
    lo = q < 1
    if low == 'nat':
        for k in KEYS[:-1]: G[k] = np.where(lo, F[k] * q, F[k])
    elif low == 'anth':
        Fa = sum(F[k] for k in ANTH); Fa_s = np.where(Fa > 0, Fa, 1)
        for k in KEYS:
            add = (1 - q) * F[k] / Fa_s if k in ANTH else 0
            G[k] = np.where(lo, F[k] * q + add, F[k])
    if high:
        hi = q > 1
        E = np.where(hi, 1 - 1 / q, 0)                       # excess as fraction of the model total
        bbcut = np.minimum(G['BB'], E); rest = E - bbcut
        G['BB'] = G['BB'] - bbcut
        others = [k for k in KEYS if k != 'BB']; so = sum(G[k] for k in others)
        for k in others: G[k] = G[k] - rest * G[k] / np.where(so > 0, so, 1)
        for k in KEYS: G[k] = np.where(hi, G[k] * q, G[k])   # renormalise to the ACAG total (M/A = q)
    if low == 'nat':
        G['Oth'] = 1 - sum(G[k] for k in KEYS[:-1])
    return G

k = (sid > 0)
def national(G): return {x: 100 * (D * G[x])[k].sum() / D[k].sum() for x in KEYS}
runs = {'Baseline': scen(), 'S1 natural': scen('nat'), 'S2 anthropogenic': scen('anth'), 'S3 fire excess': scen(high=True),
        'S1+S3': scen('nat', True), 'S2+S3': scen('anth', True)}
T = pd.DataFrame({n: national(G) for n, G in runs.items()}).T
T['sum'] = T[KEYS].sum(1)
pd.set_option('display.width', 200); print(T.round(2).to_string())
d = T.drop('Baseline') - T.loc['Baseline']
print('\nchange vs baseline (pp): min / max'); print(pd.DataFrame({'min': d[KEYS].min(), 'max': d[KEYS].max()}).round(2).T.to_string())
w_ = D[k]; qq = q[k]
print('\ndeaths in cells with q<0.8: %.1f%%, q>1.25: %.1f%%, deaths-weighted median q %.2f' % (100*w_[qq<0.8].sum()/w_.sum(), 100*w_[qq>1.25].sum()/w_.sum(), np.median(np.repeat(qq, 1))))
T.round(2).to_csv('sensitivity_national_shares.csv')

# ---- per-department diagnostics and department-level variant (bias ratio averaged over each department, deaths-weighted)
import sys as _s
from common import names
Mm = np.where(inside[None, :], np.nanmean(M, 0)[:, j], np.nan); Am = np.where(inside[None, :], np.nanmean(A, 0)[:, j], np.nan)
rows = []
qd = np.ones_like(q)
for i, n in sorted(names.items()):
    m = (sid == i) & np.isfinite(Mm) & np.isfinite(Am)
    if D[m].sum() == 0: continue
    qi = (D[m] * Mm[m]).sum() / (D[m] * Am[m]).sum(); qd[sid == i] = qi
    rows.append((n, round(D[sid == i].sum()), round(qi, 2)))
print(pd.DataFrame(rows, columns=['dept', 'deaths', 'q_dept']).to_string(index=False))
q = qd
runs = {'Baseline': scen(), 'S1 natural': scen('nat'), 'S2 anthropogenic': scen('anth'), 'S3 fire excess': scen(high=True),
        'S1+S3': scen('nat', True), 'S2+S3': scen('anth', True)}
T2 = pd.DataFrame({n: national(G) for n, G in runs.items()}).T
print('\nDEPARTMENT-LEVEL bias ratio'); print(T2.round(2).to_string())
T2.round(2).to_csv('sensitivity_national_shares_dept.csv')

# ---- best-guess regional scenario (reviewer's reasoning): east/Amazon deficit = natural; elsewhere deficit = anthropogenic; excess = fires first
EAST = {'AMAZONAS','CAQUETA','PUTUMAYO','GUAINIA','GUAVIARE','VAUPES','VICHADA','META','CASANARE','ARAUCA'}
east = np.isin(sid, [i for i, n in names.items() if n in EAST])
Gn = scen('nat', True); Ga = scen('anth', True)
Gr = {x: np.where(east, Gn[x], Ga[x]) for x in KEYS}
T2.loc['Regional (best guess)'] = pd.Series(national(Gr))
T2['sum'] = T2[KEYS].sum(axis=1)
print('\nwith regional scenario'); print(T2.round(2).to_string())
d = T2.drop('Baseline') - T2.loc['Baseline']
print(pd.DataFrame({'min': d[KEYS].min(), 'max': d[KEYS].max(), 'regional': d.loc['Regional (best guess)', KEYS]}).round(1).T.to_string())
T2.round(2).to_csv('sensitivity_national_shares_dept.csv')
# ranking check
for n, r in T2[KEYS].iterrows(): print(n, list(r.drop('Oth').sort_values(ascending=False).index[:4]))
