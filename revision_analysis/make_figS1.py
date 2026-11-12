#!/usr/bin/env python3
"""Figure S1: ACAG V6 annual mean PM2.5 (2019) vs ground stations.
MB = A - O, NMB = 100 (A - O) / O, RMSE of the annual means (= |MB|), with
O the annual mean of the monthly station means and A the ACAG annual mean of the
same months at the grid cell (0.01 deg) nearest to each city."""
import numpy as np, pandas as pd, netCDF4 as nc, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import sys
ROOT = (sys.argv[1] if len(sys.argv) > 1 else ".").rstrip("/") + "/"  # path to the MPIC_project folder
CITIES = {'Itagui': (6.181957, -75.613805), 'Caldas': (6.085799, -75.635587), 'Medellin': (6.241542, -75.584248),
          'Cali': (3.419928, -76.510929), 'Envigado': (6.160965, -75.548687), 'Bogota': (4.649475, -74.108089),
          'Bucaramanga': (7.114798, -73.117709), 'Barbosa': (6.426635, -75.360970), 'Copacabana': (6.343290, -75.517163)}
LAB = {'Itagui': 'Itagüí', 'Medellin': 'Medellín', 'Bogota': 'Bogotá'}
A = {c: [] for c in CITIES}
for m in range(1, 13):
    d = nc.Dataset(f"{ROOT}ACAG/PM25_NWSA_2019{m:02d}_V6.nc"); la, lo = d['lat'][:], d['lon'][:]
    v = np.squeeze(d['PM25'][:])
    for c, (y, x) in CITIES.items():
        A[c].append(float(v[np.argmin(abs(la - y)), np.argmin(abs(lo - x))]))
A = pd.DataFrame(A, index=range(1, 13))
o = pd.read_csv(f"{ROOT}Datos_estaciones/Pollutant_trend_data/Hourmean/PM25_hourmean_Colombia.csv",
                parse_dates=[0], index_col=0).loc['2019'].replace([0, 9999.0], np.nan)
O = o.groupby(o.index.month).mean()[list(CITIES)]
rows = []
for c in CITIES:
    idx = O[c].dropna().index
    a, ob = A[c].loc[idx].mean(), O[c].loc[idx].mean()
    rows.append({'city': c, 'obs': ob, 'acag': a, 'MB': a - ob, 'NMB': 100 * (a - ob) / ob, 'RMSE': abs(a - ob)})
df = pd.DataFrame(rows).set_index('city'); df.round(2).to_csv('FigS1_ACAG_stations.csv'); print(df.round(1))
print('mean MB %.1f  mean NMB %.1f' % (df.MB.mean(), df.NMB.mean()))
fig = plt.figure(figsize=(8, 8)); gs = GridSpec(3, 1, left=0.09, right=0.95, hspace=0.2, top=0.9, bottom=0.08)
x = np.arange(1, len(df) + 1)
spec = [('MB', 'PM$_{2.5}$ MB (μg m$^{-3}$)', [10, -10]), ('NMB', 'NMB (%)', [30, -30]), ('RMSE', 'PM$_{2.5}$ RMSE (μg m$^{-3}$)', [15])]
for k, (key, yl, bench) in enumerate(spec):
    ax = plt.subplot(gs[k])
    ax.bar(x, df[key], color='darkcyan', edgecolor='k', label='PM$_{2.5}$')
    for b in bench: ax.axhline(b, color='red', linestyle='--', linewidth=1)
    if key != 'RMSE': ax.axhline(0, linestyle=':', color='k')
    ax.set_xticks(x); ax.grid(linestyle=':', axis='x'); ax.set_ylabel(yl, fontsize=9); ax.tick_params(labelsize=9)
    ax.set_title('(%s)' % 'abc'[k], fontweight='bold', fontsize=9, loc='left')
    if k < 2: ax.set_xticklabels([])
    else: ax.set_xticklabels([LAB.get(c, c) for c in df.index], rotation=90, fontsize=9)
    if k == 0: ax.legend(frameon=False, bbox_to_anchor=(0.5, 1.2), loc='upper center', fontsize=9)
plt.savefig('Errors_yearly_ACAG_stations_V6.pdf', bbox_inches='tight', dpi=300)
plt.savefig('Errors_yearly_ACAG_stations_V6.jpg', bbox_inches='tight', dpi=300)
