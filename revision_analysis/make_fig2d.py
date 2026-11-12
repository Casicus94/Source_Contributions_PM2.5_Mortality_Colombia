#!/usr/bin/env python3
"""
Rebuild Figure 2d (WRF-Chem vs ground stations, PM2.5, 2019) from a monthly model file.

Metrics per city, using the monthly means of 2019 (months present in both model and observations):
  MB   = mean(M - O)                         [ug m-3]
  NMB  = sum(M - O) / sum(O) * 100           [%]
  RMSE = sqrt(mean((M - O)^2))               [ug m-3]   (monthly, so RMSE >= |MB|)
Observations: hourly station file -> 0 and 9999 set to NaN -> monthly means (same as the original notebook).
Model: nearest grid cell to each city; optional scale factor, with optional exempt cities.

Example (ISTA):
  python make_fig2d.py --model Results/TOTALS/PM2_5_TOTALS_ymonmean.nc \
      --obs Datos_estaciones/Pollutant_trend_data/Hourmean/PM25_hourmean_Colombia.csv \
      --factor 3.6 --no-factor Bogota --out Figures/Errors_PM25_stations_revised
"""
import argparse
import numpy as np
import pandas as pd
import netCDF4 as nc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CITIES = {'Barbosa': (6.426635, -75.360970), 'Bogota': (4.649475, -74.108089), 'Bucaramanga': (7.114798, -73.117709),
          'Caldas': (6.085799, -75.635587), 'Cali': (3.419928, -76.510929), 'Copacabana': (6.343290, -75.517163),
          'Envigado': (6.160965, -75.548687), 'Itagui': (6.181957, -75.613805), 'Medellin': (6.241542, -75.584248)}
LABELS = {'Bogota': 'Bogotá', 'Itagui': 'Itagüí', 'Medellin': 'Medellín'}
BENCH = {'MB': 10.0, 'NMB': 30.0, 'RMSE': 15.0}   # same benchmark lines as the submitted figure


def model_series(path, var):
    d = nc.Dataset(path)
    pm = np.ma.filled(d[var][:].astype(float), np.nan)
    if "XLAT" in d.variables:
        la, lo = d["XLAT"][:], d["XLONG"][:]
        if la.ndim == 3:
            la, lo = la[0], lo[0]
        pick = lambda y, x: np.unravel_index(np.argmin((la - y) ** 2 + (lo - x) ** 2), la.shape)
    else:
        la, lo = d["lat"][:], d["lon"][:]
        pick = lambda y, x: (int(np.argmin(abs(la - y))), int(np.argmin(abs(lo - x))))
    tname = [k for k in d.variables if k.upper() in ("XTIME", "TIME")][0]
    t = nc.num2date(d[tname][:], d[tname].units, only_use_cftime_datetimes=False)
    months = [int(x.month) for x in t]
    out = {}
    for c, (y, x) in CITIES.items():
        i, j = pick(y, x)
        out[c] = pd.Series(pm[:, i, j], index=months)
    return pd.DataFrame(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True)
    ap.add_argument("--obs", required=True)
    ap.add_argument("--var", default="PM2_5_DRY")
    ap.add_argument("--factor", type=float, default=1.0, help="multiply model by this factor")
    ap.add_argument("--no-factor", nargs="*", default=[], help="cities NOT multiplied by --factor")
    ap.add_argument("--year", default="2019")
    ap.add_argument("--out", default="Errors_PM25_stations_revised")
    a = ap.parse_args()

    M = model_series(a.model, a.var)
    for c in M.columns:
        if c not in a.no_factor:
            M[c] = M[c] * a.factor
    o = pd.read_csv(a.obs, parse_dates=[0], index_col=0).loc[a.year].replace([0, 9999.0], np.nan)
    O = o.groupby(o.index.month).mean()[list(CITIES)]

    rows = []
    for c in CITIES:
        m, ob = M[c], O[c]
        idx = m.dropna().index.intersection(ob.dropna().index)
        m, ob = m.loc[idx], ob.loc[idx]
        rows.append({"city": c, "months": len(idx), "model_mean": m.mean(), "obs_mean": ob.mean(),
                     "MB": (m - ob).mean(), "NMB": 100 * (m - ob).sum() / ob.sum(),
                     "RMSE": np.sqrt(((m - ob) ** 2).mean())})
    df = pd.DataFrame(rows).set_index("city")
    df.round(3).to_csv(a.out + ".csv")
    print(df.round(1).to_string())
    print(f"mean MB {df.MB.mean():.1f}  mean NMB {df.NMB.mean():.1f}%  mean |NMB| {df.NMB.abs().mean():.1f}%  "
          f"mean RMSE {df.RMSE.mean():.1f}")

    fig, axs = plt.subplots(3, 1, figsize=(6.2, 7.2), sharex=True)
    x = np.arange(len(df))
    units = {"MB": "MB (µg m$^{-3}$)", "NMB": "NMB (%)", "RMSE": "RMSE (µg m$^{-3}$)"}
    for ax, k in zip(axs, ["MB", "NMB", "RMSE"]):
        ax.bar(x, df[k], color="teal", edgecolor="k", width=0.6)
        b = BENCH[k]
        ax.axhline(b, color="red", ls="--", lw=1)
        if k != "RMSE":
            ax.axhline(-b, color="red", ls="--", lw=1)
            ax.axhline(0, color="k", lw=0.6)
        ax.set_ylabel(units[k])
        ax.grid(axis="y", ls=":", lw=0.6)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    axs[-1].set_xticks(x)
    axs[-1].set_xticklabels([LABELS.get(c, c) for c in df.index], rotation=90)
    fig.tight_layout()
    fig.savefig(a.out + ".png", dpi=300)
    fig.savefig(a.out + ".pdf")
    print("wrote", a.out + ".csv/.png/.pdf")


if __name__ == "__main__":
    main()
