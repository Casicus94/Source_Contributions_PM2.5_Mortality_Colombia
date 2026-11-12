#!/usr/bin/env python3
"""Reviewer 3, Major comment 5: AWB (EDGAR) vs FINN double counting, death-weighted, and the effect of masking.
Also: spatial correlation of the simulated agricultural contribution with AWB and AGS+MNM emissions.
Run from the MPIC_project folder. Needs Emissions/{PM2.5}/{AWB,AGS,MNM}/v8.1_FT2022_*.nc, Fires/FINN_2019_D02_domain.txt,
World_Countries/{COL_no_dptos,GDBIDEAM_departamentos}.shp, Table1_fixed_fractions.csv (Table 1 values),
Agriculture/PM2_5_Agriculture_ymean.nc and TOTALS_No_FIRES/PM2_5_TOTALS_No_FIRES_ymean.nc.
FINN GENVEG = 9 are cropland fires."""
import sys, numpy as np, pandas as pd, netCDF4 as nc, shapefile, unicodedata, shapely
from shapely.geometry import shape
from scipy.ndimage import uniform_filter
T1 = sys.argv[1] if len(sys.argv) > 1 else "Table1_fixed_fractions.csv"
norm = lambda s: "".join(c for c in unicodedata.normalize("NFD", str(s).upper()) if unicodedata.category(c) != "Mn").strip()

def edgar(sp, sec):
    d = nc.Dataset(f"Emissions/{sp}/{sec}/v8.1_FT2022_AP_{sp}_2019_{sec}_emi.nc")
    k = [v for v in d.variables if v not in ("lat", "lon", "time")][0]
    z = np.ma.filled(d[k][:], 0).astype(float); z = z.reshape(-1, z.shape[-2], z.shape[-1]).sum(0)
    return np.asarray(d["lat"][:]), np.asarray(d["lon"][:]), z

la, lo, awb = edgar("PM2.5", "AWB"); _, _, ags = edgar("PM2.5", "AGS"); _, _, mnm = edgar("PM2.5", "MNM")
rest = ags + mnm
f = pd.read_csv("Fires/FINN_2019_D02_domain.txt", skipinitialspace=True, index_col=False); f.columns = [c.strip().upper() for c in f.columns]
ii = np.round((f.LATI - la[0]) / (la[1] - la[0])).astype(int); jj = np.round((f.LONGI - lo[0]) / (lo[1] - lo[0])).astype(int)
ok = (ii >= 0) & (ii < len(la)) & (jj >= 0) & (jj < len(lo))
def grid(m):
    g = np.zeros_like(awb); np.add.at(g, (ii[ok & m], jj[ok & m]), f.PM25[ok & m] / 1e3); return g   # kg/day -> t
finn_all, finn_crop = grid(np.ones(len(f), bool)), grid(f.GENVEG == 9)
LO, LA = np.meshgrid(lo, la)
col = shapely.unary_union([shape(s.__geo_interface__) for s in shapefile.Reader("World_Countries/COL_no_dptos.shp", encoding="latin-1").shapes()])
inC = shapely.contains_xy(col, LO, LA)
print(f"Colombia (t): AWB {awb[inC].sum():.0f}, AGS+MNM {rest[inC].sum():.0f}, FINN {finn_all[inC].sum():.0f}, FINN cropland {finn_crop[inC].sum():.0f}")
print(f"overlap min(AWB,FINN cropland) {np.minimum(awb, finn_crop)[inC].sum():.0f} t; min(AWB,FINN all) {np.minimum(awb, finn_all)[inC].sum():.0f} t; "
      f"AWB in cells with any FINN fire (masking) {awb[inC & (finn_all > 0)].sum():.0f} t")

# department fractions of local agricultural primary PM2.5, weighted by agricultural deaths (Table 1)
sf = shapefile.Reader("World_Countries/GDBIDEAM_departamentos.shp", encoding="latin-1"); fl = [x[0] for x in sf.fields[1:]]
fr = {}
for rec, shp in zip(sf.records(), sf.shapes()):
    n = norm(dict(zip(fl, rec))["NMG"]); m = shapely.contains_xy(shape(shp.__geo_interface__), LO, LA); agri = (awb + rest)[m].sum()
    fr[n] = {k: (v[m].sum() / agri if agri > 0 else 0) for k, v in
             {"crop": np.minimum(awb, finn_crop), "all": np.minimum(awb, finn_all), "mask": np.where(finn_all > 0, awb, 0)}.items()}
fr["VALLE DEL CAUCA"] = fr.pop("VALLE"); fr["BOGOTA"] = fr["CUNDINAMARCA"]          # Bogota D.C. lies within the Cundinamarca polygon
T = pd.read_csv(T1).dropna(subset=["Agr"]); T = T[T.dept != "COLOMBIA"]; T["key"] = T.dept.map(norm)
T["agr_deaths"] = T.Agr / 100 * T.Tot; tot, ad = T.Tot.sum(), T.agr_deaths.sum()
for k in ("crop", "all", "mask"):
    lost = (T.agr_deaths * T.key.map(lambda d: fr[d][k])).sum()
    print(f"{k:5s}: {100*lost/ad:5.1f}% of the agricultural contribution; national agricultural share {100*ad/tot:.2f}% -> {100*(ad-lost)/tot:.2f}%")

# spatial correlation of the simulated agricultural contribution with the emission fields
A = np.squeeze(nc.Dataset("Agriculture/PM2_5_Agriculture_ymean.nc")["PM2_5_DRY"][:]).astype(float)
F = np.squeeze(nc.Dataset("TOTALS_No_FIRES/PM2_5_TOTALS_No_FIRES_ymean.nc")["PM2_5_DRY"][:]).astype(float)
g = nc.Dataset("Agriculture/PM2_5_Agriculture_ymean.nc"); lat, lon = np.asarray(g["XLAT"][:]), np.asarray(g["XLONG"][:])
a_i = np.clip(np.searchsorted(la, lat.ravel()) - 1, 0, len(la) - 1); a_j = np.clip(np.searchsorted(lo, lon.ravel()) - 1, 0, len(lo) - 1)
s = np.zeros_like(awb); n = np.zeros_like(awb); np.add.at(s, (a_i, a_j), (F - A).ravel()); np.add.at(n, (a_i, a_j), 1)
m = n > 0; C = s[m] / n[m]
for name, e in (("AWB", awb), ("AGS+MNM", rest)):
    print(f"r(agricultural contribution, {name} emissions, 9x9 smoothing) = {np.corrcoef(uniform_filter(e, 9)[m], C)[0, 1]:.2f}")
