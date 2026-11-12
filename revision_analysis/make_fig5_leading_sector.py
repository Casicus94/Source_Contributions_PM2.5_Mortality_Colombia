#!/usr/bin/env python3
"""Figure 5 (revised): leading sector of PM2.5-attributable mortality per department (excluding 'Other natural sources'),
from the corrected sector fractions (Table 1)."""
import sys, unicodedata
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import cartopy.crs as ccrs, cartopy.feature as cf, shapefile
from shapely.geometry import shape
csv, shp_dir, out = sys.argv[1:4]
norm = lambda s: "".join(c for c in unicodedata.normalize("NFD", str(s).upper()) if unicodedata.category(c) != "Mn").strip()
T = pd.read_csv(csv).dropna(subset=["Agr"]); T["key"] = T.dept.map(norm)
cols = ["Agr", "Avi", "Ene", "Ind", "BB", "ORT", "RT", "Shp", "Wst"]
T["lead"] = T[cols].astype(float).idxmax(axis=1); lead = dict(zip(T.key, T.lead)); lead["VALLE"] = lead["VALLE DEL CAUCA"]
COL = {"Agr": ("darkcyan", "Agriculture"), "Ene": ("lightsteelblue", "Energy and residential"), "BB": (matplotlib.colors.to_rgba("indianred", 0.8), "Biomass burning"),
       "Ind": ("#386cb0", "Industry"), "RT": ("#f0027f", "Road transport")}
crs = ccrs.PlateCarree()
fig, ax = plt.subplots(figsize=(7, 8), subplot_kw={"projection": crs})
sf = shapefile.Reader(f"{shp_dir}/World_Countries.shp", encoding="latin-1")
for s in sf.shapes():
    ax.add_geometries([shape(s.__geo_interface__)], crs, facecolor="#efeedc", edgecolor="k", linewidth=0.6)
used = set()
sf = shapefile.Reader(f"{shp_dir}/GDBIDEAM_departamentos.shp", encoding="latin-1")
for r in sf.shapeRecords():
    k = lead[norm(r.record["NMG"])]; used.add(k)
    ax.add_geometries([shape(r.shape.__geo_interface__)], crs, facecolor=COL[k][0], edgecolor="k", linewidth=0.3)
sf = shapefile.Reader(f"{shp_dir}/COL_no_dptos.shp", encoding="latin-1")
for s in sf.shapes():
    ax.add_geometries([shape(s.__geo_interface__)], crs, facecolor="none", edgecolor="k", linewidth=1.2)
ax.set_extent([-80, -67, -2, 12], crs=crs)
gl = ax.gridlines(crs=crs, draw_labels=True, linewidth=0.4, color="grey", alpha=0.5, linestyle="--")
gl.top_labels = gl.right_labels = False
ax.legend(handles=[Patch(facecolor=COL[k][0], edgecolor="k", label=COL[k][1]) for k in ["Agr", "Ene", "BB"] if k in used],
          loc="upper center", bbox_to_anchor=(0.5, -0.06), ncol=3, frameon=False, fontsize=11)
fig.savefig(out, dpi=300, bbox_inches="tight"); print("wrote", out, sorted(used))
