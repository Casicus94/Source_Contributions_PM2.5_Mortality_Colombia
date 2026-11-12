#!/usr/bin/env python3
"""
Figure 2 (revised): WRF-Chem PM2.5 evaluation for 2019. Maps are masked over the ocean.
(a) MB, (b) RMSE, (c) NMB of WRF-Chem against the ACAG V6 product (monthly means, conservative regrid to the WRF grid);
(d) MB, NMB and RMSE against ground stations (monthly means of hourly data).
MB = mean(M-O); NMB = sum(M-O)/sum(O)*100; RMSE = sqrt(mean((M-O)^2)), all over the months available in both datasets.
The model is multiplied by --factor (3.6 corrects a unit-conversion error in the post-processing of the model output).

Run from the MPIC_project folder (paths relative to it) or pass --root.
"""
import argparse, os
import numpy as np
import pandas as pd
import netCDF4 as nc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
from cartopy.feature import ShapelyFeature
import shapely
from shapely.ops import unary_union

np.bool = np.bool_
CITIES = {'Itagui': (6.181957, -75.613805), 'Caldas': (6.085799, -75.635587), 'Medellin': (6.241542, -75.584248),
          'Cali': (3.419928, -76.510929), 'Envigado': (6.160965, -75.548687), 'Bogota': (4.649475, -74.108089),
          'Bucaramanga': (7.114798, -73.117709), 'Barbosa': (6.426635, -75.360970), 'Copacabana': (6.343290, -75.517163)}
MARKERS = ['o', '*', 's', 'v', '^', '>', '<', 'p', 'P']
LABELS = {'Bogota': 'Bogotá', 'Itagui': 'Itagüí', 'Medellin': 'Medellín'}
BENCH = {'MB': 10.0, 'NMB': 30.0, 'RMSE': 15.0}


def metrics(M, O, axis=0):
    mb = np.nanmean(M - O, axis=axis)
    nmb = 100 * np.nansum(M - O, axis=axis) / np.nansum(O, axis=axis)
    rmse = np.sqrt(np.nanmean((M - O) ** 2, axis=axis))
    return mb, nmb, rmse


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--model-grid", default="WRF-CHEM/PM2_5_TOTALS_ymonmean_lonlat_grid.nc", help="WRF on the ACAG lon-lat grid")
    ap.add_argument("--model-native", default="WRF-CHEM/PM2_5_TOTALS_ymonmean.nc", help="WRF on its native grid (XLAT/XLONG)")
    ap.add_argument("--acag", default="ACAG/PM25_NWSA_2019_regrid_all_V6.nc")
    ap.add_argument("--obs", default="Datos_estaciones/Pollutant_trend_data/Hourmean/PM25_hourmean_Colombia.csv")
    ap.add_argument("--topo", default="GMRTv3_6_20190709topo.grd")
    ap.add_argument("--shp", default="World_Countries")
    ap.add_argument("--factor", type=float, default=3.6)
    ap.add_argument("--out", default="Figures/Errors_yearly_V6_only_pm_revised")
    a = ap.parse_args()
    R = lambda p: os.path.join(a.root, p)

    # ---------------- maps vs ACAG ----------------
    w = nc.Dataset(R(a.model_grid))
    M = np.ma.filled(w["PM2_5_DRY"][:].astype(float), np.nan) * a.factor
    lat, lon = w["lat"][:], w["lon"][:]
    A = np.ma.filled(nc.Dataset(R(a.acag))["PM25"][:].astype(float), np.nan)
    n = min(len(M), len(A)); M, A = M[:n], A[:n]
    A = np.where(A > 0, A, np.nan)
    mb_map, nmb_map, rmse_map = metrics(M, A)
    # land mask (ocean cells set to NaN in all three maps)
    box = shapely.box(-85, -6, -60, 15)
    land = unary_union([g.intersection(box) for g in shpreader.Reader(os.path.join(R(a.shp), "World_Countries.shp")).geometries() if g.intersects(box)])
    LON, LAT = np.meshgrid(np.asarray(lon), np.asarray(lat))
    ocean = ~shapely.contains_xy(land, LON, LAT)
    for f in (mb_map, nmb_map, rmse_map):
        f[ocean] = np.nan

    # ---------------- stations ----------------
    d = nc.Dataset(R(a.model_native))
    Mn = np.ma.filled(d["PM2_5_DRY"][:].astype(float), np.nan) * a.factor
    la, lo = d["XLAT"][:], d["XLONG"][:]
    t = nc.num2date(d["XTIME"][:], d["XTIME"].units, only_use_cftime_datetimes=False)
    months = [x.month for x in t]
    o = pd.read_csv(R(a.obs), parse_dates=[0], index_col=0).loc["2019"].replace([0, 9999.0], np.nan)
    Om = o.groupby(o.index.month).mean()
    rows = []
    for c, (y, x) in CITIES.items():
        i, j = np.unravel_index(np.argmin((la - y) ** 2 + (lo - x) ** 2), la.shape)
        ms = pd.Series(Mn[:, i, j], index=months)
        ob = Om[c].reindex(ms.index)
        k = ms.notna() & ob.notna()
        mb, nmb, rmse = metrics(ms[k].values, ob[k].values)
        rows.append({"city": c, "months": int(k.sum()), "model": ms[k].mean(), "obs": ob[k].mean(), "MB": mb, "NMB": nmb, "RMSE": rmse})
    st = pd.DataFrame(rows).set_index("city").sort_index()
    st.round(3).to_csv(R(a.out) + "_stations.csv")
    print(st.round(1).to_string())
    print(f"station mean: MB {st.MB.mean():.1f}  NMB {st.NMB.mean():.1f}%  |NMB| {st.NMB.abs().mean():.1f}%  RMSE {st.RMSE.mean():.1f}")

    # ---------------- figure ----------------
    crs = ccrs.PlateCarree()
    v = nc.Dataset(R(a.topo))
    xt = np.arange(v["x_range"][0], v["x_range"][1] + 0.5 * v["spacing"][0], v["spacing"][0])
    yt = np.arange(v["y_range"][0], v["y_range"][1] + 0.5 * v["spacing"][1], v["spacing"][1])
    X, Y = np.meshgrid(xt, yt); Z = np.reshape(v["z"][:], X.shape)[::-1, :]
    s = (slice(None, None, 3), slice(None, None, 3))
    shps = [os.path.join(R(a.shp), f) for f in ("World_Countries.shp", "GDBIDEAM_departamentos.shp", "COL_no_dptos.shp")]
    lws = [0.1, 0.15, 0.3]
    bi = LinearSegmentedColormap.from_list("", ["royalblue", "lightblue", "white", "orange", "red"])
    seq = LinearSegmentedColormap.from_list("", ["white", "lightyellow", "orange", "red", "darkred"])
    panels = [(mb_map, "(a) MB", np.arange(-25, 51, 5), bi, "µg m$^{-3}$", "both"),
              (rmse_map, "(b) RMSE", np.arange(0, 71, 5), seq, "µg m$^{-3}$", "max"),
              (nmb_map, "(c) NMB", np.arange(-100, 251, 25), bi, "%", "both")]

    fig = plt.figure(figsize=(9, 11))
    gs = GridSpec(8, 2, left=0.09, right=0.95, hspace=0.25, wspace=0.2, top=0.9, bottom=0.08,
                  height_ratios=[1, 0.05, 0.1, 0.32, 0.33, 0.33, 0.01, 0.05])
    pos = [(gs[0, 0], gs[1, 0]), (gs[0, 1], gs[1, 1]), (gs[3:6, 0], gs[7, 0])]
    lon2, lat2 = np.meshgrid(lon, lat)
    for k, (fld, title, lev, cmap, unit, ext) in enumerate(panels):
        ax = plt.subplot(pos[k][0], projection=crs)
        nrm = TwoSlopeNorm(vmin=lev[0], vcenter=0, vmax=lev[-1]) if lev[0] < 0 else None
        im = ax.contourf(lon2, lat2, fld, levels=lev, cmap=cmap, norm=nrm, extend=ext, transform=crs)
        ax.contour(X[s], Y[s], Z[s], levels=[1500, 3000], colors="grey", linewidths=0.6, transform=crs)
        for sh, lw in zip(shps, lws):
            ax.add_feature(ShapelyFeature(shpreader.Reader(sh).geometries(), crs, facecolor="none", edgecolor="k", linewidth=lw))
        ax.set_extent([-80, -66, -2, 13], crs=crs)
        for (c, (y, x)), mk in zip(CITIES.items(), MARKERS):
            ax.scatter(x, y, marker=mk, s=45, edgecolor="purple", facecolor="none", linewidth=1.2, transform=crs,
                       label=LABELS.get(c, c), zorder=5)
        if k == 0:
            ax.legend(frameon=False, ncol=5, bbox_to_anchor=(1.1, 1.38), loc="upper center", fontsize=10)
        gl = ax.gridlines(crs=crs, draw_labels=True, linewidth=0.6, color="k", alpha=0.5, linestyle="--")
        gl.xlocator = mticker.FixedLocator([-76, -70]); gl.ylocator = mticker.FixedLocator([2, 6, 10])
        gl.top_labels = gl.right_labels = False
        gl.xlabel_style = gl.ylabel_style = {"size": 9}
        ax.set_title(title, loc="left", fontweight="bold", fontsize=10)
        cax = plt.subplot(pos[k][1])
        cb = plt.colorbar(im, cax=cax, orientation="horizontal"); cb.set_label(unit, fontsize=9); cb.ax.tick_params(labelsize=8)

    sub = GridSpec(3, 1, left=0.57, right=0.95, top=0.49, bottom=0.12, hspace=0.12)
    xs = np.arange(len(st))
    for k, key in enumerate(["MB", "NMB", "RMSE"]):
        ax = fig.add_subplot(sub[k])
        ax.bar(xs, st[key], color="lightsteelblue", edgecolor="k", width=0.6)
        ax.axhline(BENCH[key], color="red", ls="--", lw=1)
        if key != "RMSE":
            ax.axhline(-BENCH[key], color="red", ls="--", lw=1); ax.axhline(0, color="k", lw=0.6)
        ax.set_ylabel({"MB": "MB (µg m$^{-3}$)", "NMB": "NMB (%)", "RMSE": "RMSE (µg m$^{-3}$)"}[key], fontsize=9)
        ax.tick_params(labelsize=8); ax.grid(axis="y", ls=":", lw=0.5)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.set_xticks(xs)
        ax.set_xticklabels([LABELS.get(c, c) for c in st.index] if k == 2 else [], rotation=90, fontsize=9)
        if k == 0:
            ax.set_title("(d)", loc="left", fontweight="bold", fontsize=10)
    fig.savefig(R(a.out) + ".jpg", dpi=300, bbox_inches="tight")
    fig.savefig(R(a.out) + ".pdf", dpi=300, bbox_inches="tight")
    np.savez(R(a.out) + "_maps.npz", mb=mb_map, nmb=nmb_map, rmse=rmse_map, lat=lat, lon=lon)
    print("wrote", R(a.out) + ".jpg/.pdf")


if __name__ == "__main__":
    main()
