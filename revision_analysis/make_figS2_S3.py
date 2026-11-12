#!/usr/bin/env python3
"""
Figures S2 and S3 (revised): monthly evaluation of WRF-Chem PM2.5 (2019) against ACAG V6 (S2, maps) and ground
stations (S3, bars). Model multiplied by --factor (3.6 corrects the unit conversion in the model post-processing).

Metrics over the months available in both datasets:
  NMB  = sum(M-O)/sum(O)*100 ; RMSE = sqrt(mean((M-O)^2)) ; R2 = Pearson r^2 ;
  FAC2 = % of months with 0.5 <= M/O <= 2 ; IOA = 1 - sum((M-O)^2)/sum((|M-Obar|+|O-Obar|)^2) ; MB = mean(M-O)
Run from the MPIC_project folder or pass --root.
"""
import argparse, os
import numpy as np
import pandas as pd
import netCDF4 as nc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, BoundaryNorm
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


def metrics(M, O, axis=0):
    k = np.isfinite(M) & np.isfinite(O)
    M = np.where(k, M, np.nan); O = np.where(k, O, np.nan)
    n = np.sum(k, axis=axis)
    mb = np.nanmean(M - O, axis=axis)
    nmb = 100 * np.nansum(M - O, axis=axis) / np.nansum(O, axis=axis)
    rmse = np.sqrt(np.nanmean((M - O) ** 2, axis=axis))
    Mb = np.nanmean(M, axis=axis, keepdims=True); Ob = np.nanmean(O, axis=axis, keepdims=True)
    cov = np.nansum((M - Mb) * (O - Ob), axis=axis)
    r2 = (cov / np.sqrt(np.nansum((M - Mb) ** 2, axis=axis) * np.nansum((O - Ob) ** 2, axis=axis))) ** 2
    ratio = M / O
    fac2 = 100 * np.nansum((ratio >= 0.5) & (ratio <= 2), axis=axis) / n
    ioa = 1 - np.nansum((M - O) ** 2, axis=axis) / np.nansum((np.abs(M - Ob) + np.abs(O - Ob)) ** 2, axis=axis)
    return {"NMB": nmb, "RMSE": rmse, "R2": r2, "FAC2": fac2, "IOA": ioa, "MB": mb}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--model-grid", default="WRF-CHEM/PM2_5_TOTALS_ymonmean_lonlat_grid.nc")
    ap.add_argument("--model-native", default="WRF-CHEM/PM2_5_TOTALS_ymonmean.nc")
    ap.add_argument("--acag", default="ACAG/PM25_NWSA_2019_regrid_all_V6.nc")
    ap.add_argument("--obs", default="Datos_estaciones/Pollutant_trend_data/Hourmean/PM25_hourmean_Colombia.csv")
    ap.add_argument("--topo", default="GMRTv3_6_20190709topo.grd")
    ap.add_argument("--shp", default="World_Countries")
    ap.add_argument("--factor", type=float, default=3.6)
    ap.add_argument("--out-s2", default="Manuscript/Figures_SM/Errors_PM25_WRF_ACAG_V6_revised")
    ap.add_argument("--out-s3", default="Manuscript/Figures_SM/Errors_PM25_WRF_stations_revised")
    a = ap.parse_args()
    R = lambda p: os.path.join(a.root, p)

    # ------------------------------ S2: maps ------------------------------
    w = nc.Dataset(R(a.model_grid))
    M = np.ma.filled(w["PM2_5_DRY"][:].astype(float), np.nan) * a.factor
    lat, lon = w["lat"][:], w["lon"][:]
    A = np.ma.filled(nc.Dataset(R(a.acag))["PM25"][:].astype(float), np.nan)
    n = min(len(M), len(A)); M, A = M[:n], np.where(A[:n] > 0, A[:n], np.nan)
    with np.errstate(all="ignore"):
        maps = metrics(M, A)
    # land mask: ocean cells set to NaN in all maps
    box = shapely.box(-85, -6, -60, 15)
    land = unary_union([g.intersection(box) for g in shpreader.Reader(os.path.join(R(a.shp), "World_Countries.shp")).geometries() if g.intersects(box)])
    LON, LAT = np.meshgrid(np.asarray(lon), np.asarray(lat))
    ocean = ~shapely.contains_xy(land, LON, LAT)
    for key in maps:
        maps[key] = np.where(ocean, np.nan, maps[key])

    crs = ccrs.PlateCarree()
    v = nc.Dataset(R(a.topo))
    xt = np.arange(v["x_range"][0], v["x_range"][1] + 0.5 * v["spacing"][0], v["spacing"][0])
    yt = np.arange(v["y_range"][0], v["y_range"][1] + 0.5 * v["spacing"][1], v["spacing"][1])
    X, Y = np.meshgrid(xt, yt); Z = np.reshape(v["z"][:], X.shape)[::-1, :]
    s = (slice(None, None, 3), slice(None, None, 3))
    shps = [os.path.join(R(a.shp), f) for f in ("World_Countries.shp", "GDBIDEAM_departamentos.shp", "COL_no_dptos.shp")]
    lws = [0.1, 0.15, 0.3]
    bi = LinearSegmentedColormap.from_list("", ["royalblue", "lightblue", "white", "orange", "red"])
    seq = LinearSegmentedColormap.from_list("", ["lightblue", "royalblue", "orange", "red"])
    fac = LinearSegmentedColormap.from_list("", ["royalblue", "lightblue", "azure", "orange", "salmon"])
    lev_nmb = np.arange(-100, 251, 25)
    from matplotlib.colors import TwoSlopeNorm
    nmb_cmap = bi
    spec = [("NMB", "(a) NMB", lev_nmb, nmb_cmap, "%", "max", TwoSlopeNorm(vmin=-100, vcenter=0, vmax=250)),
            ("RMSE", "(b) RMSE", np.arange(0, 71, 5), seq, "µg m$^{-3}$", "max", None),
            ("R2", "(c) R$^2$", np.arange(0, 1.01, 0.1), seq, "", "neither", None),
            ("FAC2", "(d) FAC2", np.arange(0, 100.1, 5), fac, "%", "neither", None),
            ("IOA", "(e) IOA", np.arange(0, 1.01, 0.1), bi, "", "neither", None),
            ("MB", "(f) MB", np.arange(-25, 51, 5), bi, "µg m$^{-3}$", "both", TwoSlopeNorm(vmin=-25, vcenter=0, vmax=50))]
    lon2, lat2 = np.meshgrid(lon, lat)
    fig, axs = plt.subplots(2, 3, figsize=(14, 10.8), subplot_kw={"projection": crs})
    plt.subplots_adjust(hspace=0.12, wspace=0.12, top=0.88)
    for k, (key, title, lev, cmap, unit, ext, norm) in enumerate(spec):
        ax = axs.flat[k]
        im = ax.contourf(lon2, lat2, maps[key], levels=lev, cmap=cmap, norm=norm, extend=ext, transform=crs)
        ax.contour(X[s], Y[s], Z[s], levels=[1500, 3000], colors="grey", linewidths=0.6, transform=crs)
        for sh, lw in zip(shps, lws):
            ax.add_feature(ShapelyFeature(shpreader.Reader(sh).geometries(), crs, facecolor="none", edgecolor="k", linewidth=lw))
        ax.set_extent([-80, -66, -2, 13], crs=crs)
        for (c, (y, x)), mk in zip(CITIES.items(), MARKERS):
            ax.scatter(x, y, marker=mk, s=45, edgecolor="purple", facecolor="none", linewidth=1.2, transform=crs,
                       label=LABELS.get(c, c), zorder=5)
        gl = ax.gridlines(crs=crs, draw_labels=True, linewidth=0.6, color="k", alpha=0.5, linestyle="--")
        gl.xlocator = mticker.FixedLocator([-76, -70]); gl.ylocator = mticker.FixedLocator([2, 6, 10])
        gl.top_labels = gl.right_labels = False
        ax.set_title(title, loc="left", fontweight="bold")
        cb = fig.colorbar(im, ax=ax, orientation="horizontal", pad=0.07, shrink=0.95,
                          ticks=(np.arange(-100, 251, 50) if key == "NMB" else np.arange(-20, 51, 10)) if norm is not None else None, spacing="proportional")
        cb.set_label(unit)
        if key == "NMB":
            cb.ax.tick_params(labelsize=8)
    h, l = axs.flat[0].get_legend_handles_labels()
    fig.legend(h, l, frameon=False, ncol=5, loc="upper center", bbox_to_anchor=(0.5, 0.985), fontsize=11)
    fig.savefig(R(a.out_s2) + ".jpg", dpi=300, bbox_inches="tight"); fig.savefig(R(a.out_s2) + ".pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote", R(a.out_s2))

    # ------------------------------ S3: stations ------------------------------
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
        ms = pd.Series(Mn[:, i, j], index=months); ob = Om[c].reindex(ms.index)
        with np.errstate(all="ignore"):
            r = {k: float(v) for k, v in metrics(ms.values, ob.values).items()}
        r["city"] = c; rows.append(r)
    st = pd.DataFrame(rows).set_index("city")
    st.round(3).to_csv(R(a.out_s3) + ".csv"); print(st.round(2).to_string())
    bench = {"NMB": [30, -30], "RMSE": [15], "R2": [0.5], "FAC2": [50], "IOA": [0.5], "MB": [10, -10]}
    units = {"NMB": "%", "RMSE": "µg m$^{-3}$", "R2": "", "FAC2": "%", "IOA": "", "MB": "µg m$^{-3}$"}
    titles = {"NMB": "(a) NMB", "RMSE": "(b) RMSE", "R2": "(c) R$^2$", "FAC2": "(d) FAC2", "IOA": "(e) IOA", "MB": "(f) MB"}
    fig, axs = plt.subplots(2, 3, figsize=(15.5, 6.2), sharex=True)
    xs = np.arange(len(st))
    for ax, key in zip(axs.flat, ["NMB", "RMSE", "R2", "FAC2", "IOA", "MB"]):
        ax.bar(xs, st[key], color="lightsteelblue", edgecolor="k")
        for b in bench[key]:
            ax.axhline(b, color="red", ls=":", lw=1.3)
        if key in ("NMB", "MB"):
            ax.axhline(0, color="k", lw=0.6)
        if key in ("R2", "IOA"):
            ax.set_ylim(0, 1)
        if key == "FAC2":
            ax.set_ylim(0, 105)
        ax.set_ylabel(units[key]); ax.set_title(titles[key], loc="left", fontweight="bold")
        ax.set_xticks(xs); ax.set_xticklabels([LABELS.get(c, c) for c in st.index], rotation=90)
    fig.tight_layout()
    fig.savefig(R(a.out_s3) + ".pdf", bbox_inches="tight"); fig.savefig(R(a.out_s3) + ".jpg", dpi=300, bbox_inches="tight")
    print("wrote", R(a.out_s3))


if __name__ == "__main__":
    main()
