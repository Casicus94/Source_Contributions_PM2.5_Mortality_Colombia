#!/usr/bin/env python3
"""
Is double counting between EDGAR agricultural waste burning (AWB) and FINN fires small?
(Reviewer 3, major comment 5.)

For Colombia, the model domain and two groups of departments, this script compares 2019 annual emitted mass of
PM2.5, OC, BC and NH3 from
  * EDGAR v8.1 FT2022 AWB                    (agricultural waste burning, inside our "Agriculture" group)
  * EDGAR v8.1 FT2022 AGS + ENF + MNM        (the rest of the "Agriculture" group)
  * FINN v1.5 daily fire emissions           (our "Biomass burning" sector, all fires switched off together)
and writes a CSV, a printed table and a bar figure.

Assumptions you should check (edit the flags if they differ):
  * EDGAR files:  <edgar-root>/<species>/<sector>/v8.1_FT2022_AP_<species>_2019_<sector>_emi.nc
                  variable 'emissions', tonnes per cell per year, coordinates 'lat' and 'lon'.
  * FINN file(s): FINN v1.5 daily text file, one row per fire pixel, header line with columns such as
                  DAY, LATI, LONGI, PM25, OC, BC, NH3 ...   (column names are matched case-insensitively)
                  aerosols in kg/day, gases in mol/day (FINN v1.5 convention). Use --aerosol-unit / --gas-unit to change.
Usage (on the ISTA cluster):
  python AWB_vs_FINN_check.py --finn-file /path/to/FINNv1.5_2019.*.txt \
      --edgar-root /nfs/scistore16/mullegrp/acasalla/MPI-C/Emissions \
      --shp-dir /nfs/scistore16/mullegrp/acasalla/MPI-C/World_Countries --out awb_check
"""
import argparse, glob, os, sys, unicodedata
import numpy as np
import pandas as pd
import netCDF4 as nc
import shapefile                                   # pyshp (installed with cartopy)
from shapely.geometry import box, shape
from shapely.ops import unary_union
try:
    from shapely import contains_xy as _contains_xy          # shapely >= 2
except ImportError:                                            # shapely 1.x
    from shapely.vectorized import contains as _contains_xy

SPECIES = ["PM2.5", "OC", "BC", "NH3"]
EDGAR_AWB = ["AWB"]
EDGAR_REST = ["AGS", "ENF", "MNM"]
FINN_NAMES = {"PM2.5": ["PM25", "PM2.5", "PM25_TOT"], "OC": ["OC"], "BC": ["BC"], "NH3": ["NH3"]}
LAT_NAMES, LON_NAMES = ["LATI", "LAT", "LATITUDE"], ["LONGI", "LON", "LONGITUDE"]
NH3_MOLAR_MASS = 17.031  # g mol-1
DOMAIN_D02 = box(-80.0, -2.0, -66.0, 13.0)   # inner WRF domain (lon_min, lat_min, lon_max, lat_max)
FIRE_DEPTS = ["ARAUCA", "CASANARE", "PUTUMAYO", "HUILA"]              # biomass burning is the leading source
AGRI_DEPTS = ["CALDAS", "QUINDIO", "RISARALDA", "TOLIMA"]              # agriculture is the leading source


def norm(s):
    return "".join(c for c in unicodedata.normalize("NFD", str(s).upper()) if unicodedata.category(c) != "Mn").strip()


def read_shp(path):
    sf = shapefile.Reader(path, encoding="latin-1")
    fields = [f[0] for f in sf.fields[1:]]
    return [(dict(zip(fields, sr.record)), shape(sr.shape.__geo_interface__)) for sr in sf.shapeRecords()]


def load_regions(shp_dir):
    col = [g for _, g in read_shp(os.path.join(shp_dir, "COL_no_dptos.shp"))]
    dep_geoms = {norm(a["NMG"]): g for a, g in read_shp(os.path.join(shp_dir, "GDBIDEAM_departamentos.shp"))}
    regions = {"Colombia": unary_union(col), "Model domain D02": DOMAIN_D02}
    for label, names in (("Fire-dominated departments", FIRE_DEPTS), ("Agriculture-dominated departments", AGRI_DEPTS)):
        missing = [n for n in names if n not in dep_geoms]
        if missing:
            print(f"WARNING: departments not found in shapefile: {missing}")
        regions[label] = unary_union([dep_geoms[n] for n in names if n in dep_geoms])
    return regions


def inside(geom, lon, lat):
    return np.asarray(_contains_xy(geom, np.asarray(lon, float), np.asarray(lat, float)))


def edgar_total(root, species, sector, geoms, year=2019):
    f = os.path.join(root, species, sector, f"v8.1_FT2022_AP_{species}_{year}_{sector}_emi.nc")
    if not os.path.exists(f):
        print(f"  missing EDGAR file (counted as 0): {f}")
        return {k: 0.0 for k in geoms}
    with nc.Dataset(f) as d:
        var = d.variables["emissions"]
        z = np.ma.filled(var[:], 0.0).astype(float)
        z = z.reshape(-1, z.shape[-2], z.shape[-1]).sum(axis=0)       # drop/sum any time dimension
        lat, lon = d.variables["lat"][:], d.variables["lon"][:]
    lon2, lat2 = np.meshgrid(lon, lat)
    return {k: float(z[inside(g, lon2, lat2)].sum()) for k, g in geoms.items()}     # tonnes per year


def finn_totals(files, geoms, aerosol_unit, gas_unit, sep=None, chunk=2_000_000):
    """Sum FINN emissions (tonnes) per region, for the year and for Jan-Apr (DAY<=120)."""
    box_ = DOMAIN_D02.bounds
    tot = {k: {s: [0.0, 0.0] for s in SPECIES} for k in geoms}
    for f in files:
        with open(f) as fh:
            head = fh.readline()
        sep_ = sep or ("," if head.count(",") >= head.count(" ") else r"\s+")
        cols = [c.strip().upper() for c in head.strip().split("," if sep_ == "," else None)]
        pick = {}
        for s in SPECIES:
            for cand in FINN_NAMES[s]:
                if cand in cols:
                    pick[s] = cand
                    break
        lat_c = next(c for c in LAT_NAMES if c in cols)
        lon_c = next(c for c in LON_NAMES if c in cols)
        day_c = "DAY" if "DAY" in cols else None
        print(f"FINN {os.path.basename(f)}: using lat={lat_c}, lon={lon_c}, day={day_c}, species={pick}")
        use = [lat_c, lon_c] + list(pick.values()) + ([day_c] if day_c else [])
        reader = pd.read_csv(f, sep=sep_, usecols=lambda c: c.strip().upper() in use, chunksize=chunk, index_col=False,  # FINN rows end with a trailing comma
                             engine="c" if sep_ == "," else "python", skipinitialspace=True)
        for df in reader:
            df.columns = [c.strip().upper() for c in df.columns]
            df = df[(df[lon_c] >= box_[0]) & (df[lon_c] <= box_[2]) & (df[lat_c] >= box_[1]) & (df[lat_c] <= box_[3])]
            if df.empty:
                continue
            early = (df[day_c] <= 120).to_numpy() if day_c else np.zeros(len(df), bool)
            for k, g in geoms.items():
                m = inside(g, df[lon_c].to_numpy(), df[lat_c].to_numpy())
                for s, c in pick.items():
                    v = df[c].to_numpy(float)
                    if s == "NH3":     # mol/day -> tonnes
                        v = v * NH3_MOLAR_MASS / 1e6 if gas_unit == "mol/day" else v / 1e3
                    else:              # kg/day -> tonnes
                        v = v / 1e3 if aerosol_unit == "kg/day" else v * 1.0
                    tot[k][s][0] += float(v[m].sum())
                    tot[k][s][1] += float(v[m & early].sum())
    return tot


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--finn-file", nargs="+", required=True, help="FINN v1.5 daily text file(s) for 2019 (globs allowed)")
    ap.add_argument("--edgar-root", required=True)
    ap.add_argument("--shp-dir", required=True, help="folder with COL_no_dptos.shp and GDBIDEAM_departamentos.shp")
    ap.add_argument("--out", default="awb_check")
    ap.add_argument("--aerosol-unit", default="kg/day", choices=["kg/day", "t/day"])
    ap.add_argument("--gas-unit", default="mol/day", choices=["mol/day", "kg/day"])
    a = ap.parse_args()
    files = sorted(sum([glob.glob(p) for p in a.finn_file], []))
    if not files:
        sys.exit("no FINN file found")
    geoms = load_regions(a.shp_dir)
    rows = []
    finn = finn_totals(files, geoms, a.aerosol_unit, a.gas_unit)
    for s in SPECIES:
        awb = {k: 0.0 for k in geoms}
        rest = {k: 0.0 for k in geoms}
        for sec in EDGAR_AWB:
            for k, v in edgar_total(a.edgar_root, s, sec, geoms).items():
                awb[k] += v
        for sec in EDGAR_REST:
            for k, v in edgar_total(a.edgar_root, s, sec, geoms).items():
                rest[k] += v
        for k in geoms:
            fy, fe = finn[k][s]
            rows.append({"region": k, "species": s, "EDGAR_AWB_t": awb[k], "EDGAR_AGS+ENF+MNM_t": rest[k],
                         "FINN_year_t": fy, "FINN_JanApr_t": fe,
                         "AWB_over_FINN_%": 100 * awb[k] / fy if fy else np.nan,
                         "AWB_over_FINN_JanApr_%": 100 * awb[k] / fe if fe else np.nan,
                         "AWB_share_of_agri_%": 100 * awb[k] / (awb[k] + rest[k]) if (awb[k] + rest[k]) else np.nan})
    df = pd.DataFrame(rows)
    df.to_csv(a.out + ".csv", index=False)
    pd.set_option("display.width", 220)
    print(df.round(2).to_string(index=False))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    d = df[df.region == "Colombia"].set_index("species").loc[SPECIES]
    x = np.arange(len(SPECIES)); w = 0.27
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(x - w, d["EDGAR_AWB_t"], w, label="EDGAR AWB", color="#4c78a8")
    ax.bar(x, d["EDGAR_AGS+ENF+MNM_t"], w, label="EDGAR AGS+ENF+MNM", color="#54a24b")
    ax.bar(x + w, d["FINN_year_t"], w, label="FINN (all fires)", color="#e08a2c")
    ax.set_yscale("log"); ax.set_xticks(x); ax.set_xticklabels([s.replace("PM2.5", "PM$_{2.5}$").replace("NH3", "NH$_3$") for s in SPECIES])
    ax.set_ylabel("2019 emissions over Colombia (t yr$^{-1}$)"); ax.legend(frameon=False)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout(); fig.savefig(a.out + ".png", dpi=300)
    print(f"wrote {a.out}.csv and {a.out}.png")


if __name__ == "__main__":
    main()
