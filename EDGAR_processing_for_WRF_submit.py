
from __future__ import annotations

import math
import re
import sys
import warnings
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
import requests
import tqdm
import xarray as xr

# ─────────────────── user settings ───────────────────
YEAR = 2019
OUT_BASE = Path(r"D:\DATA_NATURE_SIF\EDGAR\stacked_0.1_fixed") / str(YEAR)

BASE_URL = (
    "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/EDGAR/"
    "datasets/v81_FT2022_AP_new/monthly/PM10/"
)
# ──────────────────────────────────────────────────────

SECTOR_RE = re.compile(r'href="(bkl_[^/]+)/"', re.I)
ZIP_RE = re.compile(r'href="([^"]+_emi_nc\.zip)"', re.I)


def html(url: str) -> str:
    with requests.get(url, timeout=30) as r:
        r.raise_for_status()
        return r.text


def sector_dirs() -> List[str]:
    return sorted({*SECTOR_RE.findall(html(BASE_URL))})


def emi_zip_name(sector: str) -> str:
    page = html(f"{BASE_URL}{sector}/")
    zips = ZIP_RE.findall(page)
    if not zips:
        raise RuntimeError(f"No *_emi_nc.zip in {sector}")
    return zips[0]


# ───────────── download (resumable) ─────────────
def stream_dl(url: str, dest: Path, chunk: int = 1_048_576) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    done = dest.stat().st_size if dest.exists() else 0
    headers = {"Range": f"bytes={done}-"} if done else {}
    with requests.get(url, stream=True, headers=headers, timeout=60) as r:
        if r.status_code == 416:
            return
        r.raise_for_status()
        total = done + int(r.headers.get("content-length", 0))
        mode = "ab" if done else "wb"
        with open(dest, mode) as fh, tqdm.tqdm(
            total=total,
            unit="B",
            unit_scale=True,
            initial=done,
            desc=f"⬇ {dest.name}",
        ) as bar:
            for blk in r.iter_content(chunk):
                fh.write(blk)
                bar.update(len(blk))


# ───────────── unzip only the 2019 NetCDF ─────────────
def unzip_nc_year(zip_path: Path, out_dir: Path, year: int) -> Path:
    import zipfile
    with zipfile.ZipFile(zip_path) as zf:
        member = next(
            (m for m in zf.namelist() if m.endswith(".nc") and f"_{year}_" in m), None
        )
        if member is None:
            raise RuntimeError(f"No NetCDF with _{year}_ inside {zip_path.name}")
        out = out_dir / Path(member).name
        out.parent.mkdir(parents=True, exist_ok=True)
        if not out.exists():
            zf.extract(member, out.parent)
    return out


# ───────────── build date / datesec ─────────────
def wrf_dates(time: xr.DataArray) -> tuple[xr.DataArray, xr.DataArray]:
    # pandas ≥2.2: .to_timestamp() already gives 1st-of-month 00:00
    starts = pd.to_datetime(time.values).to_period("M").to_timestamp()
    date_int = np.asarray([int(d.strftime("%Y%m%d")) for d in starts], np.int32)
    return xr.DataArray(date_int, dims="Time"), xr.DataArray(
        np.zeros_like(date_int, np.int32), dims="Time"
    )


# ───────────── PM10 variable search ─────────────
def find_pm10_var(ds: xr.Dataset) -> str | None:
    for v in ds.data_vars:
        if str(ds[v].attrs.get("substance", "")).upper() == "PM10":
            return v
    pats = [
        r"^emi[sl]?_?pm10$",
        r"^emi.*pm10",
        r"pm10.*emi",
        r"pm10$",
    ]
    for pat in pats:
        for v in ds.data_vars:
            if re.match(pat, v, re.I):
                return v
    return None


# ───────────── cell area (lat only) ─────────────
def cell_area_lat(lat_deg: np.ndarray, dlat=0.1, dlon=0.1) -> np.ndarray:
    R = 6_371_000.0
    dphi = math.radians(dlat)
    dlambda = math.radians(dlon)
    return (R**2) * dphi * dlambda * np.cos(np.radians(lat_deg))


# ───────────── rename dims & coords safely ─────────────
def rename_dims_to_wrf(ds: xr.Dataset) -> xr.Dataset:
    coord_map = {"time": "Time", "lat": "south_north", "lon": "west_east"}
    ds = ds.rename({k: v for k, v in coord_map.items() if k in ds.coords})
    ds = ds.rename_dims({k: v for k, v in coord_map.items() if k in ds.dims})
    return ds


# ───────────── per-sector conversion ─────────────
def convert_nc(src: Path, dst: Path, sector: str) -> None:
    if dst.exists():
        return
    dst.parent.mkdir(parents=True, exist_ok=True)
    print(f"↪ converting {src.name}")
    ds = xr.open_dataset(src)
    pm10 = find_pm10_var(ds)
    if pm10 is None:
        print("PM10 variable not found. Vars/attrs:")
        for v in ds.data_vars:
            print("   •", v, ds[v].attrs)
        raise RuntimeError(f"No PM10 var in {src.name}")
    print("PM10 var:", pm10)
    # ---- Tonnes → kg m⁻² s⁻¹ -------------------------------------------
    if str(ds[pm10].attrs.get("units", "")).lower() == "tonnes":
        print("   • converting Tonnes → kg m⁻² s⁻¹")
        kg = ds[pm10] * 1_000.0
        month_sec = np.array(
            [pd.Period(t, "M").end_time.day * 86_400 for t in pd.to_datetime(ds["time"])]
        )
        kg_s = kg / xr.DataArray(month_sec, dims="time")
        area = xr.DataArray(
            cell_area_lat(ds["lat"].values), dims="lat", coords={"lat": ds["lat"]}
        )
        ds[pm10] = kg_s / area
        ds[pm10].attrs["units"] = "kg m-2 s-1"
    # -------------------------------------------------------------------
    ds = rename_dims_to_wrf(ds)
    print("   • dims after rename:", dict(ds.sizes))
    if "bottom_top" not in ds.dims:
        ds = ds.expand_dims(bottom_top=[0])

    ds["date"], ds["datesec"] = wrf_dates(ds["Time"])
    new_var = f"E_PM10_{sector.lower().replace('bkl_', '')}"
    ds = ds.rename({pm10: new_var})
    ds[new_var].attrs["long_name"] = f"{YEAR} PM10 {sector}"
    ds.attrs["YEAR"] = str(YEAR)
    print("writing", dst.relative_to(OUT_BASE))
    ds.to_netcdf(dst, format="NETCDF4_CLASSIC")
    ds.close()
    print("done")


# ───────────── build total file ─────────────
def make_total(anth_dir: Path, out_nc: Path) -> None:
    if out_nc.exists():
        return
    parts = []
    for p in anth_dir.glob("bkl_*/*.nc"):
        with xr.open_dataset(p) as d:
            parts.append(d[next(v for v in d.data_vars if v.startswith("E_PM10_"))])
    if not parts:
        print("No sector files → total skipped")
        return
    total = xr.combine_by_coords(parts).sum("variable")
    ds = xr.Dataset({"E_PM10_total": total})
    ds["E_PM10_total"].attrs.update(
        {"units": "kg m-2 s-1", "long_name": f"{YEAR} PM10 total emissions"}
    )
    ds["date"], ds["datesec"] = wrf_dates(total["Time"])
    ds.attrs["TITLE"] = f"EDGAR v8.1 PM10 total for {YEAR}"
    ds.to_netcdf(out_nc, format="NETCDF4_CLASSIC")
    print("total file:", out_nc.relative_to(OUT_BASE))


# ───────────── main driver ─────────────
def main() -> None:
    RAW = OUT_BASE / "raw_zip"
    NC = OUT_BASE / "extracted_nc"
    OUT = OUT_BASE / "anthro_emiss"
    for sector in sector_dirs():
        print(f"\n=== {sector} ===")
        zip_n = emi_zip_name(sector)
        zip_u = f"{BASE_URL}{sector}/{zip_n}"
        zip_p = RAW / sector / zip_n
        if not zip_p.exists():
            stream_dl(zip_u, zip_p)
        else:
            print(f"✔ {zip_p.relative_to(OUT_BASE)} exists – skip download")
        nc_in = unzip_nc_year(zip_p, NC / sector, YEAR)
        nc_out = OUT / sector / nc_in.name.replace("_emi", "_anthro_emiss")
        convert_nc(nc_in, nc_out, sector)
    make_total(OUT, OUT / f"PM10_total_{YEAR}_anthro_emiss.nc")


# ───────────── entrypoint ─────────────
if __name__ == "__main__":
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(1)

