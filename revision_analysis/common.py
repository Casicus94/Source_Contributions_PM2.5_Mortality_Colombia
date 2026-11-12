"""Shared loaders for sensitivity_bias_R3M3.py (same inputs as make_table1.py).
Run from the MPIC_project folder, or set MPIC_ROOT. Provides:
  P    root folder (with trailing '/')      B   Brendan's data folder (for_alejandro_Sep29_2026/)
  D    deaths per grid cell (mort_total.nc, mort_dist)   sid  department id per grid cell
  L    L(run) -> annual-mean PM2.5 of a WRF-Chem run     names {shape_id: department name}"""
import os
import numpy as np, pandas as pd, netCDF4 as nc

P = os.environ.get("MPIC_ROOT", ".").rstrip("/") + "/"
B = P + "for_alejandro_Sep29_2026/"
L = lambda run: np.ma.filled(nc.Dataset(B + f"PM2_5_{run}_ymean_lonlat.nc")["PM2_5_DRY"][:], np.nan).astype(float).squeeze()
D = np.ma.filled(nc.Dataset(B + "mort_total.nc")["mort_dist"][:], 0.0).astype(float)
sid = np.array(nc.Dataset(P + "Departments_gridded.nc")["shape_id"][:])
_lst = pd.read_csv(P + "departments_list.csv")
names = dict(zip(_lst.iloc[:, 1].astype(int), _lst.iloc[:, 0]))
