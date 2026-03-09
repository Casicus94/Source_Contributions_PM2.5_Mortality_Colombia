import os
import numpy as np
import xarray as xr
from datetime import datetime
import pdb 

# Constantes
R_EARTH = 6371000  # Radius in meters
SECONDS_IN_YEAR = 365 * 24 * 3600
TONNES_TO_KG = 1e3
START_DATE = "2019-01-01"
WRF_REFERENCE_DATE = "2019-01-01"
TIME_VALUE = 0  # required for interpolation
DATE_INT = int(WRF_REFERENCE_DATE.replace("-", ""))
DATESEC_INT = 0

# Directorios
path = '/ptmp/acasallasg/Sectors/'
#base_dirs = [
#    "TOTALS", "TOTALS_No_FIRES", "TOTALS_No_ANTHRO_FIRES",
#    "Agriculture", "Aviation", "Industry", "Off_Road_Transport",
#    "Road_Transport", "Shipping", "Waste", "Energy_Use"
#]

base_dirs = ["TOTALS", "Agriculture"]
species = ["BC", "CO", "NH3", "NMVOC", "NOx", "OC", "PM10", "PM2.5", "SO2"]

def compute_cell_area(latitudes, dlat_deg, dlon_deg):
    dlat_rad = np.deg2rad(dlat_deg)
    dlon_rad = np.deg2rad(dlon_deg)
    lat_rad = np.deg2rad(latitudes)
    area = R_EARTH**2 * dlat_rad * dlon_rad * np.cos(lat_rad)
    return area  # m²

def convert_to_wrf_format(file_path, folder):
    ds = xr.open_dataset(file_path)

    if "emissions" in ds:
        ds = ds.rename({"emissions": "emis_tot"})
    elif "emis_tot" not in ds:
        print(f"'emis_tot' or 'emissions' not found in: {file_path}")
        return

    emis = ds["emis_tot"]

    lat = ds["lat"]
    lon = ds["lon"]
    dlat = abs(lat[1] - lat[0]).item()
    dlon = abs(lon[1] - lon[0]).item()

    # Área dinámica por latitud
    cell_areas = compute_cell_area(lat.values, dlat, dlon)
    cell_areas_2d = np.repeat(cell_areas[:, np.newaxis], len(lon), axis=1)

    # Convertir t/año → kg/m²/s
    if folder == "TOTALS" or folder == "TOTALS_No_ANTHRO_FIRES" or folder == 'TOTALS_No_FIRES': 
        emis_kg_m2_s = (emis.values * TONNES_TO_KG) / (cell_areas_2d * SECONDS_IN_YEAR)
    else:
        emis_kg_m2_s = (emis[0].values * TONNES_TO_KG) / (cell_areas_2d * SECONDS_IN_YEAR)
    emis_kg_m2_s[np.isnan(emis_kg_m2_s)] = -9999.0
    emis_kg_m2_s = np.expand_dims(emis_kg_m2_s.astype("float32"), axis=0)
    #emis_kg_m2_s = emis_kg_m2_s.transpose(0, 2, 1)

    # Nueva variable emis_tot
    emis_da = xr.DataArray(
        emis_kg_m2_s,
        dims=("time", "lat", "lon"),
        coords={"time": [TIME_VALUE], "lat": lat, "lon": lon},
        attrs={
            "units": "kg m^-2 s^-1",
            "long_name": "Emissions of species",
            "standard_name": "tendency_of_atmosphere_mass_content_due_to_emission",
            "cell_methods": "time: mean (interval: 1 month, 31 days)",
            "comment": "Converted from EDGAR original data with dynamic cell area",
            "_FillValue": -9999.0
        }
    )

    # Variables de tiempo
    time = xr.DataArray(np.array([TIME_VALUE], dtype="int32"), dims="time", attrs={
        "units": f"days since {WRF_REFERENCE_DATE}",
        "long_name": "time",
        "calendar": "gregorian"
    })

    date = xr.DataArray(np.array([DATE_INT], dtype="int32"), dims="time", attrs={
        "units": "YYYYMMDD",
        "long_name": "Date in YYYYMMDD format"
    })

    datesec = xr.DataArray(np.array([DATESEC_INT], dtype="int32"), dims="time", attrs={
        "units": "seconds since midnight",
        "long_name": "Seconds since midnight"
    })

    # Nuevo dataset
    wrf_ds = xr.Dataset({
        "emis_tot": emis_da,
        "time": time,
        "date": date,
        "datesec": datesec,
        "lat": lat,
        "lon": lon
    })

    wrf_ds["time"].attrs = {
        "units": f"days since {WRF_REFERENCE_DATE}",
        "long_name": "time",
    }

    wrf_ds.attrs = {
        "title": "Anthropogenic Emission Data",
        "source": "Converted from EDGAR raster",
        "history": f"Converted on {datetime.now().strftime('%Y-%m-%d')}",
        "comment": "This file was created for use with WRF-CHEM."
    }
    # Fix lat/lon long_name and remove _FillValue
    wrf_ds["lat"].attrs["long_name"] = "lat"
    wrf_ds["lon"].attrs["long_name"] = "lon"
    wrf_ds["lat"].attrs.pop("_FillValue", None)
    wrf_ds["lon"].attrs.pop("_FillValue", None)

    # Optional: match comment to WRF format
    wrf_ds["emis_tot"].attrs["comment"] = "Data based on processed raster stack."

    # Guardar en mismo lugar con "_fix.nc"
    output_path = file_path.replace(".nc", "_fix.nc")
    wrf_ds.to_netcdf(output_path, format="NETCDF4_CLASSIC", unlimited_dims=["time"])
    print(f" Archivo guardado: {output_path}")

# Procesamiento
for folder in base_dirs:
    print(folder)
    for sp in species:
        print(sp)
        if folder.startswith("TOTALS"):
            filename = f"{path}{folder}/{sp}_TOTALS.nc"
        else:
            filename = f"{path}{folder}/{sp}_{folder}_subtracted_from_TOTALS.nc"

        if os.path.exists(filename):
            convert_to_wrf_format(filename, folder)
        else:
            print(f" Archivo no encontrado: {filename}")

