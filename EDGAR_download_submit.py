import os
import requests
import zipfile

def download_edgar_2019(pollutant, sector, output_folder):
    zip_name = f"v8.1_FT2022_AP_{pollutant}_2019_{sector}_emi_nc.zip"
    url = f"https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/EDGAR/datasets/v81_FT2022_AP_new/{pollutant}/{sector}/emi_nc/{zip_name}"

    os.makedirs(output_folder, exist_ok=True)
    zip_path = os.path.join(output_folder, zip_name)

    if os.path.exists(zip_path):
        print(f"Already exists: {zip_name}")
        return

    print(f"Downloading: {zip_name}")
    response = requests.get(url, stream=True)
    if response.status_code != 200:
        print(f"Failed to download {zip_name} — HTTP {response.status_code}")
        return

    with open(zip_path, "wb") as f:
        for chunk in response.iter_content(1024 * 1024):
            f.write(chunk)
    print(f"Saved: {zip_path}")

    # Extract
    try:
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(output_folder)
            print(f"Extracted: {zip_name}")
    except zipfile.BadZipFile:
        print(f"Corrupt ZIP file: {zip_name}")


# Run the download loop
# All pollutants available
pollutants = ["NOx", "SO2", "NH3", "NMVOC", "CO", "BC", "OC", "PM10", "PM2.5"]

# Sector codes (as per EDGAR v8.1 FT2022 directory structure)
sectors = [
    "AGS", "AWB", "CHE", "ENE", "ENF", "FFF", "IDE", "IND", "IRO", "NFE", "MNM", "N2O",
    "PRO_COAL", "PRO_GAS", "PRO_OIL", "PRU_SOL", "RCO", "REF_TRF", "SWD_INC", "SWD_LDF",
    "TNR_Aviation_CDS", "TNR_Aviation_CRS", "TNR_Aviation_LTO", "TNR_Other", "TNR_Ship",
    "TRO", "WWT", "TOTALS"
]

# Destination folder
base_dir = "/nfs/scistore16/mullegrp/acasalla/MPI-C/Emissions"

for pollutant in pollutants:
    for sector in sectors:
        out_folder = os.path.join(base_dir, pollutant, sector)
        download_edgar_2019(pollutant, sector, out_folder)