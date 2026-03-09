import os
import xarray as xr

# Define the grouped sectors
sector_groups = {
    "Agriculture": ["AGS", "AWB", "ENF", "MNM"],
    "Industry": ["IND", "CHE", "IRO", "NFE", "PRU_SOL"],
    "Energy_Use": ["ENE", "REF_TRF", "RCO"],
    "Fuel_Exploitation": ["PRO_COAL", "PRO_GAS", "PRO_OIL", "FFF"],
    "Waste": ["SWD_INC", "SWD_LDF", "WWT"],
    "Road_Transport": ["TRO"],
    "Aviation": ["TNR_Aviation_CDS", "TNR_Aviation_CRS", "TNR_Aviation_LTO"],
    "Shipping": ["TNR_Ship"],
    "Off_Road_Transport": ["TNR_Other"],
    "Indirect_Emissions": ["IDE", "N2O"],
}

def subtract_sector_groups(base_dir, output_dir, pollutants, year=2019, var_name='emissions'):
    os.makedirs(output_dir, exist_ok=True)

    for pollutant in pollutants:
        print(f"\nProcessing pollutant: {pollutant}")
        
        # Path to TOTALS file
        total_path = os.path.join(base_dir, pollutant, "TOTALS", f"v8.1_FT2022_AP_{pollutant}_{year}_TOTALS_emi.nc")
        if not os.path.exists(total_path):
            print(f"TOTALS file missing for {pollutant} — skipping")
            continue

        try:
            ds_total = xr.open_dataset(total_path)
        except Exception as e:
            print(f"Could not open TOTALS file for {pollutant}: {e}")
            continue

        for group_name, group_sectors in sector_groups.items():
            ds_group_sum = None
            found_any = False

            # Sum all available sector files in the group
            for sector in group_sectors:
                sector_file = os.path.join(base_dir, pollutant, sector, f"v8.1_FT2022_AP_{pollutant}_{year}_{sector}_emi.nc")
                if not os.path.exists(sector_file):
                    print(f"  Missing {sector} for {pollutant} — skipped")
                    continue

                try:
                    ds_sector = xr.open_dataset(sector_file)
                    if ds_group_sum is None:
                        ds_group_sum = ds_sector[var_name]
                    else:
                        ds_group_sum += ds_sector[var_name]
                    found_any = True
                except Exception as e:
                    print(f"  Error reading {sector}: {e}")

            if not found_any:
                print(f"No sectors found for group {group_name} and pollutant {pollutant} — skipping group")
                continue

            # Subtract group sum from TOTALS
            ds_diff = ds_total.copy()
            ds_diff[var_name] = ds_total[var_name] - ds_group_sum

            # Save output
            group_out_dir = os.path.join(output_dir, pollutant, group_name)
            os.makedirs(group_out_dir, exist_ok=True)
            out_path = os.path.join(group_out_dir, f"{pollutant}_{group_name}_subtracted_from_TOTALS.nc")

            ds_diff.to_netcdf(out_path)
            print(f"Saved: {out_path}")

# Example usage:
if __name__ == "__main__":
    base_dir = "/nfs/scistore16/mullegrp/acasalla/MPI-C/Emissions/"
    output_dir = "/nfs/scistore16/mullegrp/acasalla/MPI-C/Emissions/Subtractions"
    pollutants = ["NOx", "SO2", "NMVOC", "NH3", "CO", "BC", "OC", "PM10", "PM2.5"]
    subtract_sector_groups(base_dir, output_dir, pollutants)

import shutil

def copy_totals_only(base_dir, destination_dir, pollutants, year=2019):
    totals_root = os.path.join(destination_dir, "TOTALS")
    os.makedirs(totals_root, exist_ok=True)

    for pollutant in pollutants:
        src = os.path.join(base_dir, pollutant, "TOTALS", f"v8.1_FT2022_AP_{pollutant}_{year}_TOTALS_emi.nc")
        dst_dir = os.path.join(totals_root, pollutant)
        dst = os.path.join(dst_dir, f"{pollutant}_TOTALS.nc")

        if os.path.exists(src):
            os.makedirs(dst_dir, exist_ok=True)
            shutil.copy2(src, dst)
            print(f"Copied {pollutant} → {dst}")
        else:
            print(f"TOTALS file for {pollutant} not found, skipping")

base_dir = "/nfs/scistore16/mullegrp/acasalla/MPI-C/Emissions"
destination_dir = "/nfs/scistore16/mullegrp/acasalla/MPI-C/Emissions"

copy_totals_only(base_dir, destination_dir, pollutants)
