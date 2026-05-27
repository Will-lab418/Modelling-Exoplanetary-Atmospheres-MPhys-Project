import numpy as np
from scipy.io import netcdf_file
import sys

# This script converts the P-T profile from an netcdf file from ATMO to a text file suitable for VULCAN input. 
# This script is intended to be run from the command line as:
# python atmo_to_vulcan_pt.py input_file output_file

input_file = sys.argv[1]
output_file = sys.argv[2]

# Load NetCDF file
nc = netcdf_file(input_file, "r", mmap=False)

pressure = np.array(nc.variables["pressure"].data,    dtype=float)  # dyne/cm2
temperature = np.array(nc.variables["temperature"].data, dtype=float)  # K

# Sort descending pressure for VULCAN input
idx = np.argsort(pressure)[::-1]
pressure = pressure[idx]
temperature = temperature[idx]

# Write the VULCAN P-T input file
with open(output_file, "w") as f:
    f.write("# (dyne/cm2)    (K)\n")
    f.write("Pressure      Temp\n")
    for p, T in zip(pressure, temperature):
        f.write(f"{p:.3e}   {T:.1f}\n")

print(f"Done! Created: {output_file}")
