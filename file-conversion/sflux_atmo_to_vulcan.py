import sys
import numpy as np
from netCDF4 import Dataset

# This script converts a UV stellar flux file from photon flux units used in ATMO to energy flux units used in VULCAN. 
# The input file is expected to be in NetCDF format with variables Wavelength (in nm) and Hnu (photon flux in photons/cm^2/s/nm/steradian). 
# The output is a text file with two columns: wavelength (nm) and energy flux (ergs/cm^2/s/nm).

# This script is run from the command line as:
# python sflux_atmo_to_vulcan.py input.ncdf output.txt

input_file  = sys.argv[1]
output_file = sys.argv[2]

# physical constants
h = 6.62607015e-27 # Planck constant (erg s)
c = 2.99792458e10 # speed of light (cm s^-1)
au = 1.49597871E13  # Astronomical Unit (cm)
#r_sun = 6.957E10 # solar radius (cm)
#r_star = 0.939 # stellar radius in solar radius (WASP-39)
#orbital_radius = 0.048 # orbital radius in AU (WASP-39b)

# read input file
src = Dataset(input_file, "r")
wavelength_nm = src.variables["Wavelength"][:]
photon_flux = src.variables["Hnu"][:]
src.close()

# convert wavelength from nm to cm
wavelength_cm = wavelength_nm * 1e-7

# convert photon flux to energy flux
energy_per_photon = h * c / wavelength_cm
energy_flux = photon_flux * energy_per_photon

# convert per cm to per nm
energy_flux_per_nm = energy_flux * 1e-7

# IGNORE THIS CODE
#scale_factor = ((orbital_radius * au)/(r_star * r_sun))**2
# energy_flux_per_nm_scaled = energy_flux_per_nm * scale_factor

energy_flux_per_nm_scaled = energy_flux_per_nm * 4 * np.pi

# write text file
with open(output_file, "w") as f:
    f.write("# WL(nm)	 Flux(ergs/cm**2/s/nm)\n")
    for wavelength, flux in zip(wavelength_nm, energy_flux_per_nm_scaled):
        f.write(f"{wavelength:.10f} {flux:.4e}\n")

print(f'Created: {output_file}')
