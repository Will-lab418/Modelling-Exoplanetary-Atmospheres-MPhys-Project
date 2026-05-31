import sys
import numpy as np
from netCDF4 import Dataset

# This script converts a UV stellar flux file using units of nm for wavelength and
# photon flux in s^-1cm^-2 nm^-1 to energy flux units.

# This script is run from the command line as:
# python photon_to_energy_flux.py input.ncdf output.ncdf

input_file = sys.argv[1]
output_file = sys.argv[2]

# Constants
h = 6.62607015e-27 # erg s
c = 2.99792458e10 # cm s^-1

# read input file
src = Dataset(input_file, 'r')
wavelength_nm = src.variables['Wavelength'][:]
photon_flux = src.variables['Hnu'][:]
src.close()

# convert wavelength from nm to cm
wavelength_cm = wavelength_nm * 1e-7

# convert photon flux to energy flux
energy_per_photon = (h*c)/wavelength_cm
energy_flux = photon_flux * energy_per_photon

# convert per nm to per cm
energy_flux_per_cm = energy_flux * 1e7

# convert wavelength grid to wavenumber grid
nu = 1/wavelength_cm
hnu = energy_flux_per_cm * wavelength_cm**2

# sort by ascending wavenumber
order = np.argsort(nu)
nu = nu[order]
hnu = hnu[order]
wavelength_nm = wavelength_nm[order]

# write output netCDF file
dst = Dataset(output_file, 'w', format='NETCDF4')
dst.createDimension('nu', len(nu))

nu_var = dst.createVariable("nu", "f8", ("nu",))
hnu_var = dst.createVariable("hnu", "f8", ("nu",))
wave_var = dst.createVariable("wavelength_nm", "f8", ("nu",))

nu_var[:] = nu
hnu_var[:] = hnu
wave_var[:] = wavelength_nm

dst.close()

print('Created:', output_file)



