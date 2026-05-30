from numpy import *
from pylab import *
from scipy.interpolate import *
import netCDF4 as nc
import sys

# This script converts BT-Settl spectra to a format that can be read by ATMO. 
# This script has been modified from the original phoenix_convert_to_atmo.py conversion script.

def write_phoenix(fin=None,fout=None):
    
	# check input and output files
	if fin== None or fout== None:
		print('Error: input and output files must be specified')
		sys.exit()

	# read file and use first and second columns
	ff = loadtxt(fin, usecols=(0, 1))

	# sort wavelengths
	ff = ff[ff[:, 0].argsort()]

	# BT-Settl file already has linear flux. so there's no need to convert from log.
	# ff[:,1] = 10**(ff[:,1] + (-8.))

	# convert wavelength from angstrom to cm
	ff[:, 0] = ff[:, 0] * 1.0e-8
	ff[:, 1] = ff[:, 1] / 1.0e-8

	# convert erg/s/cm to erg/s/cm^2/steradian
	ff[:, 1] = ff[:, 1] / (4.0 * pi)

	wavelength = ff[:, 0]
	flux = ff[:, 1]

	# construct variables to write
	fint = interp1d(1./(wavelength[::-1]),flux[::-1]*wavelength[::-1]**2,fill_value=0.,bounds_error=False)

	nfreq = 500000 # number of points
	# maximum wavenumber
	numax = 50000.
	# wavenumber spacing
	dnu = numax/nfreq

	# compute new wavenumber grid
	wn = linspace(0.5*dnu,numax-0.5*dnu,nfreq)
	
	# interpolate flux onto new wavenumber grid
	fn = fint(wn)

	# write to netCDF
	nout = nc.Dataset(fout,'w',format='NETCDF3_CLASSIC')
	nout.createDimension('nfreq',nfreq)
	nnu=nout.createVariable('nu','f8',('nfreq',))
	nflux=nout.createVariable('hnu','f8',('nfreq',))
	nnu.units  = 'cm-1'
	nflux.units  = 'erg s-1 cm-1 ster-1'

	nnu[:]  = wn[:]
	nflux[:]  = fn[:]
	nout.close()
