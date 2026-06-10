import sys
import numpy as np
import h5py
from netCDF4 import Dataset

# Converts .h5 output file from FRECKLL into an .ncdf file so it can be plotted in
# read_all.py on the process secure shell

# This script is run from the command line:
# python h5_to_ncdf.py input.h5 output.ncdf

input_file = sys.argv[1]
output_file = sys.argv[2]

# read .h5 file
with h5py.File(input_file, 'r') as f:
    pressure = f['pressure'][:]
    temperature = f['temperature'][:]
    vmrs = f['solution/vmrs'][:]
    times = f['solution/times'][:]
    species = f['species/species_formula'][:]

# decode species names from bytes
species = [s.decode('utf-8').strip() if isinstance(s, bytes) else str(s).strip() for s in species]

# create .ncdf file
nc = Dataset(output_file, 'w')
nc.createDimension('level', len(pressure))
nc.createDimension('species', len(species))
nc.createDimension('time', len(times))

nc.createVariable('pressure', 'f8', ('level',))[:] = pressure
nc.createVariable('temperature', 'f8', ('level',))[:] = temperature
nc.createVariable('vmrs', 'f8', ('time', 'species', 'level'))[:] = vmrs
nc.createVariable('times', 'f8', ('time',))[:] = times

species_v = nc.createVariable('species_formula', str, ('species',))
for i, s in enumerate(species):
    species_v[i] = s

nc.close()
print('Created:', output_file)
