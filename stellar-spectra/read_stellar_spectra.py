import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

# This script plots the stellar spectra from BT-Settl and UV flux .ncdf files
# UV flux file should be converted to the correct energy flux units first using the photo_to_energy_flux.py python script

# The Script has been edited to plot the flux in units that match VULCAN's sflux

# physical constants
h = 6.62607015e-27 # (erg s)
c = 2.99792458e10 # (cm/s)
k_b = 1.380649e-16 # (erg/K)

def plot_spectrum(files, tit = '', labels = None, nbin = 1, alpha = 0.6):
    
    linewidth = 0.5

    # allows either a single file or multiple files
    if isinstance(files, (list, tuple)):
        file = list(files)
    else: file = [files]

    nfiles = len(files)

    # labelling the legend
    if labels is None:
        labels = [os.path.basename(f) for f in files]
    elif not isinstance(labels, (list, tuple)):
        labels = [labels] * nfiles

    # each spectra can be given its own alpha value to help clarity
    if isinstance(alpha, (list, tuple)):
        alphas = list(alpha)
    else:
        alphas = [alpha] * nfiles

    fig, ax = plt.subplots(figsize = (10,8))

    for i, fpath in enumerate(files):
    
        # read .ncdf file
        nc = Dataset(fpath, 'r')
        nu = np.array(nc.variables['nu'][:], dtype = float)
        hnu = np.array(nc.variables['hnu'][:], dtype = float)
        nc.close()

        nfreq = len(nu)

        # binning to smooth spectra
        if nbin > 1:
            nfreq_bin = nfreq // nbin 
            half = nbin // 2
            nu_plot = np.zeros(nfreq_bin)
            hnu_plot = np.zeros(nfreq_bin)

            nu_plot[0] = nu[0]
            hnu_plot[0] = np.sum(hnu[0:half])/ half

            for j in range(1, nfreq_bin):
                nu_plot[j] = nu[nbin * j]
                hnu_plot[j] = np.sum(hnu[nbin * j - half : nbin * j + half]) / nbin 
        else:
            nu_plot = nu
            hnu_plot = hnu

        # convert wavenumber to wavelength
        wavelength_cm = (1/nu_plot)

        # convert cm to nm
        wavelength_nm = wavelength_cm * 1e7

        
        # Bt-Settl and UV flux files have units per steradian
        hnu_tot = hnu_plot * 4 * np.pi 

        # convert per wavenumber to per cm
        flambda_cm = hnu_tot * nu_plot**2 #

        # convert per cm to per nm
        flambda_nm = flambda_cm * 1e-7 # erg/s/cm^2/nm

        # BT-Settl contains zeros when in log axis, replace with nan so they can be ignored
        flambda_nm[flambda_nm <=0] = np.nan

        ax.plot(wavelength_nm, flambda_nm, label=labels[i], linewidth = linewidth, alpha = alphas[i])

    ax.set_yscale('log')
    ax.set_xlim(0, 1000)
    ax.set_ylim(1e0, 1e10)

    ax.set_xlabel('Wavelength (nm)', fontsize = 20)
    ax.set_ylabel(r'Surface Spectral Flux (erg s$^{-1}$ cm$^{-2}$ nm$^{-1}$)', fontsize = 20)
    ax.set_title(tit, fontsize = 20)
    ax.tick_params(labelsize = 16)
    ax.legend(frameon = False, fontsize = 14)

    output = tit.replace(" ", "_") + ".png"
    plt.savefig(output, bbox_inches='tight')
    plt.close()
    print('Created:', output)
