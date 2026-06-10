import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from netCDF4 import Dataset
from matplotlib.lines import Line2D

# This script contains functions to plot the P-T profile, Chemical abundances as well as
# the transmission and emission spectra outputs from ATMO

# This script is a modified version of the original read_atmo.py

mpl.rc('text',usetex=False)    
mpl.rc('font',family='serif',size=12)

# Plots P-T profiles from ATMO .ncdf files
def plot_pt(nfig=1, file='t1000_g5_H2He.ncdf', tit='', clear=False, labels=None, cond='', linewidth=1.,linestyles=None, marker=''):

    # allows single or multiple files
    if isinstance(file, (list, tuple)):
        files = list(file)
    else:
       files = [file]

    nfiles = len(files)

    # Labels for legend
    if labels is None:
        labels = [os.path.basename(f) for f in files]
    elif not isinstance(labels, (list, tuple)):
        labels = [labels] * nfiles

    # gives each file a different linestyles
    if linestyles is None:
        linestyles = ['-', '--', '-.', ':']
    linestyles = (linestyles * ((nfiles // len(linestyles)) + 1))[:nfiles]

    fig, ax = plt.subplots(figsize = (10,8))

    # extracts data from file(s)
    for i, f in enumerate(files):
       nc = Dataset(f, 'r')
       tt = np.array(nc.variables['temperature'][:], dtype = float)
       pp = np.array(nc.variables['pressure'][:], dtype = float) / 1e6 # converts pressure to bar
       nc.close()

       ax.semilogy(tt, pp, label = labels[i], linewidth = linewidth, linestyle = linestyles[i])

    # Condensation curves (use last pp grid)
    if cond.lower() in ['h2o', 'both']:
        ax.semilogy(10000./(38.84-3.83*np.log10(pp)), pp,
                 label='H2O condensation', color='red', lw=linewidth)

    if cond.lower() in ['nh3', 'both']:
        ax.semilogy(10000./(68.02-6.31*np.log10(pp)), pp,
                 label='NH3 condensation', color='magenta', lw=linewidth)
        
    ax.set_ylim(pp.max(), pp.min())
    ax.set_xlabel('Temperature [K]', fontsize = 20)
    ax.set_ylabel('Pressure (bar)', fontsize = 20)
    ax.set_title(tit, fontsize = 20)
    ax.tick_params(labelsize = 16)
    ax.legend(frameon = False, fontsize = 14)

    output = tit.replace(' ', '_') + '.png'
    plt.savefig(output)
    plt.close()


# Plots the chemical abundances from ATMO .ncdf output files
def plot_abundances(file, tit = '',linewidth = 1.2, xlimits= None, atmo_imol = [1, 4, 9, 12, 34, 62, 68], linestyles = None, labels = None, legend_loc_species = 'lower left', legend_loc_files = 'lower left'):
    
    # allows single or multiple files
    if isinstance(file, (list, tuple)):
        files = list(file)
    else:
       files = [file]

    nfiles = len(files)

    # gives each file a different linestyles
    if linestyles is None:
        linestyles = ['-', '--', '-.', ':']
    linestyles = (linestyles * ((nfiles // len(linestyles)) + 1))[:nfiles]

    nc0 = Dataset(files[0], 'r') 
    vars0 = nc0.variables 
    nmol0 = len(nc0.dimensions['nmol'])

    if atmo_imol == []: 
        imol = np.array(range(nmol0)) 
    else: imol = np.array(atmo_imol) - 1

    nm0 = vars0['molname'][:, :] 
    nc0.close() 

    molnames = []
    for idx in imol:
        name = nm0[idx, :].tobytes().decode('utf-8', errors='ignore')
        molnames.append(name.replace('\x00', '').strip())

    color_cycle = mpl.rcParams['axes.prop_cycle'].by_key().get('color', [f'C{i}' for i in range(10)])
    mol_colors  = {mol: color_cycle[i % len(color_cycle)] for i, mol in enumerate(molnames)}
 
    fig, ax = plt.subplots(figsize = (10, 8))

    # extract data from .ncdf file 
    for i, f in enumerate(files):
        nc = Dataset(f, 'r')
        ab = np.array(nc.variables['abundances'][:, :])
        pp = np.array(nc.variables['pressure'][:]) / 1e6  # convert pressure to bar 
        nc.close()
 
        # plot each species with its assigned colour and file linestyle
        for j, mol in enumerate(molnames):
            ax.plot(ab[imol[j], :], pp, linewidth = linewidth, linestyle = linestyles[i], color = mol_colors[mol])
 
    ax.set_yscale('log')
    ax.set_xscale('log')

    ax.set_ylim(pp.max(), pp.min())
    if xlimits:
        ax.set_xlim(xlimits)

    ax.set_xlabel('Abundances', fontsize=20)
    ax.set_ylabel('Pressure (bar)', fontsize=20)
    ax.set_title(tit, fontsize=20)
    ax.tick_params(labelsize=16)
 
    # create species legend
    species_handles = [Line2D([0], [0], color= mol_colors[m], linewidth = linewidth, linestyle = '-') for m in molnames]
    leg1 = ax.legend(species_handles, molnames, loc = legend_loc_species, fontsize = 14, frameon = False, borderaxespad = 0.0)
    ax.add_artist(leg1)
 
    # create file legend
    if nfiles > 1:
        if labels is None:
            labels = [os.path.basename(f) for f in files]
        file_handles = [Line2D([0], [0], color= 'k', linewidth = linewidth, linestyle = linestyles[i]) for i in range(nfiles)]
        ax.legend(file_handles, labels, loc = legend_loc_files, bbox_to_anchor = (0, 0.30), fontsize = 14, frameon = False, borderaxespad = 0.0)
 
    output = tit.replace(' ', '_') + '.png'
    plt.savefig(output)
    plt.close()
            
def plot_spectrum(nfig = 3, file = 't700_spectrum.ncdf', tit = '', clear = False,
                  label = '', coeff = 1., color = '', temp = False,
                  nbin = 1, linewidth = 1., xlimits=None):

    # OPTIONAL PARAMETERS :
    # xlimits : list like [xmin, xmax] to set x-axis limits

    ncdfFile = Dataset(file,'r')
    vars = ncdfFile.variables

    nu = vars['nu'][:]
    fnu = vars['fnu'][:]
    nfreq = size(nu)

    if (not temp) :

       if (nbin == 1) :
          #conversion
          lamb = 1./(nu * 1.e-4)
          flambda = fnu/(lamb**2)*1.e4 /1000.

          figure(nfig,figsize=(15,10))
          if clear:
              clf()
          if color == '': plot(lamb, flambda*coeff,label=label,lw=linewidth)
          else :          plot(lamb, flambda*coeff,label=label,color=color,lw=linewidth)

       else :
          # arithmetic means
          nfreq_bin = nfreq // nbin
          nur = numpy.zeros((nfreq_bin))
          fnur = numpy.zeros((nfreq_bin))
          half = nbin // 2

          for i in range(1,nfreq_bin):
               nur[i] = nu[nbin*i]
               fnur[i] = sum(fnu[nbin*i-half:nbin*i+half])/nbin
          nur[0] = nu[0]
          fnur[0] = sum(fnu[0:half])/half

          #conversion
          lambr = 1./(nur*1.e-4)
          flambdar = fnur/(lambr**2)*1.e4 /1000.

          figure(nfig,figsize=(15,10))
          if clear:
              clf()
          if color == '': plot(lambr, flambdar*coeff,label=label,lw=linewidth)
          else :          plot(lambr, flambdar*coeff,label=label,color=color,lw=linewidth)

       title(tit,fontsize=25)
       xlabel('Wavelength [$\\mu$m]', fontsize=25)
       ylabel('Spectral flux [W.m-2.$\\mu$m-1]', fontsize=25)

    else :

       #conversion into brightness temperature
       tb = hplanck*c*nu/kb/numpy.log(1. + 2.*pi*hplanck*c2*nu**3/fnu)

       if (nbin == 1) :

          lamb = 1./(nu * 1.e-4)

          figure(nfig,figsize=(15,10))
          if clear:
             clf()
          if color == '': plot(lamb, tb,label=label,lw=linewidth)
          else :          plot(lamb, tb,label=label,color=color,lw=linewidth)

       else :

          nfreq_bin = nfreq // nbin
          nur = numpy.zeros((nfreq_bin))
          tbr = numpy.zeros((nfreq_bin))
          half = nbin // 2

          for i in range(1,nfreq_bin):
               nur[i] = nu[nbin*i]
               tbr[i] = sum(tb[nbin*i-half:nbin*i+half])/nbin
          nur[0] = nu[0]
          tbr[0] = sum(tb[0:half])/half

          lambr = 1./(nur*1.e-4)

          figure(nfig,figsize=(15,10))
          if clear:
             clf()
          if color == '': plot(lambr, tbr,label=label,lw=linewidth)
          else :          plot(lambr, tbr,label=label,color=color,lw=linewidth)

       title(tit,fontsize=25)
       xlabel('Wavelength [$\\mu$m]', fontsize=25)
       ylabel('Brightness temperature [K]', fontsize=25)

    # --- APPLY XLIMITS HERE ---
    if xlimits is not None:
        xlim(xlimits[0], xlimits[1])

    if not label=='':
        legend(prop=prop).draw_frame(0)

    savefig("spectrum_plot.png")
    close()

def plot_trans_spec(nfig = 4, file = 'transmision_spectrum.ncdf', tit = '',clear = False,label = '',color = '',xunit='wn',linewidth = 1.):


    # OPTIONAL PARAMETERS :
    
    # nfig : index of the figure, default 4
    # file : input ncdf effective radii file, default 't700_radii.ncdf'
    # tit : plot title, default ''
    # clear : erase previous curves if True, default False
    # label : curve label, default ''
    # color : curve color, python's default color choice if ='', default ''
    # nbin : if > 1, calculate nbin-points arithmetic means on the frequencies and the effective radii to decrease the number of points to plot
    # linewidth : line width of the curves, default 1.
    
    
    ncdfFile = Dataset(file,'r')
    vars = ncdfFile.variables

    nu = vars['nu'][:]
    Rp = vars['transit_radius'][:]
    
    figure(nfig,figsize=(15,10))
    title(tit, fontsize=25)

    if (xunit == 'wn'):

      if clear:
        clf()
      if (color == ''):
        plot(nu,Rp,label=label,lw=linewidth)  
      else:
        plot(nu,Rp,label=label,color=color,lw=linewidth)
      
      xlabel('Wavenumber [cm$^{-1}$]', fontsize=25)
      ylabel('Effective radius [$R_*$]', fontsize=25)
    
    elif (xunit == 'wl'):
      print 
      lamb = 1./(nu*1e-4)
      figure(nfig,figsize=(15,10))
      if clear:
        clf()
      if (color == ''):
        semilogx(lamb,Rp,label=label,lw=linewidth)
      else:
        semilogx(lamb,Rp,label=label,color=color,lw=linewidth)
        
      xlabel('Wavelength [$\\mu$m]', fontsize=25)
      ylabel('Effective radius [$R_*$]', fontsize=25)
    else:
      print('Unit of x-axis not recognised')
      return

    if not label=='':
        legend(prop=prop).draw_frame(0)
    
    savefig("trans_spectrum_plot.png")
    close()
