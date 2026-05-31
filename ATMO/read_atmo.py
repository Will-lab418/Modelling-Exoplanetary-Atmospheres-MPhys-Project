import matplotlib
import numpy
matplotlib.use("Agg")
import os
from matplotlib.lines import Line2D
from netCDF4 import *
from constant import *
from pylab import *
from numpy import *
import matplotlib as mpl
import matplotlib.font_manager as mpfm

mpl.rc('text',usetex=False)    
mpl.rc('font',family='serif',size=12)

prop = mpfm.FontProperties(size=12)

def plot_pt(nfig=1, file='t1000_g5_H2He.ncdf', tit='', clear=False, labels=None, cond='', color='', linewidth=1.,linestyles=None, marker=''):

    # plots the Pressure-Temperature profile from an netCDF output

    # Allow single file or list of files
    files = list(file) if isinstance(file, (list, tuple)) else [file]
    nfiles = len(files)

    # Labels for legend
    if labels is None:
        labels = [os.path.basename(f) for f in files]
    elif not isinstance(labels, (list, tuple)):
        labels = [labels] * nfiles

    # Linestyles
    if linestyles is None:
        linestyles = ['-', '--', '-.', ':']
    linestyles = (linestyles * ((nfiles // len(linestyles)) + 1))[:nfiles]

    figure(nfig, figsize=(15, 10))
    if clear:
        clf()

    # --- Plot each file ---
    for fidx, f in enumerate(files):
        if not os.path.exists(f):
            print("Skipping missing file:", f)
            continue

        ncdfFile = Dataset(f, 'r')
        vars = ncdfFile.variables

        tt = vars['temperature'][:]
        pp = vars['pressure'][:] / 1.0E6
        ncdfFile.close()

        if color == '':
            semilogy(tt, pp, label=labels[fidx],
                     lw=linewidth, linestyle=linestyles[fidx], marker=marker)
        else:
            semilogy(tt, pp, label=labels[fidx],
                     color=color, lw=linewidth,
                     linestyle=linestyles[fidx], marker=marker)

    # --- Condensation curves (use last pp grid) ---
    if cond.lower() in ['h2o', 'both']:
        semilogy(10000./(38.84-3.83*numpy.log10(pp)), pp,
                 label='H2O condensation', color='red', lw=linewidth)

    if cond.lower() in ['nh3', 'both']:
        semilogy(10000./(68.02-6.31*numpy.log10(pp)), pp,
                 label='NH3 condensation', color='magenta', lw=linewidth)

    title(tit, fontsize=25)
    ylim((pp.max(), pp.min()))
    xlabel('Temperature [K]', fontsize=25)
    ylabel('Pressure [bar]', fontsize=25)
    legend(prop=prop).draw_frame(0)

    safe_title = tit.replace(" ", "_")
    outfile = safe_title + ".png"
    savefig(outfile)
    close()

def plot_abundances(nfig=2, file='chem_out.ncdf', tit='', clear=True, linewidth=1., xlimits=[], atmo_imol=[1,4,9,12,34,62,68], linestyles=None, labels=None, legend_loc_species='lower left',legend_loc_files='lower left'): 

    # plots chemical abundances as a function of pressure from an netCDF file

    # Allow single file or list of files
    files = list(file) if isinstance(file, (list, tuple)) else [file]
    nfiles = len(files)
 
    # Linestyles for each input file 
    if linestyles is None: 
        linestyles = ['-', '--', '-.', ':'] 
    linestyles = (linestyles * ((nfiles // len(linestyles)) + 1))[:nfiles] 
            
    figure(nfig, figsize=(15, 10)) 
    if clear: clf() 
    sp = subplot(111) 
            
    # Read species name from first file 
    nc0 = Dataset(files[0], 'r') 
    vars0 = nc0.variables 
    nmol0 = len(nc0.dimensions['nmol']) 
    
    if atmo_imol == []: 
        imol = array(range(nmol0)) 
    else: imol = array(atmo_imol) - 1 

    nm0 = vars0['molname'][:, :] 
    nc0.close() 

    # clean molecule names 
    molnames = [nm0[idx, :].tobytes().decode('utf-8', errors='ignore').replace('\x00', '').strip() for idx in imol]

    # plot abundances for each file
    for fidx, f in enumerate(files): 
        if not os.path.exists(f): 
            print("Skipping missing file:", f) 
            continue 

        ncdfFile = Dataset(f, 'r') 
        vars = ncdfFile.variables 
        
        ab = vars['abundances'][:, :] 
        pp = vars['pressure'][:] / 1.0E6 # convert to bar 
        ncdfFile.close()
        
        ls = linestyles[fidx] 
        
        for j, mol in enumerate(molnames): 
            idx = imol[j] 
            plot(ab[idx, :], pp, lw=linewidth, linestyle=ls) 
            
    # Set log-log scale and labels
    sp.set_yscale('log')
    sp.set_xscale('log')
    title(tit, fontsize=25)

    ylim((pp.max(), pp.min()))
    if xlimits:
        xlim(xlimits)

    xlabel('Abundances', fontsize=25)
    ylabel('Pressure [bar]', fontsize=25)
    sp.tick_params(axis='both', labelsize=25)
            
    species_handles = [Line2D([0], [0], lw=linewidth, linestyle='-')
                   for m in molnames]
    file_handles = [Line2D([0], [0], color='k', lw=linewidth, linestyle=linestyles[i])
                for i in range(nfiles)]

    # Species legend
    leg1 = sp.legend(species_handles, molnames,loc='lower left',fontsize=16,frameon=False,borderaxespad=0.0)
    sp.add_artist(leg1)

    # File legend (only if multiple files AND labels provided)
    if nfiles > 1:
        if labels is None:
            labels = [os.path.basename(f) for f in files]

        leg2 = sp.legend(file_handles, labels,loc='lower left',bbox_to_anchor=(0, 0.20),fontsize=16,frameon=False,borderaxespad=0.0)
        sp.add_artist(leg2)
         
    safe_title = tit.replace(" ", "_") 
    outfile = safe_title + ".png" 
    savefig(outfile) 
    close()


def plot_spectrum(nfig = 3, file = 't700_spectrum.ncdf', tit = '', clear = False, label = '', coeff = 1., color = '', temp = False, nbin = 1, linewidth = 1., xlimits=None):

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

    safe_title = tit.replace(" ", "_")
    outfile = safe_title + ".png"
    savefig(outfile)
    close()

def plot_trans_spec(nfig = 4, file = 'transmision_spectrum.ncdf', tit = '',clear = False,label = '',color = '',xunit='wn',linewidth = 1.):

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
    
    safe_title = tit.replace(" ", "_")
    outfile = safe_title + ".png"
    savefig(outfile)
    close()
