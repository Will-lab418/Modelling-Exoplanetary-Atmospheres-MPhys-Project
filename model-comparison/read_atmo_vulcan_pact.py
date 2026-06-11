import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from netCDF4 import Dataset
from matplotlib.lines import Line2D

# This script reads outputs for a P-T profile and chemical abundances from ATMO, VULCAN, and PACT
# and plots them together for easier comparison.

# Plots P-T profiles for ATMO, VULCAN, and PACT

def plot_pt(atmo_file, vulcan_file, pact_file, tit = 'HD209 P-T'):

    linewidth = 1.2

    linestyles = ['-', '--', '-.']
    files = [atmo_file, vulcan_file, pact_file]
    models = ['ATMO', 'VULCAN', 'PACT']

    # reads pressure and Temperature from each model
    def read_data(model, fpath):

        if model == 'ATMO':
            nc = Dataset(fpath, 'r')
            tt = np.array(nc.variables['temperature'][:], dtype = float)
            pp = np.array(nc.variables['pressure'][:], dtype = float) / 1e6 # converts pressure to bar

            nc.close()
            return tt, pp
        
        elif model == 'VULCAN':
            with open(fpath, 'rb') as handle:
                data = pickle.load(handle)
            tt = np.array(data['atm']['Tco'], dtype = float)
            pp = np.array(data['atm']['pco'], dtype = float) / 1e6 # converts pressure to bar

            return tt, pp 
        
        elif model == 'PACT':

            # skips lines starting with !
            with open(fpath, 'r') as fh:
                lines = fh.readlines()
            header_idx = next(i for i, line in enumerate(lines) if not line.startswith('!'))
            columns = lines[header_idx].split()
            arr = np.loadtxt(fpath, comments='!', skiprows=header_idx+1)
            col_idx = {name: i for i, name in enumerate(columns)}
            tt = arr[:, col_idx['temperature[K]']]
            pp = arr[:, col_idx['pressure[bar]']]

            return tt, pp

    fig, ax = plt.subplots(figsize=(10, 8))

    for i, (fpath, model) in enumerate(zip(files, models)):
        tt, pp = read_data(model, fpath)
        ax.semilogy(tt, pp, label = model, linestyle = linestyles[i], linewidth = linewidth)
            
    ax.set_xlim(600, 2000)
    ax.set_ylim(1e2, 1e-5)

    ax.set_xlabel('Temperature (K)', fontsize = 20)
    ax.set_ylabel('Pressure (bar)', fontsize = 20)
    ax.set_title(tit, fontsize = 20)

    ax.tick_params(labelsize = 16)
    ax.legend(frameon = False, fontsize = 16)

    plt.tight_layout()
    output = tit.replace(' ', '_') + '.png'
    plt.savefig(output, dpi=300, bbox_inches='tight')
    plt.close()
    print('Created:', output)




# Plots chemical abundances for ATMO, VULCAN, PACT, and FRECKLL

def plot_abundances(atmo_file, vulcan_file, pact_file, tit = 'HD 209 Chemical Kinetics', species = ['H2O', 'H', 'CO', 'CH4', 'NH3', 'H2', 'He']):

    linewidth = 1.2
    pact_log10 = True

    # strip b' from molecule names in ATMO
    def decode(row):
        return row.tobytes().decode('utf-8', errors='ignore').replace('\x00','').strip()
    
    # same for FRECKLL
    def decode_species(arr):
        out = []
        for x in arr:
            if isinstance(x, bytes):
                out.append(x.decode('utf-8', errors='ignore').strip())
            else:
                out.append(str(x).strip())
        return out
    
    species = list(species)

    # read ATMO
    nc_atmo = Dataset(atmo_file, 'r')
    ab_atmo = np.array(nc_atmo.variables['abundances'][:, :])
    p_atmo = np.array(nc_atmo.variables['pressure'][:]) / 1e6 # convert pressure to bar
    molenames = nc_atmo.variables['molname'][:, :]
    nc_atmo.close()

    atmo_idx = {decode(molenames[i]): i for i in range(molenames.shape[0])}

    # read VULCAN
    with open(vulcan_file, 'rb') as f:
        vulcan_data = pickle.load(f)
    ymix = np.array(vulcan_data['variable']['ymix'])
    p_vul = np.array(vulcan_data['atm']['pco']) / 1e6 # convert pressure to bar
    vul_idx = {m: i for i, m in enumerate(vulcan_data['variable']['species'])}

    # read PACT
    with open(pact_file, 'r') as fh:
        lines = fh.readlines()

    # commenting lines that start with !
    header_idx = next(i for i, line in enumerate(lines) if not line.startswith('!'))
    columns = lines[header_idx].split()
    arr = np.loadtxt(pact_file, comments='!', skiprows=header_idx+1)
    col_idx = {name: i for i, name in enumerate(columns)}
    p_pact = arr[:, col_idx['pressure[bar]']]
    meta_cols = {'height[km]', 'pressure[bar]', 'temperature[K]'}
    pact_abundances = {}
    for sp in [c for c in columns if c not in meta_cols]:
        vals = arr[:, col_idx[sp] ]
        pact_abundances[sp] = 10**vals if pact_log10 else vals

    # plotting the models
    fig, ax = plt.subplots(figsize = (10, 8 ))
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    for i, mol in enumerate(species):
        c = colors[i % len(colors)]
        if mol in atmo_idx:
            ax.plot(ab_atmo[atmo_idx[mol], :], p_atmo, linewidth = linewidth, linestyle = '-', color = c)
        if mol in vul_idx:
            ax.plot(ymix[:, vul_idx[mol]], p_vul, linewidth = linewidth, linestyle = '--', color = c)
        if mol in pact_abundances:
            ax.plot(pact_abundances[mol], p_pact, linewidth = linewidth, linestyle = '-.', color = c)

    ax.set_yscale('log')
    ax.set_xscale('log')

    ax.set_xlim(1e-9, 1e0)
    ax.set_ylim(1e2, 1e-5)

    ax.set_xlabel('Abundances', fontsize = 20)
    ax.set_ylabel('Pressure (bar)', fontsize = 20)
    ax.set_title(tit, fontsize = 20)
    ax.tick_params(labelsize = 16)

    # legend for chemical species
    species_handles = [Line2D([0], [0], color = colors[i % len(colors)], linewidth = linewidth) for i in range(len(species))]
    leg1 = ax.legend(species_handles, species, loc = 'lower left', fontsize = 14, frameon = False)
    ax.add_artist(leg1)

    # legend for the models
    model_handles = [
        Line2D([0], [0], color = 'k', linewidth = linewidth, linestyle = '-'),
        Line2D([0], [0], color = 'k', linewidth = linewidth, linestyle = '--'),
        Line2D([0], [0], color = 'k', linewidth = linewidth, linestyle = '-.'),]

    ax.legend(model_handles, ['ATMO', 'VULCAN', 'PACT'], loc = 'lower left', bbox_to_anchor = (0, 0.30), fontsize = 14, frameon = False)

    plt.tight_layout()
    output = tit.replace(' ', '_') + '.png'
    plt.savefig(output, dpi=300, bbox_inches='tight')
    plt.close()
    print('Created:', output)
