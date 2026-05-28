import pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from netCDF4 import Dataset

# This script is used to plot the chemical kinetics of ATMO and VULCAN together

def read_atmo_vulcan(atmo_file, vul_file, tit ='ATMO_vs_VULCAN'):

    species = ['H2O', 'H', 'CO', 'CH4', 'NH3', 'H2', 'He']

    # read ATMO data
    nc_atmo = Dataset(atmo_file, 'r')
    ab_atmo = np.array(nc_atmo.variables['abundances'][:,:])
    p_atmo = np.array(nc_atmo.variables['pressure'][:]) / 1e6 # convert to bar
    molnames = nc_atmo.variables['molname'][:,:]
    nc_atmo.close()

    def decode(row):
        return row.tobytes().decode('utf-8', errors = 'ignore').replace('\x00', '').strip()
    
    atmo_idx = {decode(molnames[i]): i for i in range(molnames.shape[0])}

    # read VULCAN data
    with open(vul_file, 'rb') as f:
        vul_data = pickle.load(f)
    
    ymix = np.array(vul_data['variable']['ymix'])
    p_vul = np.array(vul_data['atm']['pco']) / 1e6 # convert to bar
    vul_idx = {m: i for i, m in enumerate(vul_data["variable"]['species'])}

    # plotting
    fig, ax = plt.subplots(figsize=(10, 8))
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    for i, mol in enumerate(species):
        c = colors[i]
        if mol in atmo_idx:
            ax.plot(ab_atmo[atmo_idx[mol], :], p_atmo, linewidth = 1.2, linestyle='-', color=c)
        if mol in vul_idx:
            ax.plot(ymix[:, vul_idx[mol]], p_vul, linewidth = 1.2, linestyle='--', color=c)

    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim(1e-9, 1e0)
    ax.set_ylim(1e2, 1e-5)
    ax.set_xlabel('Abundances', fontsize=20)
    ax.set_ylabel('Pressure (bar)', fontsize=20)
    ax.set_title(tit, fontsize=20)
    ax.tick_params(labelsize=16)

    # chemical species legend
    species_handles = [Line2D([0], [0], color=colors[i], linewidth=1.2) for i in range(len(species))]
    leg1 = ax.legend(species_handles, species, loc = 'lower left', fontsize = 14, frameon = False)
    ax.add_artist(leg1)

    # model legend
    model_handles = [Line2D([0], [0], color='k', linewidth = 1.2, ls='-'), Line2D([0], [0], color='k', linewidth = 1.2, ls='--')]

    ax.legend(model_handles, ['ATMO', 'VULCAN'], loc = 'lower left', bbox_to_anchor=(0, 0.30), fontsize = 14, frameon = False)

    plt.tight_layout()
    output = tit.replace(" ", "_") + ".png"
    plt.savefig(output, dpi=300, bbox_inches='tight')
    plt.close()
    print('Created:', output)
