#!/usr/bin/env python
# coding: utf-8
import numpy as np

from matplotlib.pyplot import cm, colorbar, tight_layout
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import ImageGrid
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['figure.figsize'] = ((3+3/8), (3+3/8))
plt.rc('text.latex', preamble=r'\usepackage{bm}')

import matplotlib.pyplot as plt
import pandas as pd

import os

find_str = '*.csv'
path_in, indicate_bond = './output/all_optimizations', True
# path_in, indicate_bond = './output/once_per_sweep', False
path_out = './figs'
fns = []
for fn in os.popen("find " + str(path_in) + " -path "
                        + '"' + str(find_str) + '"').read().split('\n')[0:-1]:
    fns.append(fn)
dfs = [pd.read_csv(fn) for fn in fns]
fig, ax = plt.subplots()
for (df,fn) in zip(dfs,fns):
    bonds = np.sort(np.unique(df['bond']))
    sweeps = np.sort(np.unique(df['sweep']))
    half_sweeps = np.sort(np.unique(df['half_sweep']))

    s,hs,b = sweeps[0],half_sweeps[0],bonds[0]
    mx = 1.175*np.max(np.abs(df['X']))
    data = df[df['sweep']==s]
    data = data[data['half_sweep']==hs]
    data = data[data['bond']==b]
    if len(data.iloc[:]) < 1:
        continue
    energies = np.unique(data['E'])
    mz = data['Sz'].mean()
    im = ax.scatter(data['X'], data['Y'], c=data['Sz'], marker='h', s=600, cmap='RdBu_r', vmin=-0.5, vmax=0.5)
    if indicate_bond:
        ax.scatter(data['X'].iloc[0], data['Y'].iloc[0], c='black', s=600, marker='h', cmap='RdBu_r', vmin=-0.5, vmax=0.5, alpha=0.6)
    ax.quiver(data['X'], data['Y'], data['Sx'], data['Sy'], units="xy", width=0.07, scale=1, pivot="middle", color="white")
    # plt.colorbar(im)
    ax.set_xlim(-mx,mx)
    ax.set_ylim(-mx,mx)
    s_str = f'{s}'.zfill(4)
    sweep_step = 0
    b_str = f'{sweep_step}'.zfill(3)
    ax.axis('off')
    ene_str = '{:.12f}'.format(energies[0]).zfill(16)
    mz_str = '\overline{m}_z'
    mz_val = '{:.12f}'.format(mz)
    ax.text(0,4.25,f'$E={ene_str}$',ha='center',va='center')
    ax.text(0,-4.35,f'${mz_str}={mz_val}$',ha='center',va='center')
    fn_rpl = fn.replace('.csv','')
    if not os.path.isdir(f'{path_out}/{fn_rpl}/'):
        os.makedirs(f'{path_out}/{fn_rpl}/')
    plt.tight_layout()
    plt.savefig(f'{path_out}/{fn_rpl}/{s_str}_{b_str}.png',dpi=600)
    plt.cla()

    bonds = bonds[1:]
    sweeps = sweeps[1:]
    half_sweeps = half_sweeps[1:]
    mx = 1.175*np.max(np.abs(df['X']))
    for s in sweeps:
        for hs in half_sweeps:
            for b in bonds:
                print(s, hs, b, sweep_step)
                s_str = f'{s}'.zfill(4)
                fn_rpl = fn.replace('.csv','')
                sweep_step = int((hs-1)*121+((-1)**(hs-1)*b))
                b_str = f'{sweep_step}'.zfill(3)
                if os.path.exists(f'{path_out}/{fn_rpl}/{s_str}_{b_str}.png'):
                    continue
                data = df[df['sweep']==s]
                data = data[data['half_sweep']==hs]
                data = data[data['bond']==b]
                if len(data.iloc[:]) < 1:
                    continue
                energies = np.unique(data['E'])
                mz = data['Sz'].mean()
                im = ax.scatter(data['X'], data['Y'], c=data['Sz'], marker='h', s=600, cmap='RdBu_r', vmin=-0.5, vmax=0.5)
                if indicate_bond:
                    ax.scatter(data['X'].iloc[b-(hs-1)], data['Y'].iloc[b-(hs-1)], c='black', s=600, marker='h', cmap='RdBu_r', vmin=-0.5, vmax=0.5, alpha=0.6)
                ax.quiver(data['X'], data['Y'], data['Sx'], data['Sy'], units="xy", width=0.07, scale=1, pivot="middle", color="white")
                # plt.colorbar(im)
                ax.set_xlim(-mx,mx)
                ax.set_ylim(-mx,mx)
                ax.axis('off')
                ene_str = '{:.12f}'.format(energies[0]).zfill(16)
                mz_str = '\overline{m}_z'
                mz_val = '{:.12f}'.format(mz)
                ax.text(0,4.25,f'$E={ene_str}$',ha='center',va='center')
                ax.text(0,-4.35,f'${mz_str}={mz_val}$',ha='center',va='center')
                if not os.path.isdir(f'{path_out}/{fn_rpl}/'):
                    os.makedirs(f'{path_out}/{fn_rpl}/')
                plt.tight_layout()
                plt.savefig(f'{path_out}/{fn_rpl}/{s_str}_{b_str}.png',dpi=1200)
                plt.cla()