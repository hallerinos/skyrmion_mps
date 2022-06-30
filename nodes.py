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
path_in = './output/slow'
path_out = './figs'

indicate_bond = True
fns = []
for fn in os.popen("find " + str(path_in) + " -path "
                        + '"' + str(find_str) + '"').read().split('\n')[0:-1]:
    fns.append(fn)
fns = fns[0:1]
dfs = [pd.read_csv(fn) for fn in fns]
fig, ax = plt.subplots()
for (df,fn) in zip(dfs,fns):
    # ax.scatter(df['X'], df['Y'], c='white', marker='h', edgecolors='black', s=600)
    energy = np.unique(df['E'])[0]
    df = df[df['E']==energy]
    for (id,(x,y)) in enumerate(zip(df['X'], df['Y'])):
        print(id)
        ax.scatter(x, y, c='white', marker='h', edgecolors='black', s=600)
        # ax.text(x+0.04,y-0.02,f'${id+1}$',va='center',ha='center')
    # plt.colorbar(im)
    mx = 1.175*np.max(np.abs(df['X']))
    ax.set_xlim(-mx,mx)
    ax.set_ylim(-mx,mx)
    fn_rpl = fn.replace('.csv','')
    if not os.path.isdir(f'{path_out}/{fn_rpl}/'):
        os.makedirs(f'{path_out}/{fn_rpl}/')
    ax.axis('off')
    plt.tight_layout()
    plt.savefig(f'{path_out}/{fn_rpl}/nodes.png',dpi=300)
    plt.cla()