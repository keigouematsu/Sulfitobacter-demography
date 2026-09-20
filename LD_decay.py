import numpy as np
import pandas as pd
import itertools as iter
import matplotlib as mpl
import matplotlib.pyplot as plt
from decimal import Decimal
import os, glob
from fractions import Fraction

# fit with ms simulated data
from scipy import interpolate
import scipy.optimize

data_genome = pd.read_table('./sulfito_ld_biallel_rev.txt',header=None)
data_002 = pd.read_table('./msPro/recomb_002/merge_tree.txt',header=None)
data_005 = pd.read_table('./msPro/recomb_005/merge_tree.txt',header=None)
data_01 = pd.read_table('./msPro/recomb_01/merge_tree.txt',header=None)
data_02 = pd.read_table('./msPro/recomb_02/merge_tree.txt',header=None)
data_03 = pd.read_table('./msPro/recomb_03/merge_tree.txt',header=None)


def monoExp(x, m, t, b):

    return  m * np.exp((x/1000) / (t+x/1000)  )+ b

data2=[ data_12, data_2, data_3, data_32, data_4]
color2=['b','c','k','y','m']
label2=['G=0.02','G=0.05','G=0.1','G=0.2','G=0.3']

for j in range(5):
    dtm=data2[j]
    col=color2[j]
    lab=label2[j]
    dtm['dist10000']= dtm[1].values*10000
    dtm.columns=['tr_comp','dist2','dist']
    width2_20 = np.linspace(0, 1, 501)
    dtm['key_100'] = pd.cut(dtm.dist2, width2_20)
    dtm_100 = dtm.groupby('key_100').mean()
    cubic = interpolate.interp1d(dtm_100['dist'], dtm_100['tr_comp'],kind='cubic')
    p0 = (1, 1, 0.38) # start with values near those we expect
    params, cv = scipy.optimize.curve_fit(monoExp, dtm_100['dist'], dtm_100['tr_comp'],p0)
    m, t, b = params
    print(params)
    x = np.linspace(5, 9995, num=10001, endpoint=True)
    plt.plot(x, monoExp(x, m, t, b), '--',linewidth = 0.5, alpha= 1,color=col,label=label2[j])


dist_genome = data_genome[1].values
value_genome = data_genome[0].values

width_50 = range(0,10000,50)

data=[data_genome]
color=['r']
label=['Observed']
for i in range(1):
    dt=data[i]
    col=color[i]
    lab=label[i]
    dt.columns=['tr_comp','dist']
    print(dt)
    dt['key_50'] = pd.cut(dt.dist, width_50)
    dt_50 = dt.groupby('key_50').mean()
    dt_50_size = dt.groupby('key_50').size()
    plt.scatter(x='dist', y='tr_comp', s=30,color=col,data= dt_50,label=lab,alpha=0.9, marker="o")
    
plt.legend(bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0, fontsize=8)
plt.xlabel("Distance")
plt.ylabel("LD ($r^{2}$)")

plt.xlim(-5, 1000)
plt.ylim(0.1,0.7)
plt.hlines([Fraction(1, 3)], 0, 10000, "black", linestyles='dashed') 
plt.savefig('LD_decay_50bp.pdf')
plt.show()