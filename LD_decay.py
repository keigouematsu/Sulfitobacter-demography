import numpy as np
import pandas as pd
import itertools as iter
import matplotlib as mpl
import matplotlib.pyplot as plt
from decimal import Decimal
import os, glob

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

data2=[ data_002, data_005, data_01, data_02, data_03]
color2=['b','c','k','y','m']
label2=['G=0.02','G=0.05','G=0.1','G=0.2','G=0.3']

for j in range(5):
    dtm=data2[j]
    col=color2[j]
    lab=label2[j]
    dtm['dist10000']= dtm[1].values*10000
    dtm.columns=['tr_comp','dist2','dist']
    print(dtm)
    width2_20 = np.linspace(0, 1, 501)
    dtm['key_100'] = pd.cut(dtm.dist2, width2_20)
    dtm_100 = dtm.groupby('key_100').mean()
    #plt.scatter(x='dist', y='tr_comp', color=col,data= dtm_100,s=10,marker='v', alpha=0.5,label=lab)
    cubic = interpolate.interp1d(dtm_100['dist'], dtm_100['tr_comp'],kind='cubic')
    p0 = (1, 1, 0.38) # start with values near those we expect
    params, cv = scipy.optimize.curve_fit(monoExp, dtm_100['dist'], dtm_100['tr_comp'],p0)
    m, t, b = params
    print(params)
    x = np.linspace(5, 9995, num=10001, endpoint=True)
    #plt.plot(x, cubic(x),color=col,label=lab)
    plt.plot(x, monoExp(x, m, t, b), '--',linewidth = 0.5, alpha= 1,color=col,label=label2[j])


dist_genome = data_genome[1].values
value_genome = data_genome[0].values


width_100 = range(0,10000,100)
width_50 = range(0,10000,50)
width_10 = range(0,10000,10)
width_4 = range(0,10000,4)
width_1 = range(0,10000,1)

#data=[data_fs, data_third,data_genome]
#color=['r','b','k']
#label=['1st+2nd','3rd','genome']
data=[data_genome]
color=['r']
label=['Observed']
for i in range(1):
    dt=data[i]
    col=color[i]
    lab=label[i]
    dt.columns=['tr_comp','dist']
    print(dt)

    dt['key_100'] = pd.cut(dt.dist, width_100)
    dt['key_50'] = pd.cut(dt.dist, width_50)
    dt['key_10'] = pd.cut(dt.dist, width_10)
    dt['key_4'] = pd.cut(dt.dist, width_4)
    dt['key_1'] = pd.cut(dt.dist, width_1)

    dt_100 = dt.groupby('key_100').mean()
    dt_100_size = dt.groupby('key_100').size()
    dt_50 = dt.groupby('key_50').mean()
    dt_50_size = dt.groupby('key_50').size()
    dt_10 = dt.groupby('key_10').mean()
    dt_10_size = dt.groupby('key_10').size()
    dt_4 = dt.groupby('key_4').mean()
    dt_4_size = dt.groupby('key_100').size()
    dt_1 = dt.groupby('key_1').mean()
    dt_size = dt.groupby('key_1').size()

    plt.scatter(x='dist', y='tr_comp', s=30,color=col,data= dt_50,label=lab,alpha=0.9, marker="o")
plt.legend(bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0, fontsize=8)
plt.xlabel("Distance")
plt.ylabel("LD ($r^{2}$)")

plt.xlim(-5, 1000)
plt.ylim(0.1,0.7)
plt.hlines([0.335], 0, 10000, "black", linestyles='dashed') 
plt.savefig('LD_decay_50bp.pdf')
plt.show()