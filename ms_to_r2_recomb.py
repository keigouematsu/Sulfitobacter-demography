import numpy as np
import pandas as pd
import os, glob
import random as rd

folder_path = './msPro/recomb'
file_paths = glob.glob(os.path.join(folder_path, 'x[a-z][a-z][a-z][a-z]'))

for j in range(1):
    rate=rates[j]

    for file_path in file_paths:
        with open(file_path) as f: 
            data = f.readlines()
            dist = data[0].split()
            c12 = list(data[1].strip()) #改行コードを除去するためにstrip()を私用
            c9 = list(data[2].strip())
            c1 = list(data[3].strip())
            c2 = list(data[4].strip())
            tabl = pd.DataFrame({'dist': dist,'clade12': c12,
                    'clade9': c9,
                    'clade2': c2,'clade1': c1},dtype='float').T
            #tabl = tabl.loc[:, tabl[1:4].apply(sum)>=2] # 1を２個以上含む
            tabl.columns = range(tabl.shape[1])
            
        for i, j in list(iter.combinations(range(tabl.shape[1]), 2)):
    
            tr3 =tabl.iloc[1:5][[i, j]].astype(int).values
            
            z1=((tr3[:,0]==0)).sum()/4
            #print(z1)
            z2=((tr3[:,1]==0)).sum()/4
            #print(z2)
            o1=((tr3[:,0]==1)).sum()/4
            #print(o1)
            o2=((tr3[:,1]==1)).sum()/4
            #print(o2)
        
            zz = (tr3 == np.array([0,0])).all(axis=1).sum()/4 
            zo = (tr3 == np.array([0,1])).all(axis=1).sum()/4 
            oz = (tr3 == np.array([1,0])).all(axis=1).sum()/4 
            oo = (tr3 == np.array([1,1])).all(axis=1).sum()/4 

            r2 = ((zz*oo-zo*oz)**2)/(z1*z2*o1*o2)
            #print(r2)

            dt=abs(Decimal(tabl[i][0])-Decimal(tabl[j][0]))
            with open(os.path.splitext(file_path)[0] + '_LD_r2_ms.txt', 'a') as f:
                f.write(str(r2) + '\t' + str(dt)+ '\n')