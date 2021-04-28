import numpy as np
import pandas as pd
import os, glob
import random as rd

folder_path = './msPro/recomb'
file_paths = glob.glob(os.path.join(folder_path, 'x[a-z][a-z][a-z][a-z]'))

atgc = ('A','T','G','C')
def bintoseq(inseq):
    bases = rd.sample(atgc,2)
    x=bases[0]
    y=bases[1]
    inseq2 = inseq.replace({'0':x,'1':y})
    return inseq2

for file_path in file_paths:

    with open(file_path) as f: 
    	data = f.readlines()
    	l = str(len(data[1]) - 1)
    	out= list(str(0)* int(l))
    	pop1 = list(data[0].strip()) #remove newline codes
    	pop2 = list(data[1].strip())
    	pop3 = list(data[2].strip())
    	tabl = pd.DataFrame({'out':out, 'pop1': pop1,'pop2': pop2,'pop3': pop3}).astype('str')
    	tabl2 = tabl.apply(bintoseq,axis=1)
    	tabl3 = tabl2.astype(str).apply(''.join)
    with open(os.path.splitext(file_path)[0] + '.rep', 'a') as f:
        	f.write('\n'.join(tabl3))

file_paths_reps = glob.glob(os.path.join(folder_path, 'x[a-z][a-z][a-z][a-z].rep'))
for file_path in file_paths_reps:
   
    with open(file_path) as f: 
        data = f.readlines()
        data = data[0:5]
        l = str(len(data[1]) - 1)
    # to phylip file
        data.insert(0, 'out\t')
        data.insert(2, 'c3\t')
        data.insert(4, 'c1\t')
        data.insert(6, 'c2\t')
        data.insert(0,'4 '+ l + '\n')

    with open(os.path.splitext(file_path)[0] + '.phy', 'a') as f:
        f.writelines(data)
        