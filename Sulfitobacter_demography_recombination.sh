#!/bin/sh


#estimating demography with dadi
python ./dadi_sulfito.py


#no recombination 967bp 
cd ./msPro
mkdir norecomb
./msPro 4 10000 -t 39.856  -I 4 1 1 1 1  -n 1 0.00432333 -n 2 0.000177381 -n 3 0.00532747 -n 4 0.00981703 -eN 0.01629025 0.00845891 -eN 0.01693564 0.00479238 -eN 0.017024354 1 -ej 0.01629025 4 3 -ej 0.01693564 3 2 -ej 0.017024354 2 1 -c 0 1000 -b 0 dist.txt > ./norecomb/msPro.out
cd ./norecomb
gsed -i -e '1,3d' msPro.out
gsed -i -e '/segsites: 0/,+3d' msPro.out
split -a 4  -l 7 msPro.out
find . -name "x*" -print0 | xargs -0 gsed -i -e '1,2d'
find . -name "x*" -print0 | xargs -0 gsed -i -e 's/positions: //g'

#recombination G=0.1 967bp  tract length=1000
cd ../
mkdir recomb
./msPro 4 10000 -t 39.856  -I 4 1 1 1 1  -n 1 0.00432333 -n 2 0.000177381 -n 3 0.00532747 -n 4 0.00981703 -eN 0.01629025 0.00845891 -eN 0.01693564 0.00479238 -eN 0.017024354 1 -ej 0.01629025 4 3 -ej 0.01693564 3 2 -ej 0.017024354 2 1 -c 96.7 1000 -b 0 dist.txt > ./recomb/msPro.out
cd ./recomb
gsed -i -e '1,3d' ./msPro.out
gsed -i -e '/segsites: 0/,+3d' ./msPro.out
split -a 4  -l 7 ./msPro.out
find . -name "x*" -print0 | xargs -0 gsed -i -e '1,2d'
find . -name "x*" -print0 | xargs -0 gsed -i -e 's/positions: //g'
cd ../


python ./ms_to_r2_norecomb.py
python ./ms_to_r2_recomb.py

find ./norecomb -name "*_LD_r2_ms.txt" -print0 | xargs -0 cat > ./norecomb/merge_tree.txt
find ./recomb -name "*_LD_r2_ms.txt" -print0 | xargs -0 cat > ./recomb/merge_tree.txt

python ./LD_decay.py
python ./LD_decay_difftractlength.py



