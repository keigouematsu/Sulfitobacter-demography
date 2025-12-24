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
mkdir recomb_01
./msPro 4 10000 -t 39.856  -I 4 1 1 1 1  -n 1 0.00432333 -n 2 0.000177381 -n 3 0.00532747 -n 4 0.00981703 -eN 0.01629025 0.00845891 -eN 0.01693564 0.00479238 -eN 0.017024354 1 -ej 0.01629025 4 3 -ej 0.01693564 3 2 -ej 0.017024354 2 1 -c 96.7 1000 -b 0 dist.txt > ./recomb_01/msPro.out
cd ./recomb_01
gsed -i -e '1,3d' ./msPro.out
gsed -i -e '/segsites: 0/,+3d' ./msPro.out
split -a 4  -l 7 ./msPro.out
find . -name "x*" -print0 | xargs -0 gsed -i -e '1,2d'
find . -name "x*" -print0 | xargs -0 gsed -i -e 's/positions: //g'
cd ../

#recombination G=0.1 967bp  tract length=1000
cd ../
mkdir recomb_01
./msPro 4 10000 -t 39.856  -I 4 1 1 1 1  -n 1 0.00432333 -n 2 0.000177381 -n 3 0.00532747 -n 4 0.00981703 -eN 0.01629025 0.00845891 -eN 0.01693564 0.00479238 -eN 0.017024354 1 -ej 0.01629025 4 3 -ej 0.01693564 3 2 -ej 0.017024354 2 1 -c 96.7 1000 -b 0 dist.txt > ./recomb_01/msPro.out
cd ./recomb_01
gsed -i -e '1,3d' ./msPro.out
gsed -i -e '/segsites: 0/,+3d' ./msPro.out
split -a 4  -l 7 ./msPro.out
find . -name "x*" -print0 | xargs -0 gsed -i -e '1,2d'
find . -name "x*" -print0 | xargs -0 gsed -i -e 's/positions: //g'
cd ../

#recombination G=0.02 967bp  tract length=1000
cd ../
mkdir recomb_002
./msPro 4 10000 -t 39.856  -I 4 1 1 1 1  -n 1 0.00432333 -n 2 0.000177381 -n 3 0.00532747 -n 4 0.00981703 -eN 0.01629025 0.00845891 -eN 0.01693564 0.00479238 -eN 0.017024354 1 -ej 0.01629025 4 3 -ej 0.01693564 3 2 -ej 0.017024354 2 1 -c 19.34 1000 -b 0 dist.txt > ./recomb_002/msPro.out
cd ./recomb_002
gsed -i -e '1,3d' ./msPro.out
gsed -i -e '/segsites: 0/,+3d' ./msPro.out
split -a 4  -l 7 ./msPro.out
find . -name "x*" -print0 | xargs -0 gsed -i -e '1,2d'
find . -name "x*" -print0 | xargs -0 gsed -i -e 's/positions: //g'
cd ../

#recombination G=0.05 967bp  tract length=1000
cd ../
mkdir recomb_005
./msPro 4 10000 -t 39.856  -I 4 1 1 1 1  -n 1 0.00432333 -n 2 0.000177381 -n 3 0.00532747 -n 4 0.00981703 -eN 0.01629025 0.00845891 -eN 0.01693564 0.00479238 -eN 0.017024354 1 -ej 0.01629025 4 3 -ej 0.01693564 3 2 -ej 0.017024354 2 1 -c 48.35 1000 -b 0 dist.txt > ./recomb_005/msPro.out
cd ./recomb_005
gsed -i -e '1,3d' ./msPro.out
gsed -i -e '/segsites: 0/,+3d' ./msPro.out
split -a 4  -l 7 ./msPro.out
find . -name "x*" -print0 | xargs -0 gsed -i -e '1,2d'
find . -name "x*" -print0 | xargs -0 gsed -i -e 's/positions: //g'
cd ../

#recombination G=0.2 967bp  tract length=1000
cd ../
mkdir recomb_02
./msPro 4 10000 -t 39.856  -I 4 1 1 1 1  -n 1 0.00432333 -n 2 0.000177381 -n 3 0.00532747 -n 4 0.00981703 -eN 0.01629025 0.00845891 -eN 0.01693564 0.00479238 -eN 0.017024354 1 -ej 0.01629025 4 3 -ej 0.01693564 3 2 -ej 0.017024354 2 1 -c 193.4 1000 -b 0 dist.txt > ./recomb_02/msPro.out
cd ./recomb_02
gsed -i -e '1,3d' ./msPro.out
gsed -i -e '/segsites: 0/,+3d' ./msPro.out
split -a 4  -l 7 ./msPro.out
find . -name "x*" -print0 | xargs -0 gsed -i -e '1,2d'
find . -name "x*" -print0 | xargs -0 gsed -i -e 's/positions: //g'
cd ../

#recombination G=0.3 967bp  tract length=1000
cd ../
mkdir recomb_03
./msPro 4 10000 -t 39.856  -I 4 1 1 1 1  -n 1 0.00432333 -n 2 0.000177381 -n 3 0.00532747 -n 4 0.00981703 -eN 0.01629025 0.00845891 -eN 0.01693564 0.00479238 -eN 0.017024354 1 -ej 0.01629025 4 3 -ej 0.01693564 3 2 -ej 0.017024354 2 1 -c 290.1 1000 -b 0 dist.txt > ./recomb_03/msPro.out
cd ./recomb_03
gsed -i -e '1,3d' ./msPro.out
gsed -i -e '/segsites: 0/,+3d' ./msPro.out
split -a 4  -l 7 ./msPro.out
find . -name "x*" -print0 | xargs -0 gsed -i -e '1,2d'
find . -name "x*" -print0 | xargs -0 gsed -i -e 's/positions: //g'
cd ../

python ./ms_to_r2_norecomb.py
python ./ms_to_r2_recomb.py

find ./norecomb -name "*_LD_r2_ms.txt" -print0 | xargs -0 cat > ./norecomb/merge_tree.txt
find ./recomb_002 -name "*_LD_r2_ms.txt" -print0 | xargs -0 cat > ./recomb_002/merge_tree.txt
find ./recomb_005 -name "*_LD_r2_ms.txt" -print0 | xargs -0 cat > ./recomb_005/merge_tree.txt
find ./recomb_01 -name "*_LD_r2_ms.txt" -print0 | xargs -0 cat > ./recomb_01/merge_tree.txt
find ./recomb_02 -name "*_LD_r2_ms.txt" -print0 | xargs -0 cat > ./recomb_02/merge_tree.txt
find ./recomb_03 -name "*_LD_r2_ms.txt" -print0 | xargs -0 cat > ./recomb_03/merge_tree.txt

python ./LD_decay.py



