#!/bin/sh


#estimating demography with dadi
python ./dadi_sulfito.py


#no recombination 967bp 
cd ./msPro
mkdir norecomb
./msPro 3 10000 967 -t 39.9371 -I 3 1 1 1  -n 1 2.18729e-05 -n 2 0.00333254 -n 3 1.915e-05 -en 0.0123134 2 0.000149668 -en 0.0123177 1 1 -ej 0.0123134 3 2 -ej 0.0123177 2 1 -c 0 10000 -b 0 dist.txt > ./norecomb/msPro.out
cd ./norecomb
gsed -i -e '1,3d' msPro.out
gsed -i -e '/segsites: 0/,+3d' msPro.out
split -a 4  -l 7 msPro.out
find . -name "x*" -print0 | xargs -0 gsed -i -e '1,2d'
find . -name "x*" -print0 | xargs -0 gsed -i -e 's/positions: //g'

#recombination G=0.5 967bp  tract length=10000
cd ../
mkdir recomb
./msPro 3 10000 967 -t 39.9371 -I 3 1 1 1  -n 1 2.18729e-05 -n 2 0.00333254 -n 3 1.915e-05 -en 0.0123134 2 0.000149668 -en 0.0123177 1 1 -ej 0.0123134 3 2 -ej 0.0123177 2 1 -c 483.5 10000 -b 0 dist.txt > ./recomb/msPro.out
cd ./recomb
gsed -i -e '1,3d' ./msPro.out
gsed -i -e '/segsites: 0/,+3d' ./msPro.out
split -a 4  -l 7 ./msPro.out
find . -name "x*" -print0 | xargs -0 gsed -i -e '1,2d'
find . -name "x*" -print0 | xargs -0 gsed -i -e 's/positions: //g'
cd ../

#making phylip file
cd ../
python ms_to_phy_norecomb.py
python ms_to_phy_recomb.py


#iqtree bootstrap from simulated data without recombination
rm -r ./iqtree_ms_norecomb
mkdir ./iqtree_ms_norecomb
cd ./iqtree_ms_norecomb
find ../msPro/norecomb -iname '*.phy' | parallel -j 8 'iqtree -s {}  -nt 1 -m MFP  -mrate E,I,G,I+G -nstop 500 -pers 0.5 -b 100 -pre {/.}'
cd ../
rm -r ./iqtree_ms_norecomb_analysis
mkdir ./iqtree_ms_norecomb_analysis
find ./iqtree_ms_norecomb -name "*.treefile" | xargs -J % cp % ./iqtree_ms_norecomb_analysis

sh ./iqtree_bootstrap_norecomb_analysis.sh

#iqtree bootstrap from simulated data with recombination
rm -r ./iqtree_ms_recomb
mkdir ./iqtree_ms_recomb
cd ./iqtree_ms_recomb
find ../msPro/recomb -iname '*.phy' | parallel -j 8 'iqtree -s {}  -nt 1 -m MFP  -mrate E,I,G,I+G -nstop 500 -pers 0.5 -b 100 -pre {/.}'
cd ../
rm -r ./iqtree_ms_recomb_analysis
mkdir ./iqtree_ms_recomb_analysis
find ./iqtree_ms_recomb -name "*.treefile" | xargs -J % cp % ./iqtree_ms_recomb_analysis

sh ./iqtree_bootstrap_analysis.sh


#iqtree bootstrap from real Sulfitobacter gene sequences
rm -r ./iqtree_realdata
mkdir ./iqtree_realdata
cd ./iqtree_realdata
find ../bacteria_phylip -iname '*.phy' | parallel -j 8 'iqtree -s {} -o out -nt 1 -m MFP  -mrate E,I,G,I+G -nstop 500 -pers 0.5 -b 100 -pre {/.}'
cd ../
rm -r ./iqtree_realdata_analysis
mkdir ./iqtree_realdata_analysis
find ./iqtree_realdata -name "*.treefile" | xargs -J % cp % ./iqtree_realdata_analysis

sh ./iqtree_bootstrap_realdata_analysis.sh




