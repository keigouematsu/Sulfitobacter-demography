#!/bin/sh

find ./iqtree_ms_recomb_analysis -iname '*.treefile' -type f -print | xargs gsed -i -E 's/[0-9]+\.[0-9]+//g'

var="./iqtree_ms_recomb_analysis/*.treefile"
for i in ${var}
do
  bar=`gsed 's/^.*)\([0-9]\|[0-9][0-9]\|100\):.*$/\1/' $i`
  if [ $bar -lt 60 ] ; then 
    echo `basename $i` $bar 
    echo `basename $i` $bar >> ./iqtree_ms_recomb_analysis/excluded.txt
  elif grep -E '(c1:,c2:)|(c2:,c1:)' $i; then 
    echo `basename $i` $bar >> ./iqtree_ms_recomb_analysis/clade3out.txt
  elif grep -E '(c1:,c3:)|(c3:,c1:)' $i; then 
    echo `basename $i` $bar >> ./iqtree_ms_recomb_analysis/clade2out.txt
  elif grep -E '(c2:,c3:)|(c3:,c2:)' $i; then 
    echo `basename $i` $bar >> ./iqtree_ms_recomb_analysis/clade1out.txt
  else
    echo "wrong"
  fi
done

echo "(c1,(c2,c3))"
cat ./iqtree_ms_recomb_analysis/clade1out.txt | wc -l 
echo "(c2,(c1,c3))"
cat ./iqtree_ms_recomb_analysis/clade2out.txt | wc -l 
echo "(c3,(c1,c2))"
cat ./iqtree_ms_recomb_analysis/clade3out.txt | wc -l 
echo "unknown"
cat ./iqtree_ms_recomb_analysis/excluded.txt | wc -l 
