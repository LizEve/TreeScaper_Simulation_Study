#!/bin/bash

cd /Users/chatnoir/Projects/TreeScaper_ExampleData/Parks_2017_gene_trees/


# Move all trees that have all species 
mkdir occ0.8sp100
cd trees.bootstrap._1.inclade1.ortho1.occ_0.8
for f in $(cat tree.summaries.occ_0.8.txt.100.txt.bootstrap.txt);
do
cp "$f" ../occ0.8sp100;
done

# Concat all 
cd ../occ0.8sp100
cat *.tre >> all.0.8.tre

# subsample 
head -n 10 *p.tre >> all.0.2.sub10.tre
head -n 10 *p.tre >> all.0.5.sub10.tre
head -n 10 *p.tre >> all.0.8.sub10.tre


# In treescaper online 
[1-320: blue]; [321-640:red]; [641-960:green]