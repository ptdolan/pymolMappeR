#!/bin/bash
for i in *.pdb;
do
echo $i
echo ${i/pdb/out}
mkdssp -i $i -o ${i/\.pdb/\.out}
done

