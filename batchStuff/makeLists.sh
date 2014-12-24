#!/bin/bash

for ex in 07 09 11 13 15 17 19 21 23 25 27 31 33 35 37 39 41 43 45 47 49 51 53 55 61 63 65 67 69 71 73
do
for res in on_resonance continuum
do
find /pic/projects/Belle/ZMDST/DATA/SKIM/HADRON_BJ/EXP$ex/$res/ -iname '*.mdst' > data$ex\_$res.list
done
done