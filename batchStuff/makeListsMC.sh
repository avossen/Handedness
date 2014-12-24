#!/bin/bash

for ex in 07 09 11 13 15 17 19 21 23 25 27 31 33 35 37 39 41 43 45 47 49 51 53 55 61 63 65 67 69 71 73
do
for res in on_resonance continuum
do
for spec in uds charm mixed charged
do
find /pic/projects/Belle/ZMDST/MC/MC_4S/MC_4S_EXP$ex/MC_4S_EXP$ex/$res/$spec/ -iname '*.mdst'> lists/mc$ex\_$res\_$spec.list
done
done
done

#/pic/projects/Belle/ZMDST/MC/MC_4S/MC_4S_EXP55/MC_4S_EXP55/continuum/charm/