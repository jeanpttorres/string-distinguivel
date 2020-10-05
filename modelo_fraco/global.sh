#!/bin/bash

DATA=data/instancias-tcc/grupo1/*.dssp
make clean
make READLINE=false ZIMPL=false EXECUTE=ONLYROOT PROBLEM=SUB HEUR=RA LPS=cpx TYPE=DSP
echo 'frac;inteiras;filename;nodes;lb;ub;root;time;status;heuristica;valor' > G1-RA.csv
for i in $DATA; do
    ./bin/dssp-NOTRACE-ONLYROOT-SUB-RA-DSP-BASIC-ROUNDING-VELHO -f $i -q >> G1-RA.csv
    echo $i
done
