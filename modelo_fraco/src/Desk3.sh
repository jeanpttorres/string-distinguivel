#!/bin/bash

DATA=data/instancias-tcc/grupo1/*.dssp
make clean
#make READLINE=false ZIMPL=false LPS=cpx EXECUTE=ONLYROOT PROBLEM=SUB HEUR=RA TYPE=CSP
make READLINE=false ZIMPL=false EXECUTE=ONLYROOT PROBLEM=SUB HEUR=VNS RA_TYPE=ROUNDING TYPE=CSP 
echo 'frac;inteiras;filename;nodes;lb;ub;root;time;status;heuristica;valor' > VNS-SUB-1-diff.csv
for i in $DATA; do
    ./bin/dssp-NOTRACE-ONLYROOT-SUB-VNS-CSP-BASIC-ROUNDING-VELHO -f $i -q >> VNS-SUB-1-diff.csv
    echo $i
done
