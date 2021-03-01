#!/bin/bash
i=256
j=4096
while [ $i -le $j ]
do
  time -o ./tempos/omp_$i.txt ./mat-mul $i $i  0
  i=$(expr $i \* 2)
done
 