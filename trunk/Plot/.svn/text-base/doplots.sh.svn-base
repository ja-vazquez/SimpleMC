#!/bin/bash

ls -1 *params > input.txt
sed -i 's/.params//g' input.txt


while read i
do
    rm "${i}".paramnames
    while read j
    do
       echo "$j $j" >> "${i}".paramnames
    done < "${i}".params
done < input.txt

while read i
do
    cp distparams_base.ini distparams.ini
    sed -i "s/base/$i/g" distparams.ini
    ./getdist distparams.ini
    python $i.py
done < input.txt

