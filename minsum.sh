#!/bin/bash

for file in $(ls ../Results/*hapcut)
do
python min_sum.py -i $file -o ../Results/ -m
done
