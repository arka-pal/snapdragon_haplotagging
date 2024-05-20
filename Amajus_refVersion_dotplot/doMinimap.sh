#!/bin/bash

t=$1
q=$2

minimap2 -PD -k19 -w19 -m200 -t48 $q $t > ${t}_${q}.paf