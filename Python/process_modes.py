#!/usr/bin/env python
import os
import sys

fn = "low_modes.dat"
fin = open(fn,'r')
lines = [l.split() for l in fin.readlines() ]
Sets = [set([round(float(x),3) for x in line]) for line in lines]
x = Sets[0]
for s in Sets:
    print s
    x = x.intersection(s)
print x


