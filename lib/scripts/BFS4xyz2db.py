#! /usr/bin/env python

# Fragment detection algorithm
# This is an iterative BFS
# Michael Marshall
#
# simply, start with the first atom, do a BFS
# when done, go to any untouched atom and start again
# iterate until all atoms belong to a fragment group

# BFS algorithm is a simple 3 colour scheme
# white (untouched)
# grey (Queued)
# black (touched and all edges discovered)

import os
import sys
import string
import math


def convert(p):
    if symbol[p] == 'H':
        d = 1.001
    if symbol[p] == 'He':
        d = 1.012
    if symbol[p] == 'Li':
        d = 0.825
    if symbol[p] == 'Be':
        d = 1.408
    if symbol[p] == 'B':
        d = 1.485
    if symbol[p] == 'C':
        d = 1.452
    if symbol[p] == 'N':
        d = 1.397
    if symbol[p] == 'O':
        d = 1.342
    if symbol[p] == 'F':
        d = 1.287
    if symbol[p] == 'Ne':
        d = 1.243
    if symbol[p] == 'Na':
        d = 1.144
    if symbol[p] == 'Mg':
        d = 1.364
    if symbol[p] == 'Al':
        d = 1.639
    if symbol[p] == 'Si':
        d = 1.716
    if symbol[p] == 'P':
        d = 1.705
    if symbol[p] == 'S':
        d = 1.683
    if symbol[p] == 'Cl':
        d = 1.639
    if symbol[p] == 'Ar':
        d = 1.595

    return d / 1.5


if len(sys.argv) > 1:
    F = open(sys.argv[1])
else:
    print 'Need a filename!'
    exit(1)

#Angs   H    C    N    B     O     F    P     S   Cl
VdW = [1.2, 1.7, 1.5, 1.55, 1.52, 1.9, 1.85, 1.8]
#VdW = [1.09,1.70,2.0,1.52,1.47,1.80,1.80,1.75]
#VdW = [0.23,0.68,0.68,0.68,0.64,1.05,1.02,0.99] # covalent radius

numatoms = int(F.readline())
LABEL = F.readline()

symbol = range(numatoms)
X = [0.0] * numatoms
Y = [0.0] * numatoms
Z = [0.0] * numatoms

Queue = []
White = []
Black = []

for f in range(0, numatoms):
    A = F.readline().split()
    symbol[f] = A[0]
    X[f] = float(A[1])
    Y[f] = float(A[2])
    Z[f] = float(A[3])
    White.append(f)

Fragment = [[] for i in range(numatoms)]  # stores fragments
#Max number of fragments has to be numatoms

start = 0  # starts with the first atom in the list
Queue.append(start)
White.remove(start)

frag = 0

while((len(White) > 0) or (len(Queue) > 0)):  # Iterates to the next fragment
    while(len(Queue) > 0):  # BFS within a fragment
        for u in Queue:  # find all nearest Neighbors
                         #   (still coloured white) to vertex u
            for i in White:
                Distance = math.sqrt((X[i] - X[u]) * (X[i] - X[u]) +
                                     (Y[i] - Y[u]) * (Y[i] - Y[u]) +
                                     (Z[i] - Z[u]) * (Z[i] - Z[u]))
                if Distance < convert(u) + convert(i):
                    Queue.append(i)  # if you find you, put it in the que
                    White.remove(i)  # and remove it from the untouched list
        Queue.remove(u)  # remove focus from Queue
        Black.append(u)
        Fragment[frag].append(int(u + 1))  # add to group (adding 1 to start
                                           #   list at one instead of zero)

    if(len(White) != 0):  # cant move White->Queue if no more exist
        Queue.append(White[0])
        White.remove(White[0])
    frag += 1

print 'Found %s fragments' % (frag,)


for f in range(0, frag):
    A = [str(x) for x in sorted(Fragment[f])]
    print 'Fragment %s with %s atoms' % (f + 1, len(A))
    A = string.join(A, ' ')
    print A
