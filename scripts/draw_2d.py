#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import numpy as np
import matplotlib.pyplot as plt

def readPoints(fileName):
    lines = [line.rstrip('\n') for line in open(fileName)]

    n = int(lines.pop(0))
    m = int(lines.pop(0))

    pointsY = []
    pointsW = []

    for line in lines:
        terms = line.split(';')[0].split(',')

        W, Y = [], []

        for i in range(n):
            Y.append(float(terms[i]))
        for i in range(m):
            W.append(float(terms[n + i]))

        pointsY.append(np.array(Y))
        pointsW.append(np.array(W))

    return np.array(pointsY), np.array(pointsW)

def main():

    Y, W = readPoints(sys.argv[1])

    plt.xlabel(r'$f_1$')
    plt.ylabel(r'$f_2$')

    plt.plot(W[:,0], W[:,1], 'ro')

    plt.grid()
    plt.savefig('plot2d.png', format = 'png', dpi = 200)
    plt.clf()

if __name__ == '__main__':
    main()
