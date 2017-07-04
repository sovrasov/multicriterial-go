#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.style.use('classic')

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

def drawPoints(Y, W, filename, ext):
    if len(W[0]) == 2:
        plt.xlabel(r'$f_1$')
        plt.ylabel(r'$f_2$')

        plt.plot(W[:,0], W[:,1], 'ro')

        plt.grid()
        plt.savefig(filename + '.' + ext, format = ext, dpi = 200)
        plt.clf()

    if len(Y[0]) == 2:
        plt.xlabel(r'$y_1$')
        plt.ylabel(r'$y_2$')

        plt.plot(Y[:,0], Y[:,1], 'go')

        plt.grid()
        plt.savefig(filename + '_points.' + ext, format = ext, dpi = 200)
        plt.clf()

def main():

    Y, W = readPoints(sys.argv[1])
    drawPoints(Y, W, sys.argv[2], 'png')

if __name__ == '__main__':
    main()
