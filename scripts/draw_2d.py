#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
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
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.xlabel(r'$f_1$', fontsize=20)
        plt.ylabel(r'$f_2$', fontsize=20)

        plt.plot(W[:,0], W[:,1], 'ro')

        plt.grid()
        plt.savefig(filename + '.' + ext, format = ext, dpi = 200)
        plt.clf()
    elif len(W[0]) == 3:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(W[:,0], W[:,1], W[:,2], c='r', marker='o')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        plt.show()

    if len(Y[0]) == 2:
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.xlabel(r'$y_1$', fontsize=20)
        plt.ylabel(r'$y_2$', fontsize=20)

        plt.plot(Y[:,0], Y[:,1], 'go')

        plt.grid()
        plt.savefig(filename + '_points.' + ext, format = ext, dpi = 200)
        plt.clf()

def main():

    Y, W = readPoints(sys.argv[1])
    drawPoints(Y, W, sys.argv[2], 'png')

if __name__ == '__main__':
    main()
