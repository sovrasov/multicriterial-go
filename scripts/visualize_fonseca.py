#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import math
from draw_2d import readPoints

matplotlib.style.use('classic')

def poloni_1(y):
    shift = 1. / pow(len(y), 0.5)
    eArg = 0
    for i in range(len(y)):
        eArg -= pow(y[i] - shift, 2.)
    return 1 - math.exp(eArg)

def poloni_2(y):
    shift = 1. / pow(len(y), 0.5)
    eArg = 0
    for i in range(len(y)):
        eArg -= pow(y[i] + shift, 2.)
    return 1 - math.exp(eArg)

def main():

    Y, W = readPoints(sys.argv[1])

    n = len(Y[0])
    X_p = np.linspace(-1. / pow(n, 0.5), 1. / pow(n, 0.5), 20)

    f1_p = [poloni_1([x, x]) for x in X_p]
    f2_p = [poloni_2([x, x]) for x in X_p]

    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlabel(r'$f_1$', fontsize=20)
    plt.ylabel(r'$f_2$', fontsize=20)

    plt.plot(W[:,0], W[:,1], 'ro', label='Numerical solution')
    plt.plot(f1_p, f2_p, 'g-', label='Slater set')
    plt.legend(loc = 'best', fontsize = 20)

    plt.grid()
    plt.savefig(sys.argv[1].split('.')[0] + '.pdf', format = 'pdf', dpi = 200)
    plt.clf()


if __name__ == '__main__':
    main()
