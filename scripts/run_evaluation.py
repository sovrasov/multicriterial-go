#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import subprocess
from draw_2d import readPoints, drawPoints

matplotlib.style.use('classic')

def run_command(cmd):
    PIPE = subprocess.PIPE
    p = subprocess.Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE,
        stderr=subprocess.STDOUT)
    p.wait()
    return p.stdout.read().decode("utf-8")

def main():

    problems = ['strongin', 'fonseca', 'poloni', 'schaffer2']
    common_params = ' -r 4 -e 0.01 -m 4 -s -f tmp.csv -l 4000'
    p_values = [1, 2, 4, 8, 16]
    times = [[] for i in range(len(problems))]
    iters = [[] for i in range(len(problems))]

    for i, problem in enumerate(problems):
        for p in p_values:
            output = run_command('./multicriterial_sample ' +
            ' -p ' + str(p) + common_params + ' -n ' + problem)
            #print(output)
            trials = int(re.findall('Iterations performed: (\d+)' , output)[0])
            time = float(re.findall('Time elapsed: (\d+\.\d+)' , output)[0])
            times[i].append(time)
            iters[i].append(trials)

    for i, problem in enumerate(problems):
        print('Problem: {}'.format(problem))
        current_times = times[i]
        current_iters = iters[i]
        print('Sequeantal time: {}'.format(current_times[0]))
        print('Sequeantal iters: {}'.format(current_iters[0]))
        for j, p in enumerate(p_values):
            if p is not 1:
                print('p = {} speedup time {} iters {}'.format(p, current_times[0] / current_times[j],
                    float(current_iters[0]) / current_iters[j]))
        print('--------------------------')

if __name__ == '__main__':
    main()
