#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import re
import subprocess

if not (len(sys.argv) is 2):
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib
    from draw_2d import readPoints, drawPoints
    matplotlib.style.use('classic')

def run_command(cmd):
    PIPE = subprocess.PIPE
    p = subprocess.Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE,
        stderr=subprocess.STDOUT)
    p.wait()
    return p.stdout.read().decode("utf-8")

def main():

    tmp_filename = 'tmp.csv'
    problems = ['strongin', 'fonseca', 'poloni', 'schaffer2']
    #problems = ['fonseca']
    problems_params_dict = {problem: '' for problem in problems}
    problems_params_dict['schaffer2'] = ' -e 0.001 '
    problems_params_dict['poloni'] = ' -l 4000 '
    problems_params_dict['fonseca'] = ' -l 10000 -d 3'
    problems_params_dict['strongin'] = ' '

    common_params = ' -r 4.5 -m 4 -s -f ' + tmp_filename + ' '

    p_values = [1, 2, 4, 8, 16]
    times = [[] for i in range(len(problems))]
    iters = [[] for i in range(len(problems))]
    s_sizes = [[] for i in range(len(problems))]

    for i, problem in enumerate(problems):
        for p in p_values:
            output = run_command('./multicriterial_sample ' +
            ' -p ' + str(p) + common_params + ' -n ' + problem + problems_params_dict[problem])

            trials = int(re.findall('Iterations performed: (\d+)', output)[0])
            time = float(re.findall('Time elapsed: (\d+\.\d+)', output)[0])
            times[i].append(time)
            iters[i].append(trials)
            s_sizes[i].append(int(re.findall('Number of weak optimal points: (\d+)', output)[0]))

            if p is 1 and not (len(sys.argv) is 2):
                Y, W = readPoints(tmp_filename)
                drawPoints(Y, W, problem, 'pdf')
        print('Problem: {}'.format(problem))
        current_times = times[i]
        current_iters = iters[i]
        current_sizes = s_sizes[i]
        print('Sequeantal time: {}'.format(current_times[0]))
        print('Sequeantal iters: {}'.format(current_iters[0]))
        iters_line = '{}({}) & '.format(current_iters[0], current_sizes[0])
        speedup_line = ''
        speedup_time_line = ''
        for j, p in enumerate(p_values):
            if p is not 1:
                speedup_line += '{0:.2f} & '.format(float(current_iters[0]) / current_iters[j])
                iters_line += '{}({}) & '.format(current_iters[j], current_sizes[j])
                speedup_time_line += '{0:.2f} & '.format(float(current_times[0]) / current_times[j])
                print('p = {} speedup time {} iters {}'.format(p, current_times[0] / current_times[j],
                float(current_iters[0]) / current_iters[j]))
        print(iters_line)
        print(speedup_line)
        print(speedup_time_line)
        print('--------------------------')

if __name__ == '__main__':
    main()
