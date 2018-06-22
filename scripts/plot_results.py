#!/usr/bin/env python3

import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D as plt3d
from scipy.optimize import curve_fit
from cycler import cycler
from numpy import log2
import os

# https://en.wikipedia.org/wiki/List_of_colors:_A%E2%80%93F
COLORS = [
    '#007FFF',  # (b) azure
    '#FD3F92',  # (r) french fuchsia
    '#84DE02',  # (g) alien ampit
    '#9966CC',  # (m) ametist
    '#FFA812',  # (y) dark tangerine
    '#B87333',  # (b) copper
    '#3B444B',  # (k) arsenic
    '#5F9EA0',  # cadet blue
]

# https://en.wikipedia.org/wiki/List_of_colors:_N%E2%80%93Z
MCOLORS = [
    '#B0E0E6',  # powder blue
    '#FDD5B1',  # light appricot
    '#F3E5AB',  # vanilla
    '#BEBEBE',  # gray
]

# file format:
# name      size      time


def get_results(filename):
    results = {}
    with open(filename) as f:
        for line in f:
            name, size, time = line.split()
            name, size, time = str(name), int(size), float(time)
            if name not in results:
                results[name] = {
                    'size': [],
                    'time': []
                }
            results[name]['size'] += [size]
            results[name]['time'] += [time]
    return results


def set_plot(title, x, y, xscale, yscale):
    plt.figure(figsize=(12, 12))
    plt.title(title, size=30, family='monospace', weight='bold', backgroundcolor=MCOLORS[0], color='k')
    plt.xscale(xscale)
    plt.yscale(yscale)
    plt.xlabel(x, size=24, family='monospace', weight='bold', backgroundcolor=MCOLORS[1], color='k')
    plt.rc('axes', prop_cycle=(cycler('color', COLORS)))
    plt.ylabel(y, size=24, family='monospace', weight='bold', backgroundcolor=MCOLORS[1], color='k')


def add_to_plot(results, name, x, y):
    return plt.plot(results[name][x], results[name][y], '-o', label=name)


def plot_benchmarks(title, results, filename, x, y, xscale='linear', yscale='log'):
    set_plot(title, x, y, xscale, yscale)
    plots = []
    for name in results:
        plots += add_to_plot(results, name, x, y)
        print('poly_fit[' + name + ']: ',
              ["%0.1f" % float(x * 1e5) for x in
               np.polyfit(results[name]['size'],
                         results[name]['time'],
                         4)])
    plt.legend(handles=plots, prop={'size': 27, 'family': 'monospace'})
    if os.path.isfile(filename):
        os.remove(filename)
    plt.savefig(filename)
    plt.clf()
    plt.cla()


if __name__ == "__main__":
    plt.switch_backend('agg')
    plot_benchmarks('test', get_results('test.txt'), 'test.png',
                    'size', 'time',
                    'log', 'log')
