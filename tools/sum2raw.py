#!/usr/bin/env python

import sys
import json

import math

import numpy as np
from numpy import fromfunction

import array
from base64 import b64encode

def sqr(x):
    return x * x

def cube(x):
    return x * x * x

if __name__ == "__main__":

    # open .kei file
    with open(sys.argv[-2], 'rb') as fd:
        data = json.load(fd)

    cut = float(data["gas"]["cut"])
    rad = int(data["gas"]["rad"])
    symmetry = data["gas"]["symmetry"]

#    print cut, rad, v, symmetry

    init_data = data["initial_conditions"]

    new = {}
    for key, value in init_data.iteritems():

        if value["type"] == "sum":

            x = np.fromfunction(lambda i, j: (i + 0.5 - rad) * cut / rad,
                                (2 * rad, rad), dtype=float)
            y = np.fromfunction(lambda i, j: (j + 0.5) * cut / rad,
                                (2 * rad, rad), dtype=float)

            f_data = np.zeros_like(x)[sqr(x) + sqr(y) < sqr(cut)]

            for macrodata in value["maxwellians"]:

                n = float(macrodata["n"])
                T = float(macrodata["T"])
                u = float(macrodata["u"])

                f = y * np.exp( - 0.5 * ( sqr(x-u) + sqr(y) ) / T )
                lin_f = f[sqr(x) + sqr(y) < sqr(cut)]
                s = np.sum(lin_f) * 2 * math.pi * sqr(cut / rad)
                lin_f = n * lin_f / s
                
                f_data += lin_f

            doubles = array.array('d', f_data)
            b64data = b64encode(doubles)

            new[key] = {'type': 'raw',
                        'raw': b64data}

    data['initial_conditions'] = new
    
    with open(sys.argv[-1], 'wb') as fd: 
#        json.dump(data, fd)
        json.dump(data, fd, indent=2, sort_keys=True)

