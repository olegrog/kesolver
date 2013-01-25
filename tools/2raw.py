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
    v   = data["gas"]["v"]
    symmetry = data["gas"]["symmetry"]

#    print cut, rad, v, symmetry

    init_data = data["initial_conditions"]

    new = {}
    for key, value in init_data.iteritems():

        if value["type"] == "maxwell":

            n = float(value["n"])
            T = float(value["T"])

            if symmetry == "Cylindrical":

                u = float(value["u"])
                v = float(v)

                x = np.fromfunction(lambda i, j: v + (i + 0.5 - rad) * cut / rad,
                                    (2 * rad, rad), dtype=float)
                y = np.fromfunction(lambda i, j: (j + 0.5) * cut / rad,
                                    (2 * rad, rad), dtype=float)

                f = y * np.exp( - 0.5 * ( sqr(x-u) + sqr(y) ) / T )
                lin_f = f[sqr(x-v) + sqr(y) < sqr(cut)]
                s = np.sum(lin_f) * 2 * math.pi * sqr(cut / rad)
                lin_f = n * lin_f / s

                doubles = array.array('d', lin_f)
                b64data = b64encode(doubles)

            elif symmetry == "Cartesian":

                ux, uy, uz = map(float, value["u"][1:-1].split())
                vx, vy, vz = map(float, v[1:-1].split())

                x = np.fromfunction(lambda i, j, k: vx + (i + 0.5 - rad) * cut / rad,
                                    (2 * rad, 2 * rad, 2 * rad), dtype=float)
                y = np.fromfunction(lambda i, j, k: vy + (j + 0.5 - rad) * cut / rad,
                                    (2 * rad, 2 * rad, 2 * rad), dtype=float)
                z = np.fromfunction(lambda i, j, k: vz + (k + 0.5 - rad) * cut / rad,
                                    (2 * rad, 2 * rad, 2 * rad), dtype=float)

                f = np.exp( - 0.5 * ( sqr(x-ux) + sqr(y-uy) + sqr(z-uz) ) / T )
                lin_f = f[sqr(x-vx) + sqr(y-vy) + sqr(z-vz) < sqr(cut)]
                s = np.sum(lin_f) * cube(cut / rad)
                lin_f = n * lin_f / s

                doubles = array.array('d', lin_f)
                b64data = b64encode(doubles)

            new[key] = {'type': 'raw',
                        'raw': b64data}

    data['initial_conditions'] = new
    
    with open(sys.argv[-1], 'wb') as fd: 
#        json.dump(data, fd)
        json.dump(data, fd, indent=2, sort_keys=True)

