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


    if symmetry == "Cylindrical":

        v = float(data["gas"].get("v", "0."))

        x = np.fromfunction(lambda i, j: v + (i + 0.5 - rad) * cut / rad,
                            (2 * rad, rad), dtype=float)
        y = np.fromfunction(lambda i, j: (j + 0.5) * cut / rad,
                            (2 * rad, rad), dtype=float)
        ax = [x, y]

        r = sqr(x-v) + sqr(y)
        
        vol = y
        d3v = 2 * math.pi * sqr(cut / rad)

    elif symmetry == "Cartesian":

        v = map(float, data["gas"].get("v", "( 0. 0. 0. )")[1:-1].split() )
        vx, vy, vz = v

        x = np.fromfunction(lambda i, j, k: vx + (i + 0.5 - rad) * cut / rad,
                            (2 * rad, 2 * rad, 2 * rad), dtype=float)
        y = np.fromfunction(lambda i, j, k: vy + (j + 0.5 - rad) * cut / rad,
                            (2 * rad, 2 * rad, 2 * rad), dtype=float)
        z = np.fromfunction(lambda i, j, k: vz + (k + 0.5 - rad) * cut / rad,
                            (2 * rad, 2 * rad, 2 * rad), dtype=float)
        ax = [x, y, z]

        r = sqr(x-vx) + sqr(y-vy) + sqr(z-vz)

        vol = np.ones_like(y)
        d3v = cube(cut / rad)

    def make_e(symm, ax, u=None):
        if symm == "Cylindrical":
            x, y = ax
            if u is None:
                u = 0.
            return sqr(x-u) + sqr(y)
        elif symm == "Cartesian":
            x, y, z = ax
            if u is None:
                ux, uy, uz = 0., 0., 0.
            else:
                ux, uy, uz = u
            return sqr(x-ux) + sqr(y-uy) + sqr(z-uz)

    def str_to_u(symm, u):
        if symm == "Cylindrical":
            return float(u)
        elif symm == "Cartesian":
            return map(float, u[1:-1].split())

    init_data = data["initial_conditions"]

    new = {}
    for key, value in init_data.iteritems():

        if value["type"] == "maxwell":

            n = float(value["n"])
            T = float(value["T"])
            u = str_to_u(symmetry, value("u"))
            e = make_e(symmetry, ax, u)

            f = np.exp( - 0.5 * e / T )
            f *= vol
            lin_f = f[r < sqr(cut)]
            s = np.sum(lin_f) * d3v
            lin_f = n * lin_f / s

        if value["type"] == "bkw":

            n = 1.
            xi = float(value["xi"])
            tau = 1. - xi

            e = make_e(symmetry, ax)
            e1 = 0.5 * e / tau

            f = 1. / math.sqrt(cube(2 * math.pi * tau)) * np.exp(-e1) * ( 1 + (1 - tau) / tau * (e1 - 1.5) )
            f *= vol

            lin_f = f[r < sqr(cut)]

        if value["type"] == "sum":

            lin_f = np.zeros_like(x)[sqr(x) + sqr(y) < sqr(cut)]

            for macrodata in value["maxwellians"]:

                n = float(macrodata["n"])
                T = float(macrodata["T"])
                u = str_to_u(symmetry, macrodata["u"])
                e = make_e(symmetry, ax, u)

                f = y * np.exp( - 0.5 * e / T )
                f_data = f[r < sqr(cut)]
                s = np.sum(lin_f) * 2 * math.pi * sqr(cut / rad)
                f_data = n * lin_f / s
                
                lin_f += f_data

        doubles = array.array('d', lin_f)
        b64data = b64encode(doubles)

        new[key] = {'type': 'raw',
                    'raw': b64data}

    data['initial_conditions'] = new
    
    with open(sys.argv[-1], 'wb') as fd: 
#        json.dump(data, fd)
        json.dump(data, fd, indent=2, sort_keys=True)

