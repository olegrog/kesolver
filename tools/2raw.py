#!/usr/bin/env python

import sys
import json

import math

import numpy as np
from numpy import fromfunction

import array
from base64 import b64encode

from kepy.ximesh import read_ximesh, make_e

def sqr(x):
    return x * x

def cube(x):
    return x * x * x

if __name__ == "__main__":

    # open .kei file
    with open(sys.argv[-2], 'rb') as fd:
        data = json.load(fd)

    symmetry, rad, circl, ax, vol, r, d3v = read_ximesh(data)

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
            lin_f = f[circl]
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

            lin_f = f[circl]

        if value["type"] == "sum":

            lin_f = np.zeros_like(x)[circl]

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

