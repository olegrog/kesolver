#!/usr/bin/env python

import sys
import json

import math

import numpy as np
from numpy import fromfunction
from itertools import chain

from kepy.io import read_nodes, read_elements, write_elements

def sqr(x):
    return x * x

def cube(x):
    return x * x * x

def merge_dicts(*dicts):
    return dict(chain(*[d.iteritems() for d in dicts]))

seps = set(('(', ')'))

def is_const(lines):
    for line in lines:
        for w in line.split():
            if w in seps:
                continue
            try: 
                float(w)
            except ValueError:
                return False
    return True

def to_const_str(s, x, y, z):
    return "%.6g" % eval(s)

def to_const_line(line, x, y, z):
    ws = []
    for w in line.split():
        if w in seps:
            ws.append(w)
        else:
            ws.append(to_const_str(w, x, y, z))
    return ' '.join(ws)

def to_const_dict(dic, x, y, z):
    return dict((key, to_const_line(value, x, y, z)) for key, value in dic.iteritems())

def center(cell, nodes):
    sum = np.array(nodes[cell.nodes[0]])
    for i in cell.nodes[1:]:
        sum += nodes[i]
    return sum / len(cell.nodes)

if __name__ == "__main__":

    # open .kei file
    with open(sys.argv[-2], 'rb') as fd:
        data = json.load(fd)

    meshdata = data['mesh']
    nodes  = read_nodes(meshdata)
    cells  = read_elements(meshdata['cells'])
    facets = read_elements(meshdata['facets'])

    ics = data["initial_conditions"]

    new_ics = dict(ics)
    for key, value in ics.iteritems():
        if value["type"] == "maxwell":
            if not is_const([value['n'], value['T'], value['u']]):
                for cell in cells:
                    if cell.phys_index == key:
                        x, y, z = center(cell, nodes)
                        newval = to_const_dict({'n': value['n'], 'T': value['T'], 'u': value['u']}, x, y, z)
                        newkey = key + ': ' + str(newval)
                        newval['type'] = 'maxwell'
                        new_ics[newkey] = newval
                        cell.phys_index = newkey

    meshdata['cells'] = write_elements(cells)
    data['initial_conditions'] = new_ics

    bcs = data["boundary_conditions"]

    new_bcs = dict(bcs)
    for key, value in bcs.iteritems():
        if value["type"] == "diffusion":
            if not is_const([value['T'], value['u']]):
                for facet in facets:
                    if facet.phys_index == key:
                        x, y, z = center(facet, nodes)
                        newval = to_const_dict({'T': value['T'], 'u': value['u']}, x, y, z)
                        newkey = key + ': ' + str(newval)
                        newval2 = dict(value)
                        newval2['T'] = newval['T']
                        newval2['u'] = newval['u']
                        new_bcs[newkey] = newval2
                        facet.phys_index = newkey

    meshdata['facets'] = write_elements(facets)
    data['boundary_conditions'] = new_bcs
    
    with open(sys.argv[-1], 'wb') as fd: 
        json.dump(data, fd)
#        json.dump(data, fd, indent=2, sort_keys=True)

