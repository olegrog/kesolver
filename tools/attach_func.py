#!/usr/bin/env python

import sys
import json

import struct
import array

import numpy as np

from base64 import b64encode, b64decode

from kepy.io import readNodesElems, write_elements

sqr = lambda x: x*x

def find_center(nodes, nodes_indexes):
    sum = np.zeros_like(nodes[nodes_indexes[0]])
    for i in nodes_indexes:
        sum += nodes[i]
    sum /= len(nodes_indexes)
    return sum

if __name__ == "__main__":
    
    # open .kei file of the calculation from where we will take a distribution function
    nodes, cells, _ = readNodesElems(sys.argv[-4])
    cells_num = len(cells)
    print cells_num

    centers = [ find_center(nodes, cell.nodes) for cell in cells ] 
    for c in centers:
        print c

    func = range(cells_num)

    # open a file with a distribution function
    with open(sys.argv[-3], 'rb') as fd:
        while True:
            ii = fd.read(8)
            if not ii: break

            i, s = struct.unpack('=ii', ii)
#            print i, s

            a = array.array('d')
            a.fromfile(fd, s)

            b = b64encode(a)

            func[i] = b

    # open .kei file of the calculation where we will put the distribution function
    keifile2 = sys.argv[-2]
    nodes2, cells2, _ = readNodesElems(keifile2)
    # open it again =| may be some optimization should be there
    with open(sys.argv[-2], 'r') as fd:
        data = json.load(fd)

    init_data = data['initial_conditions']

    b = [False for i in range(100)]

    for cell in cells2:
        if init_data[cell.phys_index]['type'] == 'from_func':

            center = find_center(nodes2, cell.nodes)

            for i, c in enumerate(centers):
                if np.abs(c - center)[0] < 1e-6:
                    name = '_from_func_%d_' % i
                    if not name in init_data: 
                        init_data[name] = {'type': 'raw', 'raw': func[i]}
                    cell.phys_index = name
                    break

    data['mesh']['cells'] = write_elements(cells2)

    # open resulting .kei file and dump the modified data 
    with open(sys.argv[-1], 'w') as fd:
        json.dump(data, fd, indent=2, sort_keys=True)
 
