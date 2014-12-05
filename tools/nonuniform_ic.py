#!/usr/bin/env python
# attach nonuniform initial conditions to kei-file
# Usage: ./nonuniform_ic <kei-file> <macro-file> <new-kei-file>

import sys
import json
import numpy as np

import out2

def vecToStr(l):
    return '(' + ' '.join( map(str, l) ) + ')'

if __name__ == "__main__":
    
    nodes, cells = out2.readNodesCells(sys.argv[1])
    #macro = out2.readMacros(sys.argv[2], len(cells))
    macro = np.array(out2.readMacros(sys.argv[2], len(cells))).transpose()
    # open kei-file again =| may be some optimization should be there
    with open(sys.argv[1], 'r') as f:
        data = json.load(f)

    init_data = data['initial_conditions']

    values = dict()
    for phys_index in init_data:
        if init_data[phys_index]['type'] == 'maxwell':
            init_data[phys_index]['type'] = 'grad13-nonuniform'
            values[phys_index] = {}
            del init_data[phys_index]['n']
            del init_data[phys_index]['T']
            del init_data[phys_index]['u']

    for i, cell in enumerate(cells):
        if init_data[cell.phys_index]['type'] == 'grad13-nonuniform':
            values[cell.phys_index][cell.ord_index] = {
                'n': str(macro[i,0]), 'u': vecToStr(macro[i,1:4]), 'T': str(macro[i,4]),
                't': vecToStr(macro[i,5:8]), 'q': vecToStr(macro[i,8:11]), 'p': vecToStr(macro[i,11:14])
            }

    for phys_index in init_data:
        init_data[phys_index]['values'] = values[phys_index]

    # open resulting .kei file and dump the modified data 
    with open(sys.argv[3], 'w') as f:
        json.dump(data, f, indent=2, sort_keys=True)
 
