import json
from array import array

import numpy
import math

from base64 import b64decode
from kepy.element import element

def read_nodes(data):
    nodes_array = array('d', b64decode(data['nodes']))
    nodes = [numpy.array(x) for x in zip(nodes_array[0::3],
                                         nodes_array[1::3],
                                         nodes_array[2::3])]
    assert data['nodes_num'] == len(nodes)
    return nodes

def read_elements(data):
    return [element(t           =  c['type'],
                    ns          =  c['nodes'],
                    ord_index   =  c['ord_index'],
                    phys_index  =  c['phys_name'],
                    part_index  =  c['part_index'],
                    neigbors    =  c['neigbors']) for c in data]

def readNodesElems(filename):
    with open(filename, "rb") as fd:
        data  = json.load(fd)
        meshdata = data['mesh']
        nodes = read_nodes(meshdata)
        cells  = read_elements(meshdata['cells'])
        facets = read_elements(meshdata['facets'])

    return nodes, cells, facets

def write_elements(cells):
    return [{'type':       cell.type,
             'nodes':      cell.nodes,
             'ord_index':  cell.ord_index,
             'phys_name':  cell.phys_index,
             'part_index': cell.part_index,
             'neigbors':   cell.neigbors} for cell in cells]

def calc_timestep(data, cells, nodes):

    def sqr(x):
        return numpy.dot(x, x)

    def lMin(cell, nodes):
        l_min = min( [ sqr(nodes[cell.nodes[i]]-nodes[cell.nodes[j]]) 
                       for i in range(len(cell.nodes)) 
                           for j in range(i+1, len(cell.nodes)) ] )
        return math.sqrt(l_min)

    l_min    = min( [ lMin(cell, nodes) for cell in cells ] )
    curnt    = data['curnt_limit']
    cut      = data["gas"]["cut"]
    timestep = curnt * l_min / cut * 2

    return timestep

