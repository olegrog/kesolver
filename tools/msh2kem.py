#!/usr/bin/env python

import sys
import string
import json
from array import array
from base64 import b64encode
from functools import reduce

from kepy.element import element, facets_of_elem

def read_string_line(f):
    return f.readline().decode();

def read_msh_node_text(fd):
    words = read_string_line(fd).split()
    transformations = [int, float, float, float]
    data = [t(w) for w, t in zip(words, transformations)]
    return tuple( data )

def read_msh_element_text(fd, ord_index=0):
    words = read_string_line(fd).split()
    data = [int(w) for w in words]

    number_of_tags = data[2]
    part_index = 0 if number_of_tags < 3 else data[6]
    elem_type = data[1]

    return element(t           =  elem_type,
                   ns          =  [x-1 for x in data[-element.number_of_nodes(elem_type):]],
                   ord_index   =  ord_index,
                   phys_index  =  data[3],
                   part_index  =  part_index)

def read_msh_physicalname(fd):
    words = read_string_line(fd).split()
    return int(words[-2]), words[-1]

def indexes_of_cells_and_facets(dimensions):
    dimension = max(dimensions)
    return ([i for i, d in enumerate(dimensions) if d == dimension],
            [i for i, d in enumerate(dimensions) if d == dimension-1])

def elem_to_key(elm, number_of_nodes):
    return reduce( lambda x, y: x*(number_of_nodes+1) + y, sorted( elm.nodes ) )

def gen_physical_map(facets, key_gen):
    key_to_physical_index = {}
    for facet in facets:
        key = key_gen(facet)
        key_to_physical_index[key] = facet.phys_index
    return key_to_physical_index

def write_nodes(nodes, data):
    data['nodes_num'] = len(nodes)
    nodes_list = [w for node in nodes for w in node]
    nodes_array = array('d', nodes_list)
    nodes_base64 = b64encode(nodes_array.tobytes()).decode()
    data['nodes'] = nodes_base64

def cell_to_dict(cell, phys_names):
    data = cell.__dict__
    phys_name = phys_names[data['phys_index']]
    # remove unnecessary quotes
    if phys_name[0] == '"' and phys_name[-1] == '"':
        phys_name = phys_name[1:-1]
    data['phys_name'] = phys_name
    del data['phys_index']
    return data

def write_elements(cells, phys_names):
    return [cell_to_dict(cell, phys_names) for cell in cells]

if __name__ == "__main__":

    read_msh_node    = read_msh_node_text
    read_msh_element = read_msh_element_text

    nodes          = []
    elems          = []
    physical_names = {}

    # open .msh file and parse it
    # the result is placed on the variables above
    with open(sys.argv[-2], 'rb') as mshfile:
        line = " "
        while line:
            line = read_string_line(mshfile)
            if line == "$Nodes\n":
                n = int( read_string_line(mshfile) )
                for i in range(n):
                    i, x1, x2, x3 = read_msh_node(mshfile)
                    nodes.append( (x1, x2, x3) )
                read_string_line(mshfile) # EndNodes

            elif line == "$Elements\n":
                n = int( read_string_line(mshfile) )
                elems = [ read_msh_element(mshfile) for i in range(n) ]
                read_string_line(mshfile) # EndElements

            elif line == "$PhysicalNames\n":
                n = int( read_string_line(mshfile) )
                for i in range(n):
                    index, name = read_msh_physicalname(mshfile)
                    physical_names[index] = name
                read_string_line(mshfile) # EndPhysicalNames

    # we should do this because in .msh file part_indexes start from 1
    number_of_partitions = max([e.part_index for e in elems])
    if number_of_partitions > 0:
        for e in elems:
            e.part_index -= 1

    # having been known the dimension of the mesh
    # we are getting the indexes of cells and facets
    indexes_of_cells, indexes_of_boundary_facets = \
        indexes_of_cells_and_facets( [e.dimension() for e in elems] )

    # auxiliary function, helps to get elements from a list by their indexes
    at = lambda array, indexes : [array[i] for i in indexes]

    # here we are spliting elements to cells and facets
    boundary_facets  =  at(elems, indexes_of_boundary_facets)
    cells            =  at(elems, indexes_of_cells)

    # generate a map of which stores physical_indexes of boundary facets by their keys
    bfk_to_physical_index = gen_physical_map(
        boundary_facets,
        lambda elm : elem_to_key( elm, len(nodes) )
    )

    # all facets of mesh, not just boundary facets
    facets = []
    # a map which stores indexes of facets on the list above by their keys
    table = {}

    # generate all facets of the mesh, result is placed in variable facets
    # these two fors mean for all facets in cells
    for i, cell in enumerate(cells):
        for facet in facets_of_elem(cell):
            key = elem_to_key(facet, len(nodes))
            if not key in table:    # the facet is not on the list yet
                                    # put it
                table[key] = len(facets)
                facet.neigbors = [i]
                facet.phys_index = (key in bfk_to_physical_index) and \
                                   bfk_to_physical_index[key]     or  \
                                   cell.phys_index
                facet.part_index = [cell.part_index]
                facets.append(facet)

            else:        # the facet is already in the list
                index_of_facet = table[key] 
                facet = facets[index_of_facet]
                if len(facet.neigbors) != 1:
                    raise
                neigbor_cell = facet.neigbors[0]
                cell.neigbors.append(neigbor_cell)
                cells[neigbor_cell].neigbors.append(i)
                facet.neigbors.append(i)
                if  facet.part_index[0] != cell.part_index:
                    facet.part_index.append(cell.part_index)

    # remove facets with only one neigbor having no boundary conditions on them
    facets = [ f for f in facets if ( len(f.neigbors) == 2 )
                                or elem_to_key(f, len(nodes)) in bfk_to_physical_index ]

    # TODO split boundary facets neigboring two cells

    # number cells and facets
    for i, c in enumerate(cells):
        c.ord_index = i
    for i, f in enumerate(facets):
        f.ord_index = i

    # dump data to in json file
    meshdata = {}

    write_nodes(nodes, meshdata)
    meshdata['cells']  = write_elements(cells,  physical_names)
    meshdata['facets'] = write_elements(facets, physical_names)
    meshdata['type'] = 'unstructured'

    with open(sys.argv[-1], 'w') as fd:
        json.dump(meshdata, fd, ensure_ascii=False)
#        json.dump(meshdata, fd, indent=2, sort_keys=True)

