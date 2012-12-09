#!/usr/bin/env python

import sys, struct, string
import numpy

def typeOfElem(elm):
    return elm[0]

def dimensionOfElem(elm):
    return dimensions[ typeOfElem(elm) ]


# [elm_type, nodes]
facets_of_tetrahedron = [[2, 1, 2, 3], [2, 1, 3, 4], [2, 1, 4, 2], [2, 2, 4, 3]]
facets_of_prism = [[2, 1, 2, 3], [3, 1, 3, 6, 4], [3, 3, 2, 5, 6], [3, 1, 4, 5, 2], [2, 4, 6, 5]]
facets_of_hexahedron = [[3, 1, 2, 3, 4], [3, 1, 5, 6, 2], [3, 2, 6, 7, 3], [3, 3, 7, 8, 4], [3, 1, 4, 8, 5], [3, 5, 8, 7, 6]]

def facets_of_elem(elm):
    dic_of_facets = {4 : facets_of_tetrahedron,
                     5 : facets_of_hexahedron,
                     6 : facets_of_prism}
    facets = dic_of_facets[elm.type]
    for facet in facets:
        yield element(facet[0], [elm.nodes[i-1] for i in facet[1:]])

def read_node_text(fd):
    words = fd.readline().split()
    transformations = [int, float, float, float]    
    data = [t(w) for w, t in zip(words, transformations)]
    return tuple( data )
    
def read_node_binary(fd):
    data = fd.read(4+8+8+8)
    return struct.unpack('=iddd', data)

class element:
    def __init__(self, t, ns, phys_index=0, part_index=0):
        self.type = t
        self.nodes = ns
        self.phys_index = phys_index
        self.part_index = part_index


    @staticmethod
    def number_of_nodes(t):
        number_of_nodes_list = [-1, 2, 3, 4, 4, 8, 6, 5, 3, 6, 9, 10, 27, 18, 14, 1]
        return number_of_nodes_list[t]
        

    def dimension(self):
        dimensions_list = [-1, 1, 2, 2, 3, 3, 3, 3, 1, 2, 2, 3, 3, 3, 3, 0] 
        return dimensions_list[self.type]

def read_element_text(fd):
    words = fd.readline().split()
    data = [int(w) for w in words]

    number_of_tags = data[2]
    part_index = 0 if number_of_tags < 3 else data[6]
    elem_type = data[1]

    return element(elem_type, 
                   [x-1 for x in data[-element.number_of_nodes(elem_type):]],
                   data[3], # phys_index
                   part_index)

def read_physicalname(fd):
    words = fd.readline().split()
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

def write_node_text(fd, node):
    s = "%.10f %.10f %.10f\n" % (node[0], node[1], node[2])
    fd.writelines([s])

def list_to_str(elm):
    return str(len(elm)) + ' ' + string.join( map(str, elm) )

def write_cell_text(fd, cell):
    s = str(cell.type) + ' ' + string.join( map(str, cell.nodes) ) + ' ' + \
        list_to_str(cell.neighbors) + ' ' + \
        str(cell.part_index) + ' ' + str(cell.phys_index) + '\n'
    fd.writelines([s])

def write_facet_text(fd, facet):
#   neighbors.reverse()
    s = str(facet.type) + ' ' + string.join( map(str, facet.nodes) ) + ' ' + \
        list_to_str(facet.neighbors) + ' ' + \
        list_to_str(facet.part_index) + ' ' + str(facet.phys_index) + '\n'
    fd.writelines([s])

if __name__ == "__main__":

    read_node = read_node_text
    read_element = read_element_text
    if sys.argv[1] == '-b':
        read_node = read_node_binary
        read_element = read_element_binary

    nodes = [] 
    physical_names = {}
    elems = []

    with open(sys.argv[-2], 'rb') as mshfile:
        line = " "
        while line:
            line = mshfile.readline()
            if line == "$Nodes\n":
                n = int( mshfile.readline() )
                for i in range(n):
                    i, x1, x2, x3 = read_node(mshfile)
                    nodes.append( numpy.array( [x1, x2, x3] ) )
                mshfile.readline() # EndNodes

            elif line == "$Elements\n":
                n = int( mshfile.readline() )
                elems = [ read_element(mshfile) for i in range(n) ]
                mshfile.readline() # EndElements

            elif line == "$PhysicalNames\n":
                n = int( mshfile.readline() )
                for i in range(n):
                    index, name = read_physicalname(mshfile)
                    physical_names[index] = name
                mshfile.readline() # EndPhysicalNames


    number_of_partitions = max([e.part_index for e in elems])
    if number_of_partitions > 0:
        for e in elems:
            e.part_index -= 1
            print e.part_index

    indexes_of_cells, indexes_of_boundary_facets = \
        indexes_of_cells_and_facets( [e.dimension() for e in elems] )

    at = lambda array, indexes : [array[i] for i in indexes] 

    boundary_facets = at(elems, indexes_of_boundary_facets)
    cells = at(elems, indexes_of_cells)

    key_to_physical_index = gen_physical_map( 
        boundary_facets,
        lambda elm : elem_to_key( elm, len(nodes) )
        )

    table = {}
    facets = []

    for i, cell in enumerate(cells):
        cell.neighbors = []
        for facet in facets_of_elem(cell):
            key = elem_to_key(facet, len(nodes))
            if key in table:
                index_of_facet = table[key] 
                facet = facets[index_of_facet]
                if len(facet.neighbors) != 1:
                    raise
                neighbor_cell = facet.neighbors[0]
                cell.neighbors.append(neighbor_cell)
                cells[neighbor_cell].neighbors.append(i)
                facet.neighbors.append(i)
                if  facet.part_index[0] != cell.part_index:
                    facet.part_index.append(cell.part_index)

            else:
                table[key] = len(facets)
                facet.neighbors = [i]
                facet.phys_index = (key in key_to_physical_index) and \
                                   key_to_physical_index[key]     or  \
                                   cell.phys_index
                facet.part_index = [cell.part_index]
                facets.append(facet)

    write_node = write_node_text
    write_cell = write_cell_text
    write_facet = write_facet_text

    with open(sys.argv[-1], 'wb') as kshfile:

        kshfile.writelines(["<nodes len=\"%d\">\n" % len(nodes)])
        for node in nodes:
            write_node(kshfile, node)
        kshfile.writelines(["</nodes>\n"])

        kshfile.writelines(["<cells len=\"%d\">\n" % len(cells)])
        for cell in cells:
            write_cell(kshfile, cell)
        kshfile.writelines(["</cells>\n"])

        kshfile.writelines(["<facets len=\"%d\">\n" % len(facets)])
        for facet in facets:
            write_facet(kshfile, facet)
        kshfile.writelines(["</facets>\n"])

        kshfile.writelines(["<physical len=\"%d\">\n" % len(physical_names)])
        for index, value in physical_names.iteritems():
            kshfile.writelines(["%d %s\n" % (index, value)])
        kshfile.writelines(["</physical>\n"])
        
