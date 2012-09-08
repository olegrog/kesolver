#!/usr/bin/env python

import sys, struct, string
import numpy

def numberOfNodes(elm_type):
	number_of_nodes = [-1, 2, 3, 4, 4, 8, 6, 5, 3, 6, 9, 10, 27, 18, 14, 1]
	return number_of_nodes[elm_type]

def typeOfElem(elm):
	return elm[0]

def dimensionOfElem(elm):
	dimensions = [-1, 1, 2, 2, 3, 3, 3, 3, 1, 2, 2, 3, 3, 3, 3, 0] 
	return dimensions[ typeOfElem(elm) ]


# [elm_type, nodes]
facets_of_tetrahedron = [[2, 1, 2, 3], [2, 1, 3, 4], [2, 1, 4, 2], [2, 2, 4, 3]]
facets_of_prism = [[2, 1, 2, 3], [3, 1, 3, 6, 4], [3, 3, 2, 5, 6], [3, 1, 4, 5, 2], [2, 4, 6, 5]]
facets_of_hexahedron = [[3, 1, 2, 3, 4], [3, 1, 5, 6, 2], [3, 2, 6, 7, 3], [3, 3, 7, 8, 4], [3, 1, 4, 8, 5], [3, 5, 8, 7, 6]]

def facetsOfElem(elm):
	dic_of_facets = {4 : facets_of_tetrahedron,	5 : facets_of_hexahedron, 6 : facets_of_prism}
	facets = dic_of_facets[typeOfElem(elm)]
	for facet in facets:
		yield [facet[0]] + [elm[i] for i in facet[1:]]


def readNode_text(fd):
	words = fd.readline().split()
	transformations = [int, float, float, float]	
	data = [t(w) for w, t in zip(words, transformations)]
	return tuple( data )
	
def readNode_binary(fd):
	data = fd.read(4+8+8+8)
	return struct.unpack('=iddd', data)

def readElement_text(fd):
	words = fd.readline().split()
	data = [int(w) for w in words]
	return data

def readElement_binary(fd):
	data1 = fd.read(4+4+4)
	elm_type, num_elm_follow, number_of_tags = struct.unpack('iii', data1)
	size = 1 + number_of_tags + numberOfNodes(elm_type)
	data2 = fd.read( 4*size )
	data = struct.unpack(str(size) + 'i', data2)
	return [ data[0], elm_type, number_of_tags ] + list( data[1:(1+number_of_tags)] ) + \
		list(data[-numberOfNodes(elm_type):])

def readPhysicalName(fd):
	words = fd.readline().split()
	return int(words[-2]), words[-1]

def indexesOfCellsAndFacets(dimensions):
	dimension = max(dimensions)
	return ([i for i, d in enumerate(dimensions) if d == dimension],
		[i for i, d in enumerate(dimensions) if d == dimension-1])

def keyOfElem(elm, number_of_nodes):
	return reduce( lambda x, y: x*(number_of_nodes+1) + y, sorted( elm[1:] ) )

def keyToPhysicalIndex(facets, physical_names, keyGenerator):
	key_to_physical_index = {}
	for facet, physical_name in zip(facets, physical_names):
		key = keyGenerator(facet)
		key_to_physical_index[key] = physical_name
	return key_to_physical_index

def writeNode_text(fd, node):
	s = "%.10f %.10f %.10f\n" % (node[0], node[1], node[2])
	fd.writelines([s])

def elemToStr(elm):
	return string.join( map(str, elm) )

def listToStr(elm):
	return str(len(elm)) + ' ' + string.join( map(str, elm) )

def writeCell_text(fd, cell, neighbors, partition_index, physical_index):
	s = elemToStr(cell) + ' ' + listToStr(neighbors) + ' ' + \
		str(partition_index) + ' ' + str(physical_index) + '\n'
	fd.writelines([s])

def writeFacet_text(fd, facet, neighbors, partition, physical_index):
#	neighbors.reverse()
	s = elemToStr(facet) + ' ' + listToStr(neighbors) + ' ' + \
		listToStr(partition) + ' ' + str(physical_index) + '\n'
	fd.writelines([s])

readNode = readNode_text
readElement = readElement_text
if sys.argv[1] == '-b':
	readNode = readNode_binary
	readElement = readElement_binary

nodes = [] 

physical_names = {}

elems = []
partion_indexes_of_elems = [] 
physical_indexes_of_elems = [] 

with open(sys.argv[-2], 'rb') as mshfile:
	line = " "
	while line:
		line = mshfile.readline()
		if line == "$Nodes\n":
			n = int( mshfile.readline() )
			for i in range(n):
				i, x1, x2, x3 = readNode(mshfile)
				nodes.append( numpy.array( [x1, x2, x3] ) )
			mshfile.readline() # EndNodes

		elif line == "$Elements\n":
			n = int( mshfile.readline() )
			for i in range(n):
				data = readElement(mshfile)
				elm_type = data[1]
				elem = [elm_type] + [i-1 for i in data[-numberOfNodes(elm_type):]]
				elems.append(elem)
				physical_indexes_of_elems.append(data[3])
				number_of_tags = data[2]
				if number_of_tags < 3:
					partion_index = 0
				else:
					partion_index = data[6]
#				partion_index = 0
				partion_indexes_of_elems.append(partion_index)
			mshfile.readline() # EndElements

		elif line == "$PhysicalNames\n":
			n = int( mshfile.readline() )
			for i in range(n):
				index, name = readPhysicalName(mshfile)
				physical_names[index] = name
			mshfile.readline() # EndPhysicalNames

print max(partion_indexes_of_elems)
if max(partion_indexes_of_elems) > 0:
	for i in range(len(partion_indexes_of_elems)):
		partion_indexes_of_elems[i] -= 1

print partion_indexes_of_elems

indexes_of_cells, indexes_of_boundary_facets = indexesOfCellsAndFacets(
	map(dimensionOfElem, elems) )

at = lambda array, indexes : [array[i] for i in indexes] 

boundary_facets = at(elems, indexes_of_boundary_facets)
physical_indexes_of_facets = at(physical_indexes_of_elems, indexes_of_boundary_facets) 

cells = at(elems, indexes_of_cells)

partion_indexes_of_cells = at(partion_indexes_of_elems, indexes_of_cells) 
physical_indexes_of_cells = at(physical_indexes_of_elems, indexes_of_cells) 

key_to_physical_index = keyToPhysicalIndex( 
	boundary_facets, physical_indexes_of_facets, lambda elm : keyOfElem( elm, len(nodes) ) )

table = {}
facets = []
neighbor_cells = []
physical_indexes_of_facets = []
partitions_of_facets = []

neighbors_of_cells = []

for i, (cell, physical_index_of_cell, partition_index_of_cell) in enumerate(zip(cells, physical_indexes_of_cells, partion_indexes_of_cells)):
	neighbors_of_cells.append([])
	for facet in facetsOfElem(cell):
		key = keyOfElem(facet, len(nodes))
		if key in table:
			index_of_facet = table[ key ]
			if len(neighbor_cells[ index_of_facet ]) != 1:
				raise
			neighbor_cell = neighbor_cells[ index_of_facet ][0]
			neighbors_of_cells[ neighbor_cell ].append(i)
			neighbors_of_cells[ i ].append(neighbor_cell)
			neighbor_cells[ index_of_facet ].append(i)
			if partitions_of_facets[index_of_facet][0] != partition_index_of_cell:
				partitions_of_facets[ index_of_facet ].append(partition_index_of_cell)

		else:
			table[ key ] = len(facets)
			facets.append(facet)
			neighbor_cells.append([i])
			if key in key_to_physical_index:
				physical_indexes_of_facets.append(key_to_physical_index[key])
			else:
				physical_indexes_of_facets.append(physical_index_of_cell)
			partitions_of_facets.append([partition_index_of_cell])

writeNode = writeNode_text
writeCell = writeCell_text
writeFacet = writeFacet_text

with open(sys.argv[-1], 'wb') as kshfile:

	kshfile.writelines(["<nodes len=\"%d\">\n" % len(nodes)])
	for node in nodes:
		writeNode(kshfile, node)
	kshfile.writelines(["</nodes>\n"])

	kshfile.writelines(["<cells len=\"%d\">\n" % len(cells)])
	for cell, neighbor, partition_index, physical_index in zip(cells, neighbors_of_cells,
		partion_indexes_of_cells, physical_indexes_of_cells ):
		writeCell(kshfile, cell, neighbor, partition_index, physical_index)
	kshfile.writelines(["</cells>\n"])

	kshfile.writelines(["<facets len=\"%d\">\n" % len(facets)])
	for facet, neighbor, physical_index, partition_index in zip(facets, neighbor_cells,
		physical_indexes_of_facets, partitions_of_facets):
		writeFacet(kshfile, facet, neighbor, partition_index, physical_index)
	kshfile.writelines(["</facets>\n"])

	kshfile.writelines(["<physical len=\"%d\">\n" % len(physical_names)])
	for index, value in physical_names.iteritems():
		kshfile.writelines(["%d %s\n" % (index, value)])
	kshfile.writelines(["</physical>\n"])
	
