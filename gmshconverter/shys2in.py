#!/usr/bin/env python

import sys, string
import numpy, math
import re


def readNode_text(str):
    return tuple( [float(w) for w in str.split()] )
readNode = readNode_text

def readCell(str):
    ws = str.split()
    return int(ws[-1]), ws[:-1]

def numberOfNodes(elm_type):
    number_of_nodes = [-1, 2, 3, 4, 4, 8, 6, 5, 3, 6, 9, 10, 27, 18, 14, 1]
    return number_of_nodes[elm_type]

def facetNeigbors(facet):
    i = numberOfNodes(int(facet[1]))+2
    number_of_facets =  int(facet[i])
    return facet[i+1:i+1+number_of_facets]

def facetData(facet):
    i = numberOfNodes(int(facet[1]))+2
    number_of_facets = int(facet[i])
    j = i+1+number_of_facets
    return facet[:i], facet[j:]

def center(cell, nodes):
    c = numpy.zeros(3)
    elm_type = int(cell[1])
    number_of_nodes = numberOfNodes(elm_type)
    for w in cell[2:2+number_of_nodes]:
        c += nodes[ int(w) ]
    return c / number_of_nodes

def newPhysicalData(elm, data, nodes, exceptions):
    c = center(elm, nodes)
    x = c[0]; y = c[1]; z = c[2]
    newdata = []
    for w in data[1:]:
        if w in exceptions:
            newdata.append(w)
        else:
            newdata.append(str(eval(w)))
    return newdata

def transformToOriginalOrder(elms_by):
    new_elms = {}
    for physical_index, elms in zip(elms_by.keys(),
            elms_by.values()):
        for elm in elms:
            original_i = int(elm[0])
            new_str = ' '.join( elm[1:] ) + ' ' + str(physical_index) + '\n'
            if original_i in new_elms:
                new_elms[ original_i ].append( new_str )
            else:
                new_elms[ original_i ] = [new_str]  
    return new_elms

def reverse(facet):
    i = numberOfNodes(int(facet[1]))
    newfacet = [n for n in facet]
    newfacet[2:2+i] = [ n for n in reversed(facet[2:2+i]) ]
    return newfacet

def writeElms(kinfile, elms_by, name):
    number_of_elms = sum([len(elms) for elms in elms_by.values()])
    kinfile.writelines( "<%s len=\"%d\">\n" % (name, number_of_elms) )
    for elms in elms_by.values():
        for elm in elms:
            kinfile.writelines(elm) 
    kinfile.writelines( "</%s>\n" % name )

def readElms(kshfile, name, elms_by_physical_index):
    line = " "
    endline = "</%s>\n" % name
    i = 0
    while line != True:
        line = kshfile.readline()
        if line == endline:
            break
        physical_index, data = readCell(line)
        i += 1
        if physical_index in elms_by_physical_index:
            elms_by_physical_index[ physical_index ].append([str(i)] + data)
        else:
            elms_by_physical_index[ physical_index ] = [ [str(i)] + data ]

def writeNodes(kinfile, nodes, m):
    kinfile.writelines( "<nodes len=\"%d\">\n" % (len(nodes)) )
    for node in nodes:
        kinfile.writelines( "%f %f %f\n" % tuple(node*m))
    kinfile.writelines( "</nodes>\n" )

def isXmlStart(key, line):
    return re.search(r"^<" + key + r"\b", line)
def isXmlEnd(key, line):
    return re.search(r"</" + key + r">$", line)
def getXmlValue(s, key):
    print s, key
    return re.search(key + r"\s*=\s*\"([^\"]+)\"", s).groups()[0]

def frac(x):
    return x - math.floor(x)

def zigzag(x):
    return abs(2*frac(x/2)-1)

def snake(x, y, R, l):
    if 0 < y < l:
        i = int((x + R) / 2 / R)
        if (i % 2) == 0:
            return y / l
        else:
            return 1 - y / l
    elif y >= l:
        x1 = x - R
        while not (-2*R <= x1 <= 2*R):
            if x1 < -2*R:
                x1 += 4*R
            else:
                x1 -= 4*R
        cos = x1 / math.sqrt(x1*x1 + (y-l)**2)
        phi = math.acos(cos) 
        return phi / math.pi
    else:
        x1 = x - 3*R
        while not (-2*R <= x1 <= 2*R):
            if x1 < -2*R:
                x1 += 4*R
            else:
                x1 -= 4*R
        cos = x1 / math.sqrt(x1*x1 + y**2)
        phi = math.acos(cos) 
        return phi / math.pi

with open(sys.argv[-1], 'wb') as kinfile:

    kinfile.writelines("<?xml version=\"1.0\"?>\n")
    kinfile.writelines("<kin>\n")

    physical_data = []
    mult = 1.
    with open(sys.argv[-3], 'rb') as kysfile:
        line = " "
        while line:
            line = kysfile.readline()
            if isXmlStart("physical", line):
                while True:
                    line = kysfile.readline()
                    if isXmlEnd("physical", line):
                        break
                    words = line.split()
                    physical_data.append(words)
            elif isXmlStart("mult", line):
                mult = float(getXmlValue(line, "value"))
            else:
                kinfile.writelines([line])

    physical_name_to_index = {}
    nodes = []
    cells_by_physical_index = {}
    number_of_cells = 0
    facets_by_physical_index = {}
    with open(sys.argv[-2], 'rb') as kshfile:
        line = " "
        while line:
            line = kshfile.readline()
            if isXmlStart("physical", line):
                while line != True:
                    line = kshfile.readline()
                    if isXmlEnd("physical", line):
                        break
                    words = line.split()
                    physical_name_to_index[ words[1] ] = int(words[0])
            elif isXmlStart("nodes", line):
                while line != True:
                    line = kshfile.readline()
                    if isXmlEnd("nodes", line):
                        break
                    nodes.append(numpy.array(readNode(line)))               
            elif isXmlStart("cells", line):
                readElms(kshfile, "cells", cells_by_physical_index)
            elif isXmlStart("facets", line):
                readElms(kshfile, "facets", facets_by_physical_index)
            else:
                kinfile.writelines([line])

    print "mult = %f\n" % mult

    facetsymbols = set(['m', 'w', 'l', 'f'])
    exceptions = set(['m', 'w', 'l', 'f', '(', ')', '[', ']'])
    print facetsymbols, exceptions

    j = 0
    for data in physical_data:
        if data[1] in facetsymbols:
            index =  physical_name_to_index[data[0]]
            for i, facet in enumerate(facets_by_physical_index[index]):
                if len(facetNeigbors(facet)) > 1:
                    j += 1
                    newfacets = []
                    for neigbor in facetNeigbors(facet):
                        before, after = facetData(facet)
                        newfacets.append(before + ['1', str(neigbor)] + after)
                    facets_by_physical_index[index][i] = newfacets[0]
                    for f in newfacets[1:]:
                        facets_by_physical_index[index].append(reverse(f))
        else:
            index =  physical_name_to_index[data[0]]
            i = 0
            while i < len(facets_by_physical_index[index]):
                if len(facetNeigbors(facets_by_physical_index[index][i])) != 2:
                    del facets_by_physical_index[index][i]
                else:
                    i += 1
            
    print j

    new_physical_data = []
    new_index = 0
    index_by_data = {}

    for data in physical_data:
        physical_name = data[0]
        index = physical_name_to_index[ physical_name ]
        print index, data
        data[0] = str(index)
        new_physical_data.append(data)
        try:
            for value in data[1:]:
                if value not in exceptions:
                    float(value)
        except ValueError:
            if index in cells_by_physical_index:
                for cell in cells_by_physical_index[index]:
                    new_data = newPhysicalData(cell, data, nodes, exceptions)
                    key_data = ' '.join(new_data)
                    if key_data in index_by_data:
                        i = index_by_data[key_data]
                        cells_by_physical_index[i].append(cell)
                    else:
                        new_index -= 1
                        index_by_data[key_data] = new_index
                        i = new_index
                        cells_by_physical_index[i] = [cell]
                        new_physical_data.append([str(i)] + new_data)   
                del cells_by_physical_index[index]
            elif index in facets_by_physical_index:
                for facet in facets_by_physical_index[index]:
                    new_data = newPhysicalData(facet, data, nodes, exceptions)
                    key_data = ' '.join(new_data)
                    if key_data in index_by_data:
                        i = index_by_data[key_data]
                        facets_by_physical_index[i].append(facet)
                    else:
                        new_index -= 1
                        index_by_data[key_data] = new_index
                        i = new_index
                        facets_by_physical_index[i] = [facet]
                        new_physical_data.append([str(i)] + new_data)   
                del facets_by_physical_index[index]

    kinfile.writelines( ["<physical len=\"%d\">\n" % len(new_physical_data)] )
    for data in new_physical_data:
        kinfile.writelines([ ' '.join(data) + '\n' ])   
    kinfile.writelines( ["</physical>\n"] )

    writeNodes(kinfile, nodes, mult)
    writeElms(kinfile, transformToOriginalOrder(cells_by_physical_index), "cells")
    writeElms(kinfile, transformToOriginalOrder(facets_by_physical_index), "facets")

    kinfile.writelines("</kin>\n")
    
