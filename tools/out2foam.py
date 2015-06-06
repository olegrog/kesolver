#!/usr/bin/env python

import sys, os, json, errno, re
import out2
from kepy.io import readNodesElems
import numpy as np

geom_params = {
    'scalar': { 'shift': 1, 'fmt': '%.6g' },
    'vector': { 'shift': 3, 'fmt': '(%.6g %.6g %.6g)' }
}
foam2kes_macro = { 'U': 'u', 'T': 'T' }
kes2foam_type = {
    'empty': 'empty',
    'mirror': 'symmetryPlane',
    'diffusion': 'fixedValue'
}
kes2foam_type2 = {
    'empty': 'empty',
    'mirror': 'symmetryPlane',
    'diffusion': 'wall'
}

def mkdir_p(path):
    try:
        os.makedirs(path)
        print 'Create %s' % path
        return True
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            return False
        else: raise

def print_header(out, fclass, flocation, fobject):
    out.write('\
/*--------------------------------*- C++ -*----------------------------------*\\\n\
| =========                 |                                                 |\n\
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n\
|  \\\\    /   O peration     | Version:  2.3.0                                 |\n\
|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |\n\
|    \\\\/     M anipulation  |                                                 |\n\
\*---------------------------------------------------------------------------*/\n\
FoamFile\n\
{\n\
    version     2.0;\n\
    format      ascii;\n\
    class       %s;\n\
    location    "%s";\n\
    object      %s;\n\
}\n\
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n\
' % (fclass, flocation, fobject))

def print_dict(out, dict, shift=0):
    for key, value in dict.iteritems():
        out.write('%s%-20s%s;\n' % (' '*shift, key, str(value)))

class Nonuniform_list:
    def __init__(self, out, key, gtype, size, shift=0):
        out.write('%s%-20snonuniform List<%s>\n' % (' '*shift, key, gtype))
        out.write('%d\n(\n' % size)
        self.out = out
    def __enter__(self):
        return self
    def dump_data(self, data, fmt):
        np.savetxt(self.out, np.transpose(data), fmt=fmt)
    def __exit__(self, type, value, traceback):
        self.out.write(')\n;\n')

class Record:
    def __init__(self, out, name, shift=0):
        out.write('%s%s\n%s{\n' % (' '*shift, name, ' '*shift))
        self.out = out
        self.shift = shift
    def __enter__(self):
        return self
    def print_dict(self, dict):
        print_dict(self.out, dict, shift=self.shift + 4)
    def child(self, name):
        return Record(self.out, name, self.shift + 4)
    def nlist(self, key, gtype, size):
        return Nonuniform_list(self.out, key, gtype, size, self.shift + 4)
    def __exit__(self, type, value, traceback):
        self.out.write('%s}\n' % (' '*self.shift))

def print_internal_field(out, gtype, data, pos):
    p = geom_params[gtype]
    with Nonuniform_list(out, 'internalField', gtype, len(data[0])) as nl:
        nl.dump_data(data[pos:pos+p['shift']], p['fmt'])

def print_boundary_field(out, gtype, kes_name, facets, bc):
    with Record(out, 'boundaryField') as bfield:
        boundary = filter(lambda f: f.phys_index != 'volume', facets)
        patches = set(map(lambda f: f.phys_index.split(':')[0], boundary))
        patches.add('defaultFaces')
        for patch in patches:
            with bfield.child(patch) as pfield:
                foam_type = kes2foam_type[bc[patch]['type']]
                pfield.print_dict({ 'type': foam_type })
                if foam_type == 'fixedValue':
                    patch_facets = filter(lambda f: f.phys_index.split(':')[0] == patch, boundary)
                    values = map(lambda f: bc[f.phys_index][kes_name], patch_facets)
                    with pfield.nlist('value', gtype, len(patch_facets)) as nl:
                        nl.dump_data(values, '%s')

def print_macro(foam_name, gtype, dim, data, pos, facets, bc):
    with open('%s/%s' % (time, foam_name), 'w') as f:
        print_header(f, 'vol%sField' % gtype.title(), time, foam_name)
        print_dict(f, { 'dimensions': str(dim).replace(',', '') })
        print_internal_field(f, gtype, data, pos)
        print_boundary_field(f, gtype, foam2kes_macro[foam_name], facets, bc)

#############################
### Start of instructions ###
#############################

_, case = sys.argv
nodes, cells, facets = readNodesElems('%s.kei' % case)
with open('%s.kei' % case, 'r') as f:
    kei = json.load(f)
bc = kei['boundary_conditions']
bc['defaultFaces'] = { 'type': 'empty' }

if mkdir_p('system'):
    with open('system/controlDict', 'w') as f:
        print_header(f, 'dictionary', 'system', 'controlDict')
        print_dict(f, {
            'deltaT':           kei['printer']['savemacro'],
            'writeInterval':    kei['printer']['savemacro'],
            'writeFormat':      'ascii'
        })
    with open('system/fvSchemes', 'w') as f:
        print_header(f, 'dictionary', 'system', 'fvSchemes')
        with Record(f, 'interpolationSchemes') as field:
            pass
        with Record(f, 'divSchemes') as field:
            pass
        with Record(f, 'gradSchemes') as field:
            pass
        with Record(f, 'laplacianSchemes') as field:
            pass
    with open('system/fvSolution', 'w') as f:
        print_header(f, 'dictionary', 'system', 'fvSolution')

if mkdir_p('constant'):
    os.system('gmshToFoam %s.msh > /dev/null' % case)
    lines = []
    cb, rb = 0, 0
    with open('constant/polyMesh/boundary') as f:
        for line in f:
            rb = {
                '(': lambda x: x + 1,
                ')': lambda x: x - 1
            }.get(line.strip(), lambda x: x)(rb)
            cb = {
                '{': lambda x: x + 1,
                '}': lambda x: x - 1
            }.get(line.strip(), lambda x: x)(cb)
            words = re.findall('\w+', line)
            if len(words):
                if rb == 1 and cb == 0:
                    patch = words[0]
                if rb == 1 and cb == 1 and words[0] == 'type':
                    line = line.replace(words[1], kes2foam_type2[bc[patch]['type']])
            lines.append(line)
    with open('constant/polyMesh/boundary', 'w') as f:
        for line in lines:
            f.write(line)

dirname = kei['printer']['dir']
pattern = lambda f: f.endswith(kei['printer']['file'].split('.')[-1]) 
for filename in filter(pattern, os.listdir(dirname)):
    time = re.search('[0-9]+', filename).group(0)
    if mkdir_p(time):
        data = out2.readMacros(os.path.join(dirname, filename), len(cells))
        print_macro('T', 'scalar', [0,0,0,1,0,0,0], data, 4, facets, bc)
        print_macro('U', 'vector', [0,1,-1,0,0,0,0], data, 1, facets, bc)

