#!/usr/bin/env python

import sys, os, json, errno, re
import out2
from kepy.io import readNodesElems
from collections import namedtuple
import numpy as np

geom_params = {
    'scalar': { 'shift': 1, 'fmt': '%.6g' },
    'vector': { 'shift': 3, 'fmt': '(%.6g %.6g %.6g)' }
}
kes2foam_type = {
    'scalar': {
        'empty': 'empty',
        'mirror': 'symmetryPlane',
        'diffusion': 'extrapolatedGradient' # or use zeroGradient
    },
    'vector': {
        'empty': 'empty',
        'mirror': 'symmetryPlane',
        'diffusion': 'slip'
    },
}
kes2foam_type2 = {
    'empty': 'empty',
    'mirror': 'symmetryPlane',
    'diffusion': 'wall'
}

Field = namedtuple('Field', ['foam_name', 'kes_name', 'gtype', 'dimensions', 'pos'])
fields = [
    Field('rho', '', 'scalar', [1,-3,0,0,0,0,0], 0),
    Field('U',  'u', 'vector', [0,1,-1,0,0,0,0], 1),
    Field('T',  'T', 'scalar', [0,0,0,1,0,0,0],  4)
]
aux_fields = [ 'p' ]

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

class List:
    def __init__(self, out, name, shift=0):
        out.write('%s%s\n%s(\n' % (' '*shift, name, ' '*shift))
        self.out = out
        self.shift = shift
    def __enter__(self):
        return self
    def append(self, name):
        return Dict(self.out, name, self.shift + 4)
    def __exit__(self, type, value, traceback):
        self.out.write('%s);\n' % (' '*self.shift))

class Dict:
    def __init__(self, out, name, shift=0):
        out.write('%s%s\n%s{\n' % (' '*shift, name, ' '*shift))
        self.out = out
        self.shift = shift
    def __enter__(self):
        return self
    def print_dict(self, dict):
        print_dict(self.out, dict, shift=self.shift + 4)
    def child(self, name):
        return Dict(self.out, name, self.shift + 4)
    def list(self, name):
        return List(self.out, name, self.shift + 4)
    def nlist(self, key, gtype, size):
        return Nonuniform_list(self.out, key, gtype, size, self.shift + 4)
    def __exit__(self, type, value, traceback):
        self.out.write('%s}\n' % (' '*self.shift))

def print_internal_field(out, field, data):
    p = geom_params[field.gtype]
    with Nonuniform_list(out, 'internalField', field.gtype, len(data[0])) as nl:
        nl.dump_data(data[field.pos:field.pos+p['shift']], p['fmt'])

def print_boundary_field(out, field, facets, bc):
    with Dict(out, 'boundaryField') as bfield:
        patches = set(map(lambda f: f.split(':')[0], bc.keys()))
        patches.add('defaultFaces')
        for patch in patches:
            with bfield.child(patch) as pfield:
                foam_type = kes2foam_type[field.gtype][bc[patch]['type']]
                pfield.print_dict({ 'type': foam_type })
                if foam_type == 'fixedValue':
                    patch_facets = filter(lambda f: f.phys_index.split(':')[0] == patch, facets)
                    values = map(lambda f: bc[f.phys_index][field.kes_name], patch_facets)
                    with pfield.nlist('value', field.gtype, len(patch_facets)) as nl:
                        nl.dump_data(values, '%s')
                if foam_type == 'extrapolatedGradient':
                    pfield.print_dict({ 'value': '$internalField' })

def print_field(field, time, data, facets, bc):
    with open('%s/%s' % (time, field.foam_name), 'w') as f:
        print_header(f, 'vol%sField' % field.gtype.title(), time, field.foam_name)
        print_dict(f, { 'dimensions': str(field.dimensions).replace(',', '') })
        print_internal_field(f, field, data)
        print_boundary_field(f, field, facets, bc)

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
            'startFrom':        'startTime',
            'startTime':        0,
            'stopAt':           'endTime',
            'endTime':          kei['num_steps'] - 1,
            'deltaT':           kei['printer']['savemacro'],
            'writeInterval':    kei['printer']['savemacro'],
            'writePrecision':   7,
            'writeFormat':      'ascii'
        })
        with Dict(f, 'functions') as funs:
            with funs.child('fieldAverage') as fa:
                fa.print_dict({
                    'type':                 'fieldAverage',
                    'functionObjectLibs':   '( "libfieldFunctionObjects.so" )',
                    'enabled':              'true',
                    'outputControl':        'timeStep'
                })
                with fa.list('fields') as fs:
                    for name in map(lambda f: f.foam_name, fields) + aux_fields:
                        with fs.append(name) as f:
                            f.print_dict({
                                'mean':         'on',
                                'prime2Mean':   'on',
                                'base':         'iteration',
                                'window':       0.15*(kei['num_steps'] - 1)/kei['printer']['savemacro']
                            })
    with open('system/fvSchemes', 'w') as f:
        print_header(f, 'dictionary', 'system', 'fvSchemes')
        with Dict(f, 'interpolationSchemes') as field:
            pass
        with Dict(f, 'divSchemes') as field:
            pass
        with Dict(f, 'gradSchemes') as field:
            field.print_dict({ 'default': 'Gauss' })
        with Dict(f, 'laplacianSchemes') as field:
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
        for field in fields:
            print_field(field, time, data, facets, bc)

