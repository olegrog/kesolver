#!/usr/bin/env python

import sys, math, array, struct, itertools
import json
import re
import numpy as np

import matplotlib.pyplot as plt

from scipy.integrate import quad

from kepy.ximesh import read_ximesh, make_e
from kepy.io import readNodesElems, calc_timestep

def sqr(x):
    if type(x) == np.ndarray:
        return np.dot(x, x)
    else:
        return x * x

def cube(x):
    return x * x * x

keifilename   = sys.argv[1]
maoutfilename = sys.argv[2]
ffilenames    = sys.argv[3:]

# open .kei file
with open(keifilename, 'rb') as fd:
    kei_data = json.load(fd)

symmetry, rad, circl, xyz, vol, r, d3v = read_ximesh(kei_data)

print "size = ", np.sum(circl)

nodes, cells, facets = readNodesElems(keifilename)
timestep = calc_timestep(kei_data, cells, nodes)
print 'time_step = ', timestep

def read_xi(kei_data):
    for k, v in kei_data['initial_conditions'].iteritems():
        if 'xi' in v:
            return v['xi']
    return 0.4 # default value

xi = read_xi(kei_data)
print "xi = ", xi

# calculate lambda

with open(maoutfilename, 'rb') as fd:
    size = int(fd.readline())
    theta, section = [], []
    for i in range(size):
        theta.append(float(fd.readline()))
    for i in range(size):
        section.append(float(fd.readline()))

#costheta = np.cos(theta)
#section = np.array(section)
#lamd = - math.pi / 2 * np.trapz(x=costheta, y=(section * (1 - sqr(costheta))), dx=1./size)

section = np.array(section)
lamd = np.trapz(x=theta, y=(section * cube(np.sin(theta))), dx=1./size)

print lamd

time  = []
mom_f = []
mom_g = []

mom_p = [2, 3, 4, 6]
mom_d = [15, 105, 945, 135135]

#for ffilename in ffilenames:
for i in range(int(ffilenames[0])):
    ffilename = ffilenames[1] % i

    try:
        with open(ffilename, 'rb') as fd:
            data = fd.read(4+4)
            i, size = struct.unpack('=ii', data)

            f = np.zeros_like(r)
            a = array.array('d')
            size = np.sum(circl)
            a.fromfile(fd, size)

            f[circl] = np.array(a)
    except IOError: 
        continue

    f = f / vol

    e = make_e(symmetry, xyz)
#    print "e = ", np.sum( vol[circl] * f[circl] * e[circl] )
#    print "n = ", np.sum( vol[circl] * f[circl] ) 
    sum2_f = np.sum( vol[circl] * f[circl] * e[circl] ) / \
         3 / np.sum( vol[circl] * f[circl] ) 
    m2 = 1 / sum2_f
    m1 = 1 / (np.sum( vol[circl] * f[circl] ) * d3v)
#    print m1, m2

    file_i = int(re.search(r"(\d+)[^/]*$", ffilename).groups()[0])
    t = file_i * timestep / 2 / math.sqrt(2) / math.pi

    tau = 1. - xi * math.exp( - lamd * t )

    sum_f = [m1 * np.sum( vol * f * (e * m2)**p ) * d3v / d for p, d in zip(mom_p, mom_d)]
#    sum_g = [np.sum( vol * g * e**p ) * d3v / d for p, d in zip(mom_p, mom_d)]
    sum_g = [ tau**(p-1) * (1 + (p-1) * xi * math.exp( - lamd * t )) for p in mom_p ]

    time.append(t)
    mom_f.append(sum_f)
    mom_g.append(sum_g)

mom_f = zip(*mom_f)
mom_g = zip(*mom_g)


d = [np.abs(np.array(mf) - np.array(mg))
        for mf, mg in zip(mom_f, mom_g)]
err = [np.trapz(x=time, y=y) for y in d]
print err

if False:
    for y in d:
        plt.plot(time, y)
else:
    for mf, mg in zip(mom_f, mom_g):
        plt.plot(time, mf, 'b')
        plt.plot(time, mg, 'g')

with open("outf.txt", "w") as out:
    for u, v in zip(time, zip(*mom_f)):
        out.writelines(str(u) + ' ' + ' '.join(map(str, v)) + '\n')

with open("outg.txt", "w") as out:
    for u, v in zip(time, zip(*mom_g)):
        out.writelines(str(u) + ' ' + ' '.join(map(str, v)) + '\n')

with open("outd.txt", "w") as out:
    for u, v in zip(time, zip(*d)):
        out.writelines(str(u) + ' ' + ' '.join(map(str, v)) + '\n')

plt.show()
