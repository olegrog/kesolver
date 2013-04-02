#!/usr/bin/env python

import sys, math, array, struct, itertools
import json
import re
import numpy as np

import matplotlib.pyplot as plt

from scipy.integrate import quad

from kepy.ximesh import read_ximesh, make_e

def sqr(x):
    return x * x

def cube(x):
    return x * x * x

keifilename   = sys.argv[1]
maoutfilename = sys.argv[2]
timestep      = float(sys.argv[3])
ffilenames    = sys.argv[4:]

# open .kei file
with open(keifilename, 'rb') as fd:
    data = json.load(fd)

symmetry, rad, cut, xyz, vol, r, d3v = read_ximesh(data)

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

mom_p = [2, 4, 6]
mom_d = [15, 945, 135135]

#for ffilename in ffilenames:
for i in range(int(ffilenames[0])):
    ffilename = ffilenames[1] % i

    try:
        with open(ffilename, 'rb') as fd:
            data = fd.read(4+4)
            i, size = struct.unpack('=ii', data)

            f = np.zeros_like(r)
            a = array.array('d')
            circl = r < cut*cut
            size = np.sum(circl)
            a.fromfile(fd, size)

            f[circl] = np.array(a)
    except IOError: 
        continue

    f = f / vol

    file_i = int(re.search(r"(\d+)[^/]*$", ffilename).groups()[0])
    t = file_i * timestep / math.sqrt(2) / math.pi

    xi = 0.4
    tau = 1. - xi * math.exp( - lamd * t )

    e = make_e(symmetry, xyz)
    e1 = 0.5 * e / tau

    g = 1. / math.sqrt(cube(2 * math.pi * tau)) * np.exp(-e1) * ( 1 + (1 - tau) / tau * (e1 - 1.5) )

    g /= np.sum(vol*g) / np.sum(vol*f)

    sum_f = [np.sum( vol * f * e**p ) * d3v / d for p, d in zip(mom_p, mom_d)]
    sum_g = [np.sum( vol * g * e**p ) * d3v / d for p, d in zip(mom_p, mom_d)]

    time.append(t)
    mom_f.append(sum_f)
    mom_g.append(sum_g)

mom_f = zip(*mom_f)
mom_g = zip(*mom_g)


d = [np.abs(np.array(mf) - np.array(mg)) / np.array(mg)
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
