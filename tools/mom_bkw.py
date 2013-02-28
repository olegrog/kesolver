#!/usr/bin/env python

import sys, math, array, struct, itertools
import json
import re
import numpy as np

import matplotlib.pyplot as plt

from scipy.integrate import quad

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

cut = float(data["gas"]["cut"])
rad = int(data["gas"]["rad"])
symmetry = data["gas"]["symmetry"]

print rad

if symmetry == "Cylindrical":
    v = float(data["gas"].get("v", "0."))
    dim = (2*rad, rad) 
    x = np.fromfunction(lambda i, j: v + (i + 0.5 - rad) * cut / rad,
                        (2 * rad, rad), dtype=float)
    y = np.fromfunction(lambda i, j: (j + 0.5) * cut / rad,
                        (2 * rad, rad), dtype=float)
    vol = y
    e = x*x + y*y
    
elif symmetry == "Cartesian":
    v = map(float, data["gas"].get("v", "( 0. 0. 0. )")[1:-1].split() )
    vx, vy, vz = v
    dim = (2*rad, 2*rad, 2*rad)
    x = np.fromfunction(lambda i, j, k: vx + (i + 0.5 - rad) * cut / rad,
                        (2 * rad, 2 * rad, 2 * rad), dtype=float)
    y = np.fromfunction(lambda i, j, k: vy + (j + 0.5 - rad) * cut / rad,
                        (2 * rad, 2 * rad, 2 * rad), dtype=float)
    z = np.fromfunction(lambda i, j, k: vz + (k + 0.5 - rad) * cut / rad,
                        (2 * rad, 2 * rad, 2 * rad), dtype=float)
    vol = np.ones_like(y)
    e = x*x + y*y + z*z

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

#for ffilename in ffilenames:
for i in range(int(ffilenames[0])):
    ffilename = ffilenames[1] % i

    with open(ffilename, 'rb') as fd:
        data = fd.read(4+4)
        i, size = struct.unpack('=ii', data)

        f = np.zeros_like(x)
        a = array.array('d')
        circl = e < cut*cut
        size = np.sum(circl)
        a.fromfile(fd, size)

        f[circl] = np.array(a)

    f = f / vol

    file_i = int(re.search(r"(\d+)[^/]*$", ffilename).groups()[0])
    t = file_i * timestep / math.sqrt(2) / math.pi

    xi = 0.4
    tau = 1. - xi * math.exp( - lamd * t )
    e1 = 0.5 * e / tau

    g = 1. / math.sqrt(cube(2 * math.pi * tau)) * np.exp(-e1) * ( 1 + (1 - tau) / tau * (e1 - 1.5) )

    g /= np.sum(vol*g) / np.sum(vol*f)

    sum_f = np.sum( vol * f * sqr(e) )
    sum_g = np.sum( vol * g * sqr(e) )

    time.append(t)
    mom_f.append(sum_f)
    mom_g.append(sum_g)

#plt.plot(time, mom_f, 'b')
#plt.plot(time, mom_g, 'g')

d = np.abs(np.array(mom_f) - np.array(mom_g)) / np.array(mom_g)
err = np.trapz(x=time, y=d)
print err

plt.plot(time, d, 'g')

with open("out1.txt", "w") as out:
    for u, v in zip(time, mom_g):
        out.writelines(str(u) + ' ' + str(v) + '\n')

with open("out2.txt", "w") as out:
    for u, v in zip(time, d):
        out.writelines(str(u) + ' ' + str(v) + '\n')

plt.show()
