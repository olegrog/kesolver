#!/usr/bin/env python

import sys, math, numpy, array, struct, itertools
import json
import re

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

dim = (2*rad, rad) 
x = numpy.fromfunction(lambda i, j: (i + 0.5 - rad) * cut / rad, dim, dtype=float)
r = numpy.fromfunction(lambda i, j: (j + 0.5) * cut / rad, dim, dtype=float)

# calculate lambda

with open(maoutfilename, 'rb') as fd:
    size = int(fd.readline())
    theta, section = [], []
    for i in range(size):
        theta.append(float(fd.readline()))
    for i in range(size):
        section.append(float(fd.readline()))

#costheta = numpy.cos(theta)
#section = numpy.array(section)
#lamd = - math.pi / 2 * numpy.trapz(x=costheta, y=(section * (1 - sqr(costheta))), dx=1./size)

section = numpy.array(section)
lamd = numpy.trapz(x=theta, y=(section * cube(numpy.sin(theta))), dx=1./size)

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

        f = numpy.zeros_like(x)
        a = array.array('d')
        circl = x*x + r*r < cut*cut
        size = numpy.sum(circl)
        a.fromfile(fd, size)

        f[circl] = numpy.array(a)

# plot function calculated by projection method

    f = f / r

# plot function calculated by theoretical formula

    file_i = int(re.search(r"(\d+)[^/]*$", ffilename).groups()[0])
    t = file_i * timestep / math.sqrt(2) / math.pi

    xi = 0.4
    tau = 1. - xi * math.exp( - lamd * t )
    e = 0.5 * ( sqr(x) + sqr(r) ) / tau

    g = 1. / math.sqrt(cube(2 * math.pi * tau)) * numpy.exp(-e) * ( 1 + (1 - tau) / tau * (e - 1.5) )

    sum_f = numpy.sum( r * f * sqr( sqr(x) + sqr(r) ) )
    sum_g = numpy.sum( r * g * sqr( sqr(x) + sqr(r) ) )

    time.append(t)
    mom_f.append(sum_f)
    mom_g.append(sum_g)

plt.plot(time, mom_f, 'b')
plt.plot(time, mom_g, 'g')

plt.show()
