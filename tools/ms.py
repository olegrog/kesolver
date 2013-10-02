#!/usr/bin/python

import numpy, re, json
from kepy.io import readNodesElems, calc_timestep
from pylab import *

def diff(x, t):
    dx = [(x[1]-x[0])/(t[1]-t[0])]
    for i in range(1, len(x)-1):
        dx.append( (x[i+1]-x[i-1])/(t[i+1]-t[i-1]) )
    dx.append( (x[-1]-x[-2])/(t[-1]-t[-2]) )
    return dx

keifilename = sys.argv[1]

# open .kei file
with open(keifilename, 'rb') as fd:
    kei_data = json.load(fd)

nodes, cells, facets = readNodesElems(keifilename)
time_step = calc_timestep(kei_data, cells, nodes)
print 'time_step = ', time_step

volumes = {}
for filename in sys.argv[2:]:
    f = file(filename)

    file_i = int(re.search(r"(\d+)[^/]*$", filename).groups()[0])
    time = file_i
    time = time_step * file_i
    print filename, file_i, time_step

    for line in f.readlines():
        words = line.split()
        volume = words[0]
        data = map(float, words[1:])
        if volume in volumes:
            volumes[volume].append( ( time, data ) )
        else:
            volumes[volume] = [ ( time, data ) ]

number_of_volumes = len(volumes)
subplottype = 2*number_of_volumes * 100 + 10 + 1;
colors = ['r', 'b', 'g']
styles = ['-', '--', '-.']

print volumes.keys()
for key, values, color in zip(volumes.keys(), volumes.values(), colors):
    subplot(subplottype)
    values.sort(key = lambda x: x[0])
    ts, ms = zip( *values )
    ms = zip(*ms)
    
    for ys, s in zip(ms, styles):
        plot(ts, ys, color+s)
    grid(True)
    ylabel(str(key))

    with open(str(key)+".ms", "w") as f:
        for t, m in zip(ts, zip(*ms)):
            f.writelines( "%.15f %s\n" % (t, ' '.join(map(str, m))) )

    subplottype += 1

    dms = [diff(ys, ts) for ys in ms]

    subplot(subplottype)
    for ys, s in zip(dms, styles):
        plot(ts, ys, color+s)
    grid(True)

    subplottype += 1

show()

