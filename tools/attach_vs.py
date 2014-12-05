#!/usr/bin/env python
# attach v,vs,vvs to kep-file
# rely on vector `q` in section `gas`
# Usage: ./attach_vs <kep-file> <new-kep-file>

import sys
import json
import numpy as np

from base64 import b64encode

def vecToStr(l):
    return '(' + ' '.join( map(str, l) ) + ')'

def h1(cut, q, R):
    return NaN if q==1 else cut*(q-1)/(q**R-1)

def i2h(i, q, R):
    i = np.abs(i-R+.5)+.5
    return q**(i-1)

def i2xi(i, q, R):
    sgn = np.sign(i-R+.5)
    i = np.abs(i-R+.5)+.5
    return np.NaN if q==1 else sgn*(q**i+q**(i-1)-2)/(q-1)/2

if __name__ == "__main__":
    
    with open(sys.argv[1], 'r') as f:
        data = json.load(f)

    # input data
    gas = data['gas']
    rad = gas['rad']
    cut = gas['cut']
    q = map(float, gas.get("q", "( 0. 0. 0. )")[1:-1].split())

    # supplementary
    Hxi, Hh = np.polynomial.hermite.hermgauss(2*rad)
    Hh *= np.sqrt(2)/np.exp(-Hxi**2)
    Hxi *= np.sqrt(2)
    idx = np.arange(2*rad)
    Ch = lambda q: i2h(idx, q, rad) * h1(cut, q, rad)
    Cxi = lambda q: i2xi(idx, q, rad) * h1(cut, q, rad)
    switch = lambda q, H, C: H if q==1 else C(q)
    h, xi = np.meshgrid(idx.astype(np.float32), q)

    for (i,q) in enumerate(q):
        h[i] = switch(q, Hh, Ch)
        xi[i] = switch(q, Hxi, Cxi)
        print xi[i,0] - h[i,0]/2, xi[i,2*rad-1] + h[i,2*rad-1]/2

    # output data
    gas['type'] = "Rect"
    gas['v'] = vecToStr([0,0,0])
    gas['vs'] = b64encode(xi.flatten())
    gas['vvs'] = b64encode(h.flatten())

    # open resulting .kep file and dump the modified data 
    with open(sys.argv[2], 'w') as f:
        json.dump(data, f, indent=2, sort_keys=True)
 
