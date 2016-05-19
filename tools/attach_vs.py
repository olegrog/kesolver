#!/usr/bin/env python
# attach v,vs,vvs to kep-file
# rely on vector `q` and `N` in section `gas`
# Usage: ./attach_vs <kep-file> <new-kep-file>

import sys, json
import numpy as np

from base64 import b64encode
from collections import namedtuple

Grid = namedtuple('Grid', ['i2h','i2xi'])

grids = {
    'uniform': Grid(
        i2h =  lambda q, cut, N: cut/N + idx(N)*0,
        i2xi = lambda q, cut, N: idx(N)*cut/N
    ),
    'geometric': Grid(
        i2h =  lambda q, cut, N: h1(q, cut, N) * symm_h(q**np.arange(N)),
        i2xi = lambda q, cut, N: h1(q, cut, N) * (1 if q==1 else symm_xi((q**np.arange(N+1)-1)/(q-1)))
    ),
    'hermite': Grid(
        i2h =  lambda q, cut, N: hermite(cut, N)[0],
        i2xi = lambda q, cut, N: hermite(cut, N)[1]
    ),
    'quadratic': Grid(
        i2h =  lambda q, cut, N: quadratic(q, cut, N)[0],
        i2xi = lambda q, cut, N: quadratic(q, cut, N)[1]
    )
}

def hermite(cut, N):
    H_xi, H_h = np.polynomial.hermite.hermgauss(2*N)
    H_h /= np.exp(-H_xi**2)
    C = 2*cut / np.sum(H_h)
    H_xi *= C
    H_h *= C
    return H_h, H_xi

def quadratic(q, cut, N):
    p=2
    if cut-N*q < 0:
       raise NameError('q = %g is too big' % q)
    A = (cut-N*q) / np.sum(np.arange(N)**p)
    h = q + A*np.arange(N)**p
    X = np.append(0, np.cumsum(h))
    return symm_h(h), symm_xi(X)

symm_h = lambda x: np.hstack((x[::-1], x))
semi_sum = lambda x: .5*(x[:-1] + x[1:])
symm_xi = lambda x: np.hstack((-semi_sum(x)[::-1], semi_sum(x)))
h1 = lambda q, cut, N: cut/N if q==1 else cut*(q-1)/(q**N-1)
idx = lambda N: np.arange(2*N) - N + .5
vecToStr = lambda l: '(' + ' '.join( map(str, l) ) + ')'

def extend_grid(xi, h, cut, N):
    idx = (N-xi.size/2)
    if idx:
        xi_new, h_new = np.full(2*N, 10*cut)*np.sign(np.arange(-N,N)), np.full(2*N, 0)
        xi_new[idx:-idx] = xi
        h_new[idx:-idx] = h
        return xi_new, h_new
    else:
        return xi, h

def make_vector(x):
    if len(x) == 1:
        x = [ x[0], x[0], x[0] ]
    return x
    

if __name__ == '__main__':
    
    with open(sys.argv[1], 'r') as f:
        data = json.load(f)

    # input data
    gas = data['gas']
    cut = gas['cut']
    grid = make_vector(gas.get('nonuniform', '( uniform )')[1:-1].split())
    q = make_vector(map(float, gas.get('q', '( 0. 0. 0. )')[1:-1].split()))
    N = make_vector(map(int, gas.get('N_R', '( 0 0 0 )')[1:-1].split()))
    rad = max(gas['rad'], max(N))
    if np.dot(N, N) == 0:
        print 'No velocity grid has been attached.'
    else:
        h, xi = np.ndarray((3,2*rad)), np.ndarray((3,2*rad))

        for i in xrange(3):
            xi[i], h[i] = extend_grid(
                grids[grid[i]].i2xi(q[i], cut, N[i]),
                grids[grid[i]].i2h (q[i], cut, N[i]),
                cut, rad)
            # sqrt(2) Tcheremissine --> Sone
            print '--- Axis %s' % chr(ord('X')+i)
            print 'xi =', xi[i]/np.sqrt(2)
            print 'h =', h[i]/np.sqrt(2)

        # output data
        gas['v'] = vecToStr([0,0,0])
        gas['vs'] = b64encode(xi.flatten())
        gas['vvs'] = b64encode(h.flatten())
        gas['rad'] = rad

    # open resulting .kep file and dump the modified data 
    with open(sys.argv[2], 'w') as f:
        json.dump(data, f, indent=2, sort_keys=True)
 
