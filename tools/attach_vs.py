#!/usr/bin/env python
# attach v,vs,vvs to kep-file
# rely on vector `q` and `N` in section `gas`
# Usage: ./attach_vs <kep-file> <new-kep-file>

import sys
import json
import numpy as np

from base64 import b64encode

def hermite(N):
    H_xi, H_h = np.polynomial.hermite.hermgauss(2*N)
    H_h *= np.sqrt(2)/np.exp(-H_xi**2)
    H_xi *= np.sqrt(2)
    if np.sum(H_h) < 2*cut:
        C = (2*cut) / np.sum(H_h)
        H_xi *= C
        H_h *= C
    return lambda i: H_xi[i], lambda i: H_h[i]

h1 = lambda q, cut, N: cut/N if q==1 else cut*(q-1)/(q**N-1)
vecToStr = lambda l: '(' + ' '.join( map(str, l) ) + ')'

def i2h(i, q, cut, N):
    j = lambda i: abs(i-N+.5)+.5
    nonuniform = lambda i: h1(q, cut, N)*q**(j(i)-1)
    uniform = lambda i: cut/N + i*0
    return {
        0: hermite(N)[1],       # hermite grid
        -1: uniform             # uniform grid
    }.get(q, nonuniform)(i)     # exp refinement

def i2xi(i, q, cut, N):
    j_y = lambda i: abs(i-N+.5)+.5
    j_x = lambda i: i-N+.5
    sgn = lambda i: np.sign(j_x(i))
    nonuniform = lambda i: sgn(i)*h1(q, cut, N) * (j_y(i)-.5 if q==1 else (q**j_y(i)+q**(j_y(i)-1)-2)/(q-1)/2)
    uniform = lambda i: j_x(i)*cut/N
    return {
        0: hermite(N)[0],       # hermite grid
        -1: uniform             # uniform grid
    }.get(q, nonuniform)(i)     # exp refinement

def extend_grid(xi, h, cut, N):
    idx = (N-xi.size/2)
    if idx:
        xi_new, h_new = np.full(2*N, 10*cut)*np.sign(np.arange(-N,N)), np.full(2*N, 0)
        xi_new[idx:-idx] = xi
        h_new[idx:-idx] = h
        return xi_new, h_new
    else:
        return xi, h


if __name__ == "__main__":
    
    with open(sys.argv[1], 'r') as f:
        data = json.load(f)

    # input data
    gas = data['gas']
    rad = gas['rad']
    cut = gas['cut']
    qi = map(float, gas.get("q", "( 0. 0. 0. )")[1:-1].split())
    Ni = map(int, 2*gas.get("N", "( 0 0 0 )")[1:-1].split())

    h, xi = np.ndarray((3,2*rad)), np.ndarray((3,2*rad))

    for (i, q) in enumerate(qi):
        idx = np.arange(2*Ni[i])
        xi[i], h[i] = extend_grid(i2xi(idx,q,cut,Ni[i]), i2h(idx,q,cut,Ni[i]), cut, rad)
        print (xi[i])/np.sqrt(2), h[i]/np.sqrt(2)

    # output data
    gas['type'] = "Rect"
    gas['v'] = vecToStr([0,0,0])
    gas['vs'] = b64encode(xi.flatten())
    gas['vvs'] = b64encode(h.flatten())

    # open resulting .kep file and dump the modified data 
    with open(sys.argv[2], 'w') as f:
        json.dump(data, f, indent=2, sort_keys=True)
 
