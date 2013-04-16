import math
import numpy as np
from base64 import b64decode
from array import array

def sqr(x):
    return x * x

def cube(x):
    return x * x * x

def read_ximesh(data):
    cut = data["gas"]["cut"]
    rad = data["gas"]["rad"]
    symmetry = data["gas"]["symmetry"]
    type = data["gas"]["type"]

    print "rad = ", rad

    if type == "Simple":

        if symmetry == "Cylindrical":
            v = float(data["gas"].get("v", "0."))
            dim = (2*rad, rad) 
            x = np.fromfunction(lambda i, j: v + (i + 0.5 - rad) * cut / rad,
                                dim, dtype=float)
            y = np.fromfunction(lambda i, j: (j + 0.5) * cut / rad,
                                dim, dtype=float)
            vol = y
            r = sqr(x-v) + sqr(y)
            ax = (x, y)
            d3v = 2 * math.pi * sqr(cut / rad)
            
        elif symmetry == "Cartesian":
            v = map(float, data["gas"].get("v", "( 0. 0. 0. )")[1:-1].split() )
            vx, vy, vz = v
            dim = (2*rad, 2*rad, 2*rad)
            x = np.fromfunction(lambda i, j, k: vx + (i + 0.5 - rad) * cut / rad,
                                dim, dtype=float)
            y = np.fromfunction(lambda i, j, k: vy + (j + 0.5 - rad) * cut / rad,
                                dim, dtype=float)
            z = np.fromfunction(lambda i, j, k: vz + (k + 0.5 - rad) * cut / rad,
                                dim, dtype=float)
            vol = np.ones_like(y)
            r = sqr(x-vx) + sqr(y-vy) + sqr(z-vz)
            ax = (x, y, z)
            d3v = cube(cut / rad)

        circle = r < cut*cut

    elif type == "Rect":

        print "type = Rect"

        if "hcenter" in data["gas"]:

            hcenter = data["gas"]["hcenter"]

            def find_q(d, n, err=1e-10):
                # first approx
                q = 2 * d / n - 1
                while True:
                    qq = math.pow(q, n-1)
                    f  = (q * qq - 1) / (q - 1)
                    if (abs(f - d) < err):
                        break
                    df = ((n-1) * q * qq - n * qq + 1) / sqr(q-1)
                    dq = (d - f) / df
                    q += dq
                return q

            q = find_q(cut / hcenter, rad)
            print "q = ", q

            vs, vols = [], []
            h = hcenter
            x = h / 2
            for i in range(rad):
                vs.append(x)
                vols.append(h)
                x += h / 2
                h *= q
                x += h / 2

            vs2   = np.array([-v for v in vs[::-1]]   + vs)
            vols2 = np.array([ v for v in vols[::-1]] + vols)
            vs    = np.array(vs)
            vols  = np.array(vols)

            if symmetry == "Cylindrical":
                dim = (2*rad, rad) 
                x = np.array( [ [ vs2[i]
                                  for j in range(rad) ] for i in range(2*rad) ] )
                y = np.array( [ [ vs[j]
                                  for j in range(rad) ] for i in range(2*rad) ] )
                vol = np.array( [ [ vols2[i]*vols[j]
                                  for j in range(rad) ] for i in range(2*rad) ] )
                vol *= y
                r = x*x + y*y 
                ax = (x, y)
                d3v = 2 * math.pi

            elif symmetry == "Cartesian":
                dim = (2*rad, 2*rad, 2*rad)
                x = np.array( [ [ [ vs2[i]
                        for k in range(2*rad) ] for j in range(2*rad) ] for i in range(2*rad) ] )
                y = np.array( [ [ [ vs2[j]
                        for k in range(2*rad) ] for j in range(2*rad) ] for i in range(2*rad) ] )
                z = np.array( [ [ [ vs2[k]
                        for k in range(2*rad) ] for j in range(2*rad) ] for i in range(2*rad) ] )
                vol = np.array( [ [ [ vols2[i]*vols2[j]*vols2[k]
                        for k in range(2*rad) ] for j in range(2*rad) ] for i in range(2*rad) ] ) 
                r = x*x + y*y + z*z
                ax = (x, y, z)
                d3v = 1. 

            circle = r < cut*cut

        else:
            
            vs_str  = data["gas"]["vs"]
            vvs_str = data["gas"]["vvs"]

            vs  = array('d', b64decode(vs_str))
            vvs = array('d', b64decode(vvs_str))

            if symmetry == "Cylindrical":
                dim = (2*rad, rad) 
                x = np.array( [ [ vs[i]
                                  for j in range(rad) ] for i in range(2*rad) ] )
                y = np.array( [ [ vs[j + 2*rad]
                                  for j in range(rad) ] for i in range(2*rad) ] )
                vol = np.array( [ [ vvs[i]*vvs[j + 2*rad]
                                  for j in range(rad) ] for i in range(2*rad) ] )
                vol *= y
                r = x*x + y*y 
                ax = (x, y)
                d3v = 2 * math.pi

                cut_x, cut_y = 0, 0
                for i in range(2*rad):
                    xi = abs(vs[i]) + 0.5 * vvs[i]
                    if cut_x < xi:
                        cut_x = xi
                for j in range(rad):
                    xi = abs(vs[j + 2*rad]) + 0.5 * vvs[j + 2*rad]
                    if cut_y < xi:
                        cut_y = xi
                circle = ((x / cut_x)**2 + (y / cut_y)**2) < 1 

            elif symmetry == "Cartesian":
                dim = (2*rad, 2*rad, 2*rad)
                x = np.array( [ [ [ vs[i]
                        for k in range(2*rad) ] 
                            for j in range(2*rad) ] 
                                for i in range(2*rad) ] )
                y = np.array( [ [ [ vs[j + 2*rad]
                        for k in range(2*rad) ] 
                            for j in range(2*rad) ] 
                                for i in range(2*rad) ] )
                z = np.array( [ [ [ vs[k + 4*rad]
                        for k in range(2*rad) ] 
                            for j in range(2*rad) ] 
                                for i in range(2*rad) ] )
                vol = np.array( [ [ [ vvs[i] * vvs[j + 2*rad] * vvs[k + 4*rad]
                        for k in range(2*rad) ] 
                            for j in range(2*rad) ] 
                                for i in range(2*rad) ] ) 
                r = x*x + y*y + z*z
                ax = (x, y, z)
                d3v = 1.0

                cut_x, cut_y, cut_z = 0.0, 0.0, 0.0
                for i in range(2*rad):
                    xi = abs(vs[i]) + 0.5 * vvs[i]
                    if cut_x < xi:
                        cut_x = xi
                for j in range(2*rad):
                    xi = abs(vs[j + 2*rad]) + 0.5 * vvs[j + 2*rad]
                    if cut_y < xi:
                        cut_y = xi
                for k in range(2*rad):
                    xi = abs(vs[k + 4*rad]) + 0.5 * vvs[k + 4*rad]
                    if cut_y < xi:
                        cut_y = xi
                circle = ((x / cut_x)**2 + (y / cut_y)**2 + (z / cut_z)**2) < 1  

    return symmetry, rad, circle, ax, vol, r, d3v

def make_e(symm, ax, u=None):
    if symm == "Cylindrical":
        x, y = ax
        if u is None:
            u = 0.
        return sqr(x-u) + sqr(y)
    elif symm == "Cartesian":
        x, y, z = ax
        if u is None:
            ux, uy, uz = 0., 0., 0.
        else:
            ux, uy, uz = u
        return sqr(x-ux) + sqr(y-uy) + sqr(z-uz)
