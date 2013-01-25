import numpy, re, math, itertools

from array import array

import json
from base64 import b64decode

from scipy import integrate
from scipy.ndimage import map_coordinates

from kepy.element import element

from kepy.io import read_nodes, read_elements, readNodesElems

def readNodesCells(filename):
    nodes, cells, facets = readNodesElems(filename)
    return nodes, cells

def isNumber(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def readMacros(filename, size):
    with open(filename, "rb") as f:
        numbers = (filter(isNumber, line.split()) for line in f.readlines())
        data = [ (int(line[0]), map(float, line[1:])) for line in numbers ]
        _, data = zip(*sorted(data, key = lambda pair: pair[0]))
        data = zip(*data)
    return data

def splitFacetsToVerges(facets):
    vdic = {}
    for facet in facets:
        n = lambda i: facet.vertexes[i]
        verges = [ sortPair( (n(i), n(j)) ) for i, j in verges_dic[facet.type] ]
        for v1, v2 in listToPairs(verges):        
            if v1 in vdic:
                if not v2 in vdic[v1]:
                    vdic[v1].append(v2)
            else:
                vdic[v1] = [v2]
    return vdic

def lineCrossPlane(a, b, u, v):
    g = gamma(u, v, b-a)
    if g != 0:
        a1 = - gamma(v, a, b-a) / g
        a2 =   gamma(u, a, b-a) / g
        a3 = - gamma(a, u, v)   / g
        return a1, a2, a3

def pointOnPlane(p, u, v):
    return gamma(p, u, v) == 0

def lineCrossFacet(O, u, p1, p2, p3):
    o = O - p1
    n1 = p2 - p1
    n2 = p3 - p1
    result = lineCrossPlane(o, o+u, n1, n2)
    if result:
        a1, a2, a3 = result
        if (a1 >= 0) & (a2 >= 0) & ((a1+a2) <= 1):
            return a3

def sortV(v):
    return v if v[0] <= v[1] else (v[1], v[0])

def cellsToVerges(cells):
    verges = dict()
    for l, cell in enumerate(cells):
        vs = verges_dic[ cell.type ]
        vs = [ sortV( (cell.vertexes[i], cell.vertexes[j]) ) for i, j in vs ]
        cell.set_verges(vs)
        for v in vs:
            if v in verges:
                verges[v].append(l)
            else:
                verges[v] = [l]
    return verges

def intersectEdges(verges, cells, nodes, O, u, v):
    points = []
    l = 0
    for (i, j), cs in verges.iteritems():
        if gamma(nodes[i] - O, u, v) * gamma(nodes[j] - O, u, v) <= 0:
            p = vergeCrossPlane(nodes[i], nodes[j], O, u, v)
            points.append(numpy.array(p))
            for cell in cs:
                cells[cell].add_point(l)
            l += 1
    return points

def cellCrossPlane(cell, nodes, O, u, v):
    verges = verges_dic[ cell.type ]
    for i, j in verges:
        p = vergeCrossPlane(nodes[cell.vertexes[i]], nodes[cell.vertexes[j]], O, u, v)
        if p:
            yield p

def intersect(O, u, cell, nodes):
    line = []
    for facet in facets_dic[cell.type]:
        a3 = lineCrossFacet(O, u, nodes[cell.nodes[facet[0]]], 
                                  nodes[cell.nodes[facet[1]]], 
                                  nodes[cell.nodes[facet[-1]]])
        if a3:
            line.append(a3)     
    if line:
        return center(line) 

def vergeCrossPlane(p1, p2, O, u, v):
    result = lineCrossPlane(p1 - O, p2 - O, u, v)
    if result:
        a1, a2, a3 = result
        if (a3 > 0) & (a3 < 1):
            return a1, a2

def vergesCrossPlane(vdic, nodes, O, u, v):
    #clean keys
    for verge in vdic.keys():
        v1, v2 = verge
        ps = vergeCrossPlane(nodes[v1], nodes[v2], O, u, v)
        if ps:
            vdic[verge] = (numpy.array(ps), vdic[verge])
        else:
            del vdic[verge]
    #clean values
    for key, value in vdic.items():
        p, neigbors = value
        newneigbors = [verge for verge in neigbors if verge in vdic]
        vdic[key] = (p, newneigbors)

    return vdic

def vergesToLines(vdic):
    while vdic.keys():
        first = fst(vdic.keys()) 
        pred = first
        line = [fst(vdic[pred])]
        curr = fst(snd(vdic[pred]))
        del vdic[pred]
        while curr != first:
            line.append(fst(vdic[curr]))
            v1, v2 = snd(vdic[curr])
            del vdic[curr]
            if v1 != pred:
                nxt = v1
            else: 
                nxt = v2
            pred = curr
            curr = nxt
        yield line

def inCell(cell, nodes, p):
    for facet in facets_dic[cell.type]:
        g = gamma( nodes[cell.vertexes[facet[0]]]  - p, 
                   nodes[cell.vertexes[facet[1]]]  - p, 
                   nodes[cell.vertexes[facet[-1]]] - p )
        if g > 0:
            return False
    return True

def angle(x, y):
    a = math.acos( numpy.dot(x, y) / 
            math.sqrt( numpy.dot(x, x) * numpy.dot(y, y) ) )
    if (numpy.cross(x, y) >= 0):
        return a
    else:
        return -a

def triangles(points, nodes):
    ps = [nodes[p] for p in points]
    c = center(ps)
    pair = zip(ps, points)
    pair.sort( key = lambda x: angle( x[0]-c, numpy.array( [1., 0.] ) ) )
    sorted_points = zip(*pair)[1]
    trs = []
    for i, j in zip(sorted_points[1:], sorted_points[2:]):
        trs.append( [sorted_points[0], i, j] )
    return trs

def rot2(x, y):
    return x[0] * y[1] - x[1] * y[0]

def inTr(x, tr, points):
    v0 = points[tr[0]] - x
    v1 = points[tr[1]] - x
    v2 = points[tr[2]] - x
    g1 = rot2(v0, v1) 
    g2 = rot2(v1, v2) 
    g3 = rot2(v2, v0) 
    return (g1 * g2 >= 0 ) and (g1 * g3 >= 0)
        
def l2Min(cell, nodes):
    return min( [ sqr(nodes[i]-nodes[j]) for i,j in listToPairs(cell.vertexes) ] )
    
def lMin(cells, nodes):
    l2_min = min( [l2Min(cell, nodes) for cell in cells] )
    return math.sqrt(l2_min)
 
def rect(points):
    xs, ys = zip(*points)
    return (min(xs), min(ys)), (max(xs), max(ys))

def rectI(p1, p2, step):
    x1, y1 = p1
    x2, y2 = p2
    return (toL(x1, step)-1, toL(y1, step)-1), (toL(x2, step)+1, toL(y2, step)+1)

def addCellToCover(cover, si, sj, step, cell, cell_i, nodes, O, u, v):
    ps = [p for p in cellCrossPlane(cell, nodes, O, u, v)]
    if ps:
        p1, p2 = rect(ps)
        i1, i2 = rectI(p1, p2, step)
#       print ps
#       print p1, p2
#       print i1, i2
        for i, j in between(i1, i2):
#           print i, j, cell_i
            p = O + toI(i, step) * u + toI(j, step) * v
            if inCell(cell, nodes, p):
                cover[i-si][j-sj].append(cell_i)

def projToPlane(c, u, v):
    return numpy.array( ( numpy.dot(c, u) / numpy.dot(u, u), 
                          numpy.dot(c, v) / numpy.dot(v, v) ) )
   
def projPointToLine(p, p1, p2):
    s1 = p1 - p
    s2 = p2 - p1
    return - numpy.dot(s1, s2) / numpy.dot(s2, s2)

def lineToIs(p1, p2, step, cover, si, sj):
    q1, q2 = rect([p1, p2])
    print q1, q2
    p1, p2 = numpy.array(p1), numpy.array(p2)
    i1, i2 = rectI(q1, q2, step)
    reserve, rs = [], []
    for i, j in between(i1, i2):
        if not cover[i-si][j-sj]:
            p = numpy.array((toI(i, step), toI(j, step)), dtype=numpy.float64)
            t = projPointToLine(p, p1, p2)
            if 0 < t < 1:
                o = p1 + t * (p2-p1)
                q = o - p
                if numpy.dot(q, q) < 2:
                    yield i, j, o
            else:
                r = numpy.dot(p-p1, p-p1)
                if r < 1:
                    rs.append(r)
                    reserve.append((i, j))
#    if rs:
#        r, (i, j) = min(zip(rs, reserve), key = fst)
#        yield i, j, p1



def addLineToCover(p1, p2, step, cover, X, Y, si, sj, cell_i):
    for i_, j_, p in lineToIs(p1, p2, step, cover, si, sj):
        i, j = i_ - si, j_ - sj
        cover[i][j] = [cell_i]
        X[i][j], Y[i][j] = p

tetrahedron_facets = [ [0, 1, 2], [0, 2, 3], [0, 3, 1], [1, 3, 2] ] 
prism_facets       = [ [0, 1, 2], [3, 5, 4],
                       [0, 2, 5], [5, 3, 0],
                       [2, 1, 4], [4, 5, 2],
                       [0, 3, 4], [4, 1, 0] ]
hexahedron_facets  = [ [0, 1, 2], [2, 3, 0],
                       [0, 4, 5], [5, 1, 0],
                       [1, 5, 6], [6, 2, 1],
                       [2, 6, 7], [7, 3, 2],
                       [0, 3, 7], [7, 4, 0],
                       [4, 7, 6], [6, 5, 4] ]
facets_dic = {4: tetrahedron_facets, 5: hexahedron_facets, 6: prism_facets}

tetrahedron_verges = [ (0, 1), (1, 2), (2, 0), (3, 0), (3, 1), (3, 2) ]
prism_verges       = [ (0, 1), (1, 2), (2, 0), 
                       (3, 4), (4, 5), (5, 3), 
                       (0, 3), (1, 4), (2, 5) ]
hexahedron_verges  = [ (0, 1), (1, 2), (2, 3), (3, 0),
                       (4, 5), (5, 6), (6, 7), (7, 4),
                       (0, 4), (1, 5), (2, 6), (3, 7) ]
tringle_verges     = [ (0, 1), (1, 2), (2, 0) ]
quadrangle_verges  = [ (0, 1), (1, 2), (2, 3), (3, 0) ]
verges_dic = { 4: tetrahedron_verges, 5: hexahedron_verges, 6: prism_verges,
               2: tringle_verges, 3: quadrangle_verges }

def gamma(c, u, v):
    a = numpy.array( [ c, u, v ] )
    return numpy.linalg.det(a)
 
def center(ps):
    return sum(ps) / len(ps)

def listToStr(elm):
    return str(len(elm)) + ' ' + ' '.join( map(str, elm) )

def list3ToStr(l):
    return ' '.join( map(str, l) )

def sortPair(pair):
    x, y = pair
    if y < x:
        return y, x
    else:
        return x, y

def listToPairs(l):
    size = len(l)
    for i in range(size):
        for j in range(size):
            if i != j:
                yield (l[i], l[j])
def between(p1, p2):
    i1, j1 = p1
    i2, j2 = p2
    for i in range(i1, i2+1):
        for j in range(j1, j2+1):
            yield (i, j)

def pairs(l):
    return zip(l, l[1:]+[l[0]])

def sqr(x):
    return numpy.dot(x, x)

def fst(l):
    return l[0]
def snd(l):
    return l[1]

def toI(x, step):
    return (x+0.5)*step

def toL(y, step):
    x = y / step
    if x >= 0:
        return int(x)
    else:
        return int(x)-1
def toG(y, step):
    x = y / step
    if x >= 0:
        return int(x)+1
    else:
        return int(x)

def pointsCross(points):
    a = numpy.array( [  points[1] - points[0],
                        points[2] - points[0],
                        points[3] - points[0] ] );  
    return numpy.linalg.det(a) / 6

def tetrahedronVolume(elm, points):
    ps = [points[i] for i in elm.vertexes]
    return pointsCross(ps) / 6.    

def prismVolume(elm, points):
    ps = [points[i] for i in elm.vertexes]
    return pointsCross([ps[0], ps[1], ps[2], ps[3]]) + \
           pointsCross([ps[1], ps[3], ps[4], ps[5]]) + \
           pointsCross([ps[1], ps[2], ps[3], ps[5]])  

def hexahedronVolume(elm, points):
    ps = [points[i] for i in elm.vertexes]
    return pointsCross([ps[0], ps[1], ps[3], ps[4]]) + \
           pointsCross([ps[1], ps[2], ps[3], ps[6]]) + \
           pointsCross([ps[1], ps[4], ps[5], ps[6]]) + \
           pointsCross([ps[3], ps[4], ps[6], ps[7]]) + \
           pointsCross([ps[1], ps[3], ps[4], ps[6]])

volumes_dic = {4: tetrahedronVolume, 5: hexahedronVolume, 6: prismVolume}

def cellVolume(cell, nodes):
    return volumes_dic[cell.type](cell, nodes)


def circle(rad):
    i = 0
    for i1 in range(-2*rad+1, 2*rad, 2):
        for i2 in range(1, 2*rad, 2):
            s = i1 * i1 + i2 * i2
            if s < 4*rad*rad:
                i += 1
    return i

sqrt2 = math.sqrt(2.)

def makeG():
    n = 100
    shape = (n, n)
    gamma = numpy.fromfunction(lambda i, j: 0.5   * (i+0.5) / n, shape)
    delta = numpy.fromfunction(lambda i, j: sqrt2 * j / n, shape)
    cos   = (1 - delta * delta) / 2 / gamma
    cos[cos < -1] = -1
    cos[cos >  1] =  1
    phi   = numpy.arccos(cos)

    def pyfunc(gamma, phi0):
        return integrate.quad(lambda x: math.sqrt(1-2*gamma*math.cos(x)), phi0, math.pi)[0]
#        return math.pi - phi0
    func = numpy.frompyfunc(pyfunc, 2, 1)
    G = func(gamma, phi)

    return gamma, delta, numpy.array(G, dtype='float64')

gamma_glob, delta_glob, G_glob = makeG()

def my_interp(x, y, z, xs, ys):

    nx, ny = z.shape
    xmin, xmax = numpy.min(x), numpy.max(x)
    ymin, ymax = numpy.min(y), numpy.max(y)

    xi = numpy.array(xs, dtype=numpy.float)
    yi = numpy.array(ys, dtype=numpy.float)

    xi = (nx - 1) * (xi - xmin) / (xmax - xmin)
    yi = (ny - 1) * (yi - ymin) / (ymax - ymin)

    zs = map_coordinates(z, [xi, yi], order=0)

    return zs


def readF(filename, length, gr, ms, rads, dvs):

    def make_xr(rad, dv):
        dim = (2*rad, rad) 
        x = numpy.fromfunction(lambda i, j: i - rad + 0.5, dim)
        r = numpy.fromfunction(lambda i, j: j       + 0.5, dim)
        x *= dv
        r *= dv

        return x, r

    def make_xrf(fd, rad, x, r, dv):
        a = array.array('d')
        size = numpy.sum(x*x + r*r < rad*rad*dv*dv)
        a.fromfile(fd, size)

        f = numpy.zeros_like(x)
        f[x*x + r*r < rad*rad*dv*dv] = numpy.array(a)

        n  = numpy.sum(f / r)
        vx = numpy.sum(f / r * x) / n 
        tx = numpy.sum(f / r * (x - vx)**2) / n
        tr = numpy.sum(f / r * r**2) / n

        fm = r * numpy.exp(- (x - vx)**2 / 2 / tx - r**2 / 2 / tr)
        fm *= n / numpy.sum(fm / r)

        return f

    xs, rs = zip(*[make_xr(rad, dv) for rad, dv in zip(rads, dvs)])

    def combs(l):
        for i, x in enumerate(l):
                for j in range(i, len(l)):
                        yield x, l[j]
    xpairs = combs(xs)
    rpairs = combs(rs)
    mpairs = combs(ms)

    ggs = [calcGG(x1, x2, r1, r2, gr / math.sqrt(2 * m1 * m2 / (m1 + m2))) for (x1, x2), (r1, r2), (m1, m2) in zip(xpairs, rpairs, mpairs)]

    l = []
    qr = []

    with open(filename, "rb") as fd:
        for i in range(length):
            data = fd.read(4+4)
            i, size = struct.unpack('=ii', data)
            print i, size, rads

            fs = [make_xrf(fd, rad, x, r, dv) for rad, x, r, dv in zip(rads, xs, rs, dvs)]

            fpairs = combs(fs)

            ffs = [numpy.outer(*fpair).flatten() for fpair in fpairs]

            l.append  ( i )
            qr.append ( [numpy.sum(ff * gg) for gg, ff in zip(ggs, ffs)] )

        _, qr = zip(*sorted(zip(l, qr), key = lambda pair: pair[0]))

    return zip(*qr)

def isXmlStart(key, line):
    return re.search(r"^<" + key + r"\b", line)
def isXmlEnd(key, line):
    return re.search(r"</" + key + r">$", line)
def getXmlValue(s, key):
    print s, key
    return re.search(key + r"\s*=\s*\"([^\"]+)\"", s).groups()[0]

def readRadsMassesDvs(filename):
    rads = []
    masses = []
    dvs = []
    with open(filename, "rb") as fd:
        line = " "
        while line:
            line = fd.readline()
            if isXmlStart("gas", line):
                cut = float(getXmlValue(line, "cut"))
                rad = float(getXmlValue(line, "rad"))
                try:
                    str = getXmlValue(line, "masses")
                    masses = map(float, str.split())
                except AttributeError:
                    masses = [1.]
                try:
                    str = getXmlValue(line, "adds")
                    adds = map(float, str.split())
                except AttributeError:
                    adds = [0. for m in masses]
                print "adds = ", adds
                try:
                    v0 = float(getXmlValue(line, "v"))
                except AttributeError:
                    v0 = 0
                rads = [rad for m in masses]
                dvs = [(cut / math.sqrt(m) + a) / r for r, m, a in zip(rads, masses, adds)] 
                break
    return rads, masses, dvs, v0

def makeXR(rad):
    dim = (2*rad, rad) 
    x = numpy.fromfunction(lambda i, j: i - rad + 0.5, dim)
    r = numpy.fromfunction(lambda i, j: j       + 0.5, dim)
    return x, r

def readFunc(fd, x, r, rad):
    f = numpy.zeros_like(x)

    a = array.array('d')
    size = numpy.sum(x*x + r*r < rad*rad)
    print "size = ", size
    a.fromfile(fd, size)

    f[x*x + r*r < rad*rad] = numpy.array(a)

    return f

def calcGG(x1, x2, r1, r2, g):
    ones1 = numpy.ones(x1.shape)
    ones2 = numpy.ones(x2.shape)
    xisqr = numpy.outer(x1*x1 + r1*r1, ones2) + numpy.outer(ones1, x2*x2 + r2*r2) - \
          2*numpy.outer(x1, x2)
    xirr  = numpy.outer(r1, r2)
    gammas = (xirr / xisqr).flatten()
    deltas = (g / numpy.sqrt(xisqr)).flatten()
    deltas[deltas > sqrt2] = sqrt2

    gg = my_interp(gamma_glob, delta_glob, G_glob, gammas, deltas) * numpy.sqrt(xisqr.flatten())
    return gg
 
def calcQr(x1, x2, r1, r2, f1, f2, g):
    ff = numpy.outer(f1, f2).flatten()
    return numpy.sum(ff * calcGG(x1, x2, r1, r2, g))

