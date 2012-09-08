#!/usr/bin/env python

import sys, math, numpy
import pylab
from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, LinearLocator, MaxNLocator

from out2 import *

nodes, cells = readNodesCells(sys.argv[1])

centers = [ center( [ nodes[i] for i in cell.vertexes ] ) for cell in cells ]

#xs = [-2, -1, 0, 1, 2, 5, 9, 14, 25]
#globtitle = r"$\mathrm{M}=5$, $m^\alpha/m^\beta = 1/33$, $n^\alpha/n^\beta = 500/1$, "

#xs = [-10, -5, -2, 2, 5]
#globtitle = r"$\mathrm{M}=5$, $m^\alpha/m^\beta = 1/33$, $n^\alpha=n^\beta$, "

xs = [-5, -4, -3, -1, 1, 5, 9, 15]
globtitle = r"$\mathrm{M}=5$, $m^\alpha/m^\beta = 1/33$, $n^\alpha/n^\beta = 50/1$, "


ps = [ numpy.array( [x, 0, 0] ) for x in xs ]
print ps

def numpy_sqr(x):
    return numpy.dot(x, x)

def i_min_dist(x, xs):
    j = 0
    d = numpy_sqr(xs[0] - x)
    for i, c in enumerate(xs[1:]):
        d1 = numpy_sqr(c - x)
        if d1 < d:
            d = d1
            j = i

    return j

ics = [ i_min_dist(x, centers) for x in ps ]

for ic in ics:
    print "ic, center[ic] = ", ic, centers[ic]


colors = ['b', 'r', 'g', 'k']

fd = file(sys.argv[2])

rads, ms, dvs, v0 = readRadsMassesDvs(sys.argv[1])
print rads, ms, dvs
xrs = [makeXR(rad) for rad in rads]

params = {
            'backend'         : 'ps',
            'axes.labelsize'  : 8,
            'text.fontsize'   : 8,
            'legend.fontsize' : 9,
            'xtick.labelsize' : 8,
            'ytick.labelsize' : 8,
            'text.usetex'     : False,
            'font.family'     : 'serif'
        }
pylab.rcParams.update(params)

def get_fig_size(columnwidth, height):
    fig_width_pt  = columnwidth
    inches_per_pt = 1.0 / 72.27
    fig_width     = fig_width_pt * inches_per_pt
    fig_height    = height * inches_per_pt
    fig_size      = [fig_width, fig_height]
    return fig_size

fig_fx = pylab.figure(1, figsize = get_fig_size(400, 400))
titles = [r'$f^\alpha$', r'$f^\beta$']
cs = len(titles)
axs_fx = [ fig_fx.add_subplot(cs, 1, 1+i) for i in range(cs) ]
fig_fx.subplots_adjust(left=0.1, right=0.95, bottom=0.05, top=0.95)

figs = [ pylab.figure(i+2, figsize = get_fig_size(400, 1.618*400)) for i, _ in enumerate(xs) ]
axs  = [ fig.add_subplot(1+cs, 1, 1, projection='3d', title=globtitle + "$x="+str(x)+"$") for fig, x in zip(figs, xs) ]
axss = [ [ fig.add_subplot(1+cs, 1, 2+i, projection='3d') for i in range(cs) ] for fig in figs ]

figslog = [ [ pylab.figure(2+len(xs)+i*cs+j, figsize = get_fig_size(200, 200)) for j in range(cs)] for i, _ in enumerate(xs) ]
axslog  = [ [ x2.add_subplot(111, title=title + ", $x="+str(x)+"$") for x2, title in zip(x1, titles)] for x1, x in zip(figslog, xs) ]

for fig in figs:
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)

while fd:
    data = fd.read(4+4)
    try:
        i, size = struct.unpack('=ii', data)
    except struct.error:
        break
    print i, size

    fs   = [readFunc(fd, x, r, rad) for (x, r), rad in zip(xrs, rads)]

    for j, ic in enumerate(ics):
        if (i == ic):
            xrs1 = [(x * dv + v0, r * dv) for (x, r), dv in zip(xrs, dvs)]

            fs  = [ f / dv**3 for f, dv in zip(fs, dvs)]

            for k, ( (x, r), f, c, ax_fx, title, fig ) in enumerate(zip(xrs1, fs, colors, axs_fx, titles, figs)):

                """
                gs  = numpy.linspace(0, 0.7 * numpy.max(r), 20)
                qrs = [calcQr(x, r, f, x, r, f, g) for g in gs]
                qrs = [q1 - q2 for q1, q2 in zip(qrs, qrs[1:]+[0])]
                ax1.plot(gs, qrs, color = c)
                """

                f = f / r

                ax_fx.plot(x[:,0], f[:,0], label = "$x="+str(xs[j])+"$")

                f = numpy.ma.masked_where(f == 0.0, f)

                paxs = [ axs[j], axss[j][k] ]

                paxs[0].set_zlabel(','.join(titles))
                paxs[1].set_zlabel(title)

                for pax in paxs:

                    pax.set_xlabel(r'$\xi_y$')
                    pax.set_ylabel(r'$\xi_x$')
                    pax.view_init(35, -15)

                    for axis in pax.w_xaxis, pax.w_yaxis, pax.w_zaxis:
                        axis.pane.set_visible(False)
                        axis.gridlines.set_visible(False)
                        axis.set_rotate_label(False)

                    for elt in pax.w_zaxis.get_ticklines() + pax.w_zaxis.get_ticklabels():
                        elt.set_visible(False)
                    pax.w_zaxis.set_visible(False)

                    pax.plot_wireframe(r, x, f, color = c, rstride=1, cstride=1)

                pax = axslog[j][k]

                f = numpy.concatenate(( f[:,::-1], f), axis=1)
                r = numpy.concatenate((-r[:,::-1], r), axis=1)
                x = numpy.concatenate(( x[:,::-1], x), axis=1)
                f = -numpy.log(f)

                pax.contour(x, r, f, 15, colors='k')

pylab.figure(1)

for ax_fx, title in zip(axs_fx, titles):
    pylab.text(.95, 0.1, r'$\xi_x/v_0$',
            horizontalalignment = 'center',
            verticalalignment   = 'center',
            size = 'small', transform = ax_fx.transAxes)

    pylab.text(.5, .9, title,
            horizontalalignment = 'center',
            verticalalignment   = 'center',
            size = 'small', transform = ax_fx.transAxes)

    ax_fx.spines['left'].set_position(('axes',0.0))
    ax_fx.spines['right'].set_color('none')
    ax_fx.spines['bottom'].set_position(('data',0.0))
    ax_fx.spines['top'].set_color('none')
    ax_fx.spines['left'].set_smart_bounds(True)
    ax_fx.spines['bottom'].set_smart_bounds(True)
    ax_fx.xaxis.set_ticks_position('bottom')
    ax_fx.yaxis.set_ticks_position('left')

    xmajorLocator = MaxNLocator(5)
    ax_fx.xaxis.set_major_locator(xmajorLocator)

    ymin, ymax = ax_fx.get_ylim()
    ax_fx.autoscale(tight = True)
    ax_fx.set_ylim( (ymin, ymax) )

    ax_fx.legend(loc = 'upper left', ncol = 2, frameon = False)

fig_fx.savefig("fx.eps")


for fig, x in zip(figs, xs):
    fig.savefig("f" + str(x) + ".eps", bbox_inches='tight')

filenames = ['fa', 'fb']
for x1, ax1, x in zip(figslog, axslog, xs):
    for x2, ax2, title, filename in zip(x1, ax1, titles, filenames):
        x2.savefig(filename + str(x) + '.eps')
        

