#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MiRP - a microtubule RELION-based pipeline for cryo-EM image processing. 
MiRP v2 is designed to function with RELION v3.1, and is not compatible with earlier versions of RELION.
"""

__author__ = 'Alexander D. Cook & Joseph Atherton'
__license__ = 'GPLv3'
__version__ = '2.0'


import microtubules
import argparse
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

parser = argparse.ArgumentParser()
parser.add_argument('-i', required=True, help='Starfile to plot microtubule euler angles and XY shifts from.')
parser.add_argument('-o', required=False, help='Give file name, if saving a copy is desired.')
parser.add_argument('-n', required=False, type=int, help='The number of microtubule to plot.')
args = parser.parse_args()

mts = microtubules.Microtubules(args.i, '.')
num_mts = mts.mt_tot
if args.o:
    plot_pdf = PdfPages('%s.pdf' % args.o)        

if args.n:
    num_mts = args.n       

for mt in mts[:num_mts]:
    psi = mt['rlnAnglePsi']
    rot = mt['rlnAngleRot']
    tilt = mt['rlnAngleTilt']
    xsh = mt['rlnOriginXAngst']
    ysh = mt['rlnOriginYAngst']
    xax = [i for i in range(1, len(psi)+1)]
    yt = [i for i in range(-180, 180 + 1, 40) ]

    fig = plt.figure()
    pax = fig.add_subplot(1,4,1)
    pax.set_title('Psi')
    pax.plot(xax, psi, 'o')
    pax.set_ylabel('Angle')
    pax.set_xticks = xax
    pax.set_xlabel('Particle Number')
    pax.set_ylim([-181, 181])
    pax.set_yticks = yt

    tax = fig.add_subplot(1,4,2)
    tax.set_title('Tilt')
    tax.plot(xax, tilt, 'o')
    tax.set_xticks = xax
    tax.set_xlabel('Particle Number')
    tax.set_ylim([-181, 181])
    tax.set_yticks = yt

    rax = fig.add_subplot(1,4,3)
    rax.set_title('Rot')
    rax.plot(xax, rot, 'o')
    rax.set_xticks = xax
    rax.set_xlabel('Particle Number')
    rax.set_ylim([-181, 181])
    rax.set_yticks = yt

    sax = fig.add_subplot(1,4,4)
    sax.set_title('X/Y shifts')
    sax.plot(xax, xsh, 'o')
    sax.plot(xax, ysh, 'o')
    sax.set_ylabel('Angstrom')
    sax.set_xticks = xax
    sax.set_xlabel('Particle Number')

    plt.tight_layout()
    if args.o:
        plot_pdf.savefig(fig)
    else:
        plt.show()

if args.o:
    plot_pdf.close()
plt.close()

