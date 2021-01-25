#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
MiRP - a microtubule RELION-based pipeline for cryo-EM image processing. 
This script is dependent on EMAN2 (tested with v2.13) and generates seam
references for 3D seam classification from a single 3D reference
"""

__author__ = 'Alexander D. Cook & Joseph Atherton'
__license__ = 'GPLv3'
__version__ = '2.0'


import EMAN2
import argparse
import microtubules
import starfileIO
import os

assert os.path.exists('default_pipeline.star'), 'default_pipeline.star not found. Please execute in a RELION directory'

parser = argparse.ArgumentParser()
parser.add_argument('-i', help='reference mrc file to be transformed')
parser.add_argument('-pf', help='protofilament number', type=int)
parser.add_argument('-r', help='helical rise (angstrom)', type=float)
parser.add_argument('-p', help='pixel size', type=float)
parser.add_argument('-o', help='output location. Will make folder if does not currently exist')
args = parser.parse_args()

assert not os.path.exists(args.o), 'Error, output path already exists!'
os.mkdir(args.o)

pftotal = args.pf
twist = -360.0 / float(pftotal)
rise = args.r
apix = args.p

outfiles = []
shifted_outfiles = []

print('Generating seam references...')
#for each possible relative seam position, create a new 3D volume 
for pf in range(1, pftotal+1):
    #read in 3D volume
    ref = EMAN2.EMData()
    ref.read_image(args.i)
    #get +/- symmetry operators
    seampos = microtubules.convert_pfnum_to_semicircle(pf, pftotal)
    #calculate change in angle required for that seam position
    Drot = seampos * twist
    #calculate change in z-axis translation for that seam position
    Dz = seampos * -rise / apix
    #use EMAN2 to transform the original 3D volume and save
    t = EMAN2.Transform({'type':'eman', 'az':Drot, 'tz':Dz})
    ref.transform(t)
    outfile = '%s/seamref%i.mrc' % (args.o, seampos)
    ref.write_image(outfile)
    outfiles.append(outfile)
    # transform and save 40 angstrom shifted volume for testing alpha/beta-tubulin register
    t=EMAN2.Transform({'type':'eman', 'tz':40/apix})
    ref.transform(t)
    outfile = '%s/seamref%i_40shift.mrc' % (args.o, seampos)    
    ref.write_image(outfile)
    shifted_outfiles.append(outfile)

#create a starfile for input into RELION
starout = starfileIO.Starfile('%s/seam_references.star' % args.o)
starout.add_datablock('data_references', {'rlnReferenceImage': outfiles})
starout.add_loop_data('data_references', {'rlnReferenceImage': shifted_outfiles})
starout.write_star('%s/seam_references.star' % args.o)
print('Done!')