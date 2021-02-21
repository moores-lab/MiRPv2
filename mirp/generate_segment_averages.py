#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
MiRP - a microtubule RELION-based pipeline for cryo-EM image processing. 
This script is dependent on EMAN2 (tested with v2.13) and generates microtubule 
particles averaged over a 7-particle window
"""

__author__ = 'Alexander D. Cook & Joseph Atherton'
__license__ = 'GPLv3'
__version__ = '2.0'


import starfileIO
import microtubules
import helper_fns
import argparse
import sys
import EMAN2 
import re
import os

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--in_parts', required=True, help='Input starfile (e.g. Class3D/PF_sorting/run_it001_data.star')
parser.add_argument('-o', '--o', required=True, help='Output folder for segment averages (e.g. Extract/seg_averages')
args = parser.parse_args()

#can only make a new directory for segment averages (to avoid overwriting something else)
assert not os.path.exists(args.o), 'Error, output path already exists!'
os.makedirs('%s/Micrographs' % args.o)

print('Reading in starfile...')
mts = microtubules.Microtubules(args.in_parts, '.')
name_regex = re.compile('.+/(.+\.mrcs)')

for ix, microtubule in enumerate(mts):
    sys.stdout.write('\rGenerating segment averages for microtubule %i for %i' % (ix, mts.mt_tot))
    sys.stdout.flush()
    #read particles from RELION stack into EMAN2
    particle_stack = microtubule['rlnImageName'][0][7:]
    spart = int(microtubule['rlnImageName'][0][:6])
    epart = int(microtubule['rlnImageName'][-1][:6])
    part_range = [i for i in range(spart-1, epart)]
    images = EMAN2.EMData.read_images(particle_stack, part_range)
    microtubule['images'] = images
    mt_len = len(microtubule['rlnImageName'])
    #prepare output 
    mname = re.search(name_regex, microtubule['rlnImageName'][0]).group(1)
    mname = mname.replace('.mrcs', '_SAs.mrcs')
    outfile = args.o + '/Micrographs/' + mname 

    #2D transform partices prior to averaging
    for ptcl_index in range(mt_len):
        psi = microtubule['rlnAnglePsi'][ptcl_index]
        xsh = microtubule['rlnOriginXAngst'][ptcl_index] / mts.apix
        ysh = microtubule['rlnOriginYAngst'][ptcl_index] / mts.apix
        t = EMAN2.Transform()
        t.set_params({'type':'2d', 'alpha':psi, 'tx':xsh, 'ty':ysh})
        microtubule['images'][ptcl_index].transform(t)

    #average along a sliding window
    for ptcl_index in range(mt_len):
        low, hi = helper_fns.get_window(ptcl_index, 3, 4, mt_len)
        avg = EMAN2.Averagers.get('mean')
        for i in range(low, hi):
            avg.add_image(microtubule['images'][i])
        avg_img = avg.finish()
        #transform segment average back to original position and save
        psi = microtubule['rlnAnglePsi'][ptcl_index]
        xsh = microtubule['rlnOriginXAngst'][ptcl_index] / mts.apix
        ysh = microtubule['rlnOriginYAngst'][ptcl_index] / mts.apix
        t = EMAN2.Transform()
        t.set_params({'type':'2d', 'alpha':-psi, 'tx':-xsh, 'ty':-ysh})
        avg_img.transform(t)
        avg_img.append_image(outfile)

#update per micrograph information in starfile
starfile_data = starfileIO.Starfile(args.in_parts)
starfile_data.read_star()    
starfile_data.sort_loop_datablock('data_particles', 'rlnMicrographName', 'rlnHelicalTubeID', 'rlnHelicalTrackLengthAngst')
particles = starfile_data.get_datablock('data_particles')
micrographs = helper_fns.group_dict_of_list(particles, 'rlnMicrographName')
image_names = []
def format_particle_number(s1, s2):
    return (s1[:-len(s2)] + s2)

for mgph in micrographs:
    for index, name in enumerate(mgph['rlnImageName']):
        mgphname = re.search(name_regex, name).group(1)
        mgphname = mgphname.replace('.mrcs', '_SAs.mrcs')
        pnum = format_particle_number('000000', str(index+1))
        image_name = pnum + '@' + args.o + '/Micrographs/' + mgphname
        image_names.append(image_name)
        
particles['rlnImageName'] = image_names
starfile_data.add_datablock('data_particles', particles)
starfile_data.write_star('%s/segment_averages.star' % args.o)
print('\nFinished! Wrote %s/segment_averages.star' % args.o)    