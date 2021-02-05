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
import helper_fns
import argparse
import sys
import EMAN2 
import re
import os

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--in_parts', required=True, help='Input starfile (e.g. Class3D/PF_sorting/run_it001_data.star')
parser.add_argument('-o', '--o', required=True, help='Output folder for segment averages (e.g. Extract/seg_averages')
parser.add_argument('-r', '--revert', required=False, help='Reverts to non-averaged particles - give the path to the original extract folder (e.g. External/job011/Micrographs)')
args = parser.parse_args()


if args.revert:
    starfile_data = starfileIO.Starfile(args.in_parts)
    starfile_data.read_star()
    npath = args.revert.split('/')[1:]
    for ix, i in enumerate(starfile_data.get_entry('data_particles', 'rlnImageName')):
        sa = i.split('/')
        reverted = sa[0] 
        for d in npath:
            reverted = reverted + '/' + d
        reverted = reverted + '/' + sa[-1]
        starfile_data['data_particles']['rlnImageName'][ix] = reverted
    starfile_data.write_star('%s/particles_reverted_data.star' % args.o)
    print('Reverted to original particles from segement averages.\nparticles_reverted_data.star saved to %s' % args.o)

else:
    assert not os.path.exists(args.o), 'Error, output path already exists!'
    os.makedirs('%s/Micrographs' % args.o)

    def format_particle_number(s1, s2):
        return (s1[:-len(s2)] + s2)

    print('Reading in starfile...')
    starfile_data = starfileIO.Starfile(args.in_parts)
    starfile_data.read_star()
    apix = int(starfile_data.get_entry('data_optics', 'rlnImagePixelSize')[0])
    starfile_data.sort_loop_datablock('data_particles', 'rlnMicrographName', 'rlnHelicalTubeID', 'rlnHelicalTrackLengthAngst')
    particles = starfile_data.get_datablock('data_particles')
    micrographs = helper_fns.group_dict_of_list(particles, 'rlnMicrographName')
    nmgphs = len(micrographs)

    image_names=[]
    name_regex = re.compile('.+/(.+\.mrcs)')


    for ix, mgph in enumerate(micrographs):
        sys.stdout.write('\rGenerating microtubule segment averages for micrograph %i for %i' % (ix, nmgphs))
        sys.stdout.flush()
        images = EMAN2.EMData.read_images('%s' % mgph['rlnImageName'][0][7:])
        mgph['particles'] = images
        mts = helper_fns.group_dict_of_list(mgph, 'rlnHelicalTubeID')

        for mt in mts:
            mt_len = len(mt['rlnHelicalTubeID'])
            mname = re.search(name_regex, mt['rlnImageName'][0]).group(1)
            mname = mname.replace('.mrcs', '_SPs.mrcs')
            outfile = args.o + '/Micrographs/' + mname 

            #transform partices prior to averaging
            for ptcl_index in range(mt_len):
                psi = mt['rlnAnglePsi'][ptcl_index]
                xsh = mt['rlnOriginXAngst'][ptcl_index] / apix
                ysh = mt['rlnOriginYAngst'][ptcl_index] / apix
                t = EMAN2.Transform()
                t.set_params({'type':'2d', 'alpha':psi, 'tx':xsh, 'ty':ysh})
                mt['particles'][ptcl_index].transform(t)

            #average along a sliding window
            for ptcl_index in range(mt_len):
                low, hi = helper_fns.get_window(ptcl_index, 3, 4, mt_len)
                avg = EMAN2.Averagers.get('mean')
                for i in range(low, hi):
                    avg.add_image(mt['particles'][i])
                avg_img = avg.finish()
                
                #transform segment average back to original position and save
                psi = mt['rlnAnglePsi'][ptcl_index]
                xsh = mt['rlnOriginXAngst'][ptcl_index] / apix
                ysh = mt['rlnOriginYAngst'][ptcl_index] / apix
                t = EMAN2.Transform()
                t.set_params({'type':'2d', 'alpha':-psi, 'tx':-xsh, 'ty':-ysh})
                avg_img.transform(t)
                avg_img.append_image(outfile)
        #update image location    
        for index, name in enumerate(mgph['rlnImageName']):
            mgphname = re.search(name_regex, name).group(1)
            mgphname = mgphname.replace('.mrcs', '_SPs.mrcs')
            pnum = format_particle_number('000000', str(index+1))
            image_name = pnum + '@' + args.o + '/Micrographs/' + mgphname
            image_names.append(image_name)
            
    particles['rlnImageName'] = image_names
    starfile_data.add_datablock('data_particles', particles)
    starfile_data.write_star('%s/segment_averages.star' % args.o)
    print('\nFinished! Wrote %s/segment_averages.star' % args.o)
