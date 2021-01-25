#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MiRP - a microtubule RELION-based pipeline for cryo-EM image processing. 
MiRP v2 is designed to function with RELION v3.1, and is not compatible with earlier versions of RELION.
"""

__author__ = 'Alexander D. Cook & Joseph Atherton'
__license__ = 'GPLv3'
__version__ = '2.0'


import starfileIO
import helper_fns
from scipy import stats
import numpy as np
import collections
import itertools
import operator
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import math
import os
import warnings

class Microtubules:

    def __init__(self, starfile_in, job_path):
        #check if in RELION directory, and setup output path and standard out
        assert os.path.exists('default_pipeline.star'), 'default_pipeline.star not found. Please execute in a RELION directory'
        self.job_path = job_path
        self.stdout = '%s/run.out' % self.job_path
        self.outfile = None
        
        #read in data_particles datablock from _data.star type file
        self.starfile_in = starfile_in
        self.starfile_data =  starfileIO.Starfile(self.starfile_in)
        self.starfile_data.read_star()
        self.apix = float(self.starfile_data.get_entry('data_optics', 'rlnImagePixelSize')[0])
        #split data into microtubule blocks (list of dictionaries)
        self._data = self._get_microtubules()
        self.mt_tot = len(self._data)
        #supress scipy warnings with linear regression
        warnings.filterwarnings('ignore')


    
    ###### I/O methods ######
    #Add text to the RELION run.out file, for display in the GUI
    def _add_stdout(self, content, replace):
        #option to replace the last line of the file - used for progress bar
        if replace:
            with open(self.stdout, 'r') as f:
                new_file = f.readlines()[:-1]
                new_file.append(content)
            with open(self.stdout, 'w') as f:
                for l in new_file:
                    f.write(l)
        #directly addend output to standard out
        else:
            with open(self.stdout, 'a') as f:
                f.write(content)

    #Take RELION v3.1 stafile, convert particle datablock (dictionary of lists) to microtubules (list of dictionaries)
    def _get_microtubules(self):
        self.starfile_data.sort_loop_datablock('data_particles', 'rlnMicrographName', 'rlnHelicalTubeID', 'rlnHelicalTrackLengthAngst')
        db = self.starfile_data.get_datablock('data_particles')
        return helper_fns.group_dict_of_list(db, 'rlnMicrographName', 'rlnHelicalTubeID')

    #convert microtubules to particle datablock and save to specified starfile. Destructive of original starfile data
    def _write_microtubules(self, name):
        data = self._microtubules_to_particles(self._data)
        self.starfile_data.add_datablock('data_particles', data)
        self.starfile_data.write_star(name)



    
    ###### Protofilament number correction ######
    def vote_pf_number(self, cutoff):
        self._add_stdout('\nMiRP - voting on protofilament number for microtubules in %s...\n\n' %  self.starfile_in, False)
        cutoff = float(cutoff)
        split_mts = []
        confidence_data = []
        plot_pdf = PdfPages('%s/protofilament_corrected.pdf' % self.job_path)                

        #for each microtubule find the most common (modal) class number 
        for ix, microtubule in enumerate(self._data):
            self._add_stdout('Correcting microtubule %i of %i' % (ix, self.mt_tot), True)
            uncorr_data = microtubule['rlnClassNumber']
            mts_to_plot = {'uncorrected':uncorr_data, 'corrected':[]}
            
            #method to split microtubules where a significant switch in class assignment occurs
            #smoothen class numbers by taking the mode for each particle over a seven particle window
            smoothened_data = self._mode_smoothen(uncorr_data)
            #if there is a change in class number, split the microtubule microtubule at the index of the change
            nmts = self._split_mt_on_change(microtubule, smoothened_data)

            #for each 'new' microtubule (however, if data is good, most microtuubles should not be split) determine protofilament number
            for mt in nmts:
                #vote on the modal class assigmnment
                confidence = self._vote_mode(mt, 'rlnClassNumber')
                #remove very short microtubules
                if self._microtubule_len(mt) < 5:
                    pass
                #remove microtubules with a low confidence in class assignment
                elif confidence < cutoff:
                    confidence_data.append(confidence)
                else:
                    confidence_data.append(confidence)
                    mts_to_plot['corrected'].append(mt['rlnClassNumber'])
                    #renumber microtubules, since split microtubules increase the total
                    mt = self._renumber_tube_id(split_mts, mt)
                    split_mts.append(mt)
            self._plot_pf_number_corrected(mts_to_plot, plot_pdf)
        plot_pdf.close()
        self._plot_confidence(confidence_data, cutoff)
        
        #sort and separate microtubules by class number, and save to separate starfiles for each protofilament number
        particles = self._microtubules_to_particles(split_mts)
        particles = helper_fns.sort_dict_of_list(particles, 'rlnClassNumber')   
        pf_seg = helper_fns.group_dict_of_list(particles, 'rlnClassNumber')
        self._add_stdout('\nWrote ', False)
        for pf in pf_seg:
            self.starfile_data.add_datablock('data_particles', pf)
            self.starfile_data.sort_loop_datablock('data_particles', 'rlnMicrographName', 'rlnHelicalTubeID', 'rlnHelicalTrackLengthAngst')
            fname = '%s/1%ipf_data.star' % (self.job_path, self.starfile_data.get_entry('data_particles', 'rlnClassNumber')[0])
            self.starfile_data.write_star(fname)
            self._add_stdout('%s, ' % fname, False)
        self._pf_number_stats(self._data, split_mts)

    #function for renumbering expanded microtubules in micrograph (where microtubules have been split, e.g. based on protofilament number changes) 
    def _renumber_tube_id(self, mts, mt):
        try:
            prev_mgph = mts[-1]['rlnMicrographName'][0]
            curr_mgph = mt['rlnMicrographName'][0]
            if curr_mgph == prev_mgph:
                prev_tubeid = int(mts[-1]['rlnHelicalTubeID'][0])
                curr_tubeid = prev_tubeid + 1
            else:
                curr_tubeid = 1
        except IndexError:
            curr_tubeid = 1
        
        mt['rlnHelicalTubeID'] = [curr_tubeid for _ in range(self._microtubule_len(mt))]
        return mt

    #Calcualte the percentage of different protofilament numbers in uncorrected and corrected data. 
    def _pf_number_stats(self, uncorr_mts, corr_mts):

        def per(frac, total):
            return (frac / total) * 100

        uncorr_total_mts = len(uncorr_mts)
        uncorr_total_ptcls = self._get_total_particle_number(uncorr_mts)
        corr_total_ptcls = self._get_total_particle_number(corr_mts)
        corr_total_mts = len(corr_mts)

        uncorr_class = collections.Counter(self._get_global_data(uncorr_mts, 'rlnClassNumber'))
        corr_class = collections.Counter(self._get_global_data(corr_mts, 'rlnClassNumber'))
        
        stats = starfileIO.Starfile('%s/pf_number_sorting_stats.star' % self.job_path)
        gen = {'rlnTotalNumberTubes': uncorr_total_mts,
               'mirpTotalNumberTubes': corr_total_mts,
               'mirpPredictedChangesInPfNumber': corr_total_mts - uncorr_total_mts}
        pf = {'mtProtofilamentNumber': [11, 12, 13, 14, 15, 16],
              'rlnClassDistribution': [per(uncorr_class[i], uncorr_total_ptcls) for i in range(1, 7)],
              'mirpClassDistribution': [per(corr_class[i], corr_total_ptcls) for i in range(1, 7)]}
        stats.add_datablock('data_general', gen)
        stats.add_datablock('data_percent_protofilament_number', pf)
        stats.write_star('%s/pf_number_sorting_stats.star' % self.job_path)

        changes = stats.get_entry('data_general', 'mirpPredictedChangesInPfNumber')        
        pfnums = stats.get_entry('data_percent_protofilament_number', 'mtProtofilamentNumber')
        p_pfnum = stats.get_entry('data_percent_protofilament_number', 'mirpClassDistribution')
        self._add_stdout('\n\nTotal predicted changes in protofilament number: %i' % changes, False)
        self._add_stdout('\nData composed of:', False)
        for pf, per in zip(pfnums, p_pfnum):
            self._add_stdout(' %i %iPF ' % (per, pf), False)
        self._add_stdout(' microtubules\n', False)

        return stats

    #make plots for each microtubule showing uncorrected and corrected data side by side
    def _plot_pf_number_corrected(self, mts_to_plot, pdfpages):
        uncorr_yax = mts_to_plot['uncorrected']
        uncorr_xax = [i for i in range(1, len(uncorr_yax)+1 )]
        
        corr_yaxs = mts_to_plot['corrected']
        corr_xaxs = [ [i for i in range(1, len(y)+1 )] for y in corr_yaxs]
        
        total_plots = len(corr_yaxs) + 1
        fig=plt.figure()
        ax1 = fig.add_subplot(1, total_plots+1, 1)
        ax1.set_title('Uncorrected')
        ax1.plot(uncorr_xax, uncorr_yax, 'o')
        ax1.set_ylim([1,6])
        ax1.set_ylabel('Protofilament Class Number')
        ax1.set_xticks = uncorr_xax

        for i in range(1, total_plots):
            ax = fig.add_subplot(1, total_plots+1, i+1, sharey=ax1)
            ax.set_title('Corrected%i' % i)
            ax.set_xlabel = 'Particle Number'
            ax.set_xticks = corr_xaxs[i-1]
            ax.plot(corr_xaxs[i-1], corr_yaxs[i-1], 'o')

        plt.tight_layout()
        pdfpages.savefig(fig)
        plt.close()



    ###### Rot angle correction ######
    def vote_on_rot(self):
        self._add_stdout('\nMiRP - voting on Rotation angle for microtubules in %s...\n\n' %  self.starfile_in, False)
        cutoff = 8
        corrected_mts = []
        mts_to_remove = 0
        confidence = []
        plot_pdf = PdfPages('%s/rotation_corrected.pdf' % self.job_path)                

        #for each microtubule, find the most commonly assigned (modal) Rot angle, whilst accounting for the slope of microtubule supertwist
        for ix, microtubule in enumerate(self._data):
            self._add_stdout('Correcting microtubule %i of %i' % (ix, self.mt_tot), True)
            rot_angles = microtubule['rlnAngleRot']
            modal_clust, outliers = self._cluster_shallow_slopes(rot_angles, cutoff)
            #remove any microtubules for which clusters cannot be find (low particle number microtubules, or with widely distributed Rot angles)
            if not modal_clust:
                mts_to_remove += 1
            else:
                #correct the particle Rot angles to follow the fitted straight line of the modal Rot angle cluster
                #calculate confidence in Rot angle assignment (will always be low)
                c = len(modal_clust) / self._microtubule_len(microtubule) * 100
                confidence.append(c)
                fit = self._fit_eulerXY(rot_angles, modal_clust)
                microtubule['rlnAngleRot'] = fit
                microtubule['rlnAngleRotPrior'] = fit
                corrected_mts.append(microtubule)
                self._plot_rot_vote(outliers, rot_angles, microtubule['rlnAngleRot'], plot_pdf)
        plot_pdf.close()

        self._plot_confidence(confidence, 0)
        self._add_stdout('\n%s microtubules could not be fitted and were removed.' % mts_to_remove, False)
        self._data = corrected_mts
        self.outfile = '%srotCorrected_data.star' % self.job_path
        self._add_stdout('\nWrote %s' % self.outfile, False)
        self._write_microtubules(self.outfile)

    #for each microtubule, plot the uncorrected Rot angles, with straight lines demonstrating the clusters found
    #then plot the corrected Rot angle
    def _plot_rot_vote(self, outliers, uncorr, corr, pdfpages):
        xax = [i for i in range(1, len(uncorr)+1)]
        yticks = [i for i in range(-180, 180 + 1, 40) ]

        fig = plt.figure()
        ax1 = fig.add_subplot(1, 2, 1)
        ax1.set_title('Uncorrected + Clusters')
        ax1.plot(xax, uncorr, 'o')
        ax1.plot(xax, corr, '-', linewidth=2)
        for o in outliers:
            fitted = self._fit_eulerXY(uncorr, o)
            ax1.plot(xax, fitted, '-', linewidth=1)
        ax1.set_ylabel('Rot (o)')
        ax1.set_yticks = yticks
        ax1.set_ylim(-181, 181)
        ax1.set_xticks = xax
        ax1.set_xlabel('Particle Number')

        ax2 = fig.add_subplot(1, 2, 2)
        ax2.set_title('Corrected')
        ax2.plot(xax, corr, 'o')
        ax2.set_ylim(-181, 181)
        ax2.set_ylabel('Rot (o)')
        ax2.set_yticks = yticks
        ax2.set_xticks = xax
        ax2.set_xlabel('Particle Number')

        plt.tight_layout()
        pdfpages.savefig(fig)
        plt.close()



    ###### X/Y shift correction ######
    def vote_on_xy(self):
        self._add_stdout('\nMiRP - voting on X/Y shifts for microtubules in %s...\n\n' %  self.starfile_in ,False)
        plot_pdf = PdfPages('%s/XY_corrected.pdf' % self.job_path)                
        corrected_mts = []

        #for each microtubule pick the most populated linear region in the X/Y-shifts, and force all shifts to follow that line
        for ix, microtubule in enumerate(self._data):
            try:
                #remove these parameters as they can work against MiRP Rot angle assignment
                del microtubule['rlnAnglePsiFlipRatio']
            except KeyError:
                pass
            self._add_stdout('Correcting microtubule %i of %i' % (ix, self.mt_tot), True)
            Xsh, Ysh = microtubule['rlnOriginXAngst'], microtubule['rlnOriginYAngst']
            #find most populated linear region (modal cluster)
            Xmodal_clust = max(self._cluster_breaks(Xsh, 4), key = len)
            Ymodal_clust = max(self._cluster_breaks(Ysh, 4), key = len)
            if Xmodal_clust == [0]:
                Xcorr = [0 for x in range(1, len(Xsh)+1)]
                Ycorr = [0 for x in range(1, len(Ysh)+1)]
            else:
                #linear regression to get slope and y-intercept of modal cluster, and correct shifts based on this
                Xcorr = self._fit_eulerXY(Xsh, Xmodal_clust)
                Ycorr = self._fit_eulerXY(Ysh, Ymodal_clust)
            microtubule['rlnOriginXAngst'] = Xcorr
            microtubule['rlnOriginYAngst'] = Ycorr
            self._plot_xy_vote(Xsh, Ysh, Xcorr, Ycorr, plot_pdf)
            corrected_mts.append(microtubule)
        plot_pdf.close()

        self._data = corrected_mts        
        self.outfile = '%sxyCorrected_data.star' %  self.job_path
        self._add_stdout('\nWrote %s' % self.outfile, False)
        self._write_microtubules(self.outfile)    
    
    #for each microtubule, plot uncorrected and corrrected X/Y-shifts
    def _plot_xy_vote(self, Xuncorr, Yuncorr, Xcorr, Ycorr, pdfpages):
        xax = [i for i in range(1, len(Xuncorr)+1)]

        fig = plt.figure()
        ax1 = fig.add_subplot(1, 2, 1)
        ax1.set_title('Uncorrected')
        ax1.plot(xax, Xuncorr, 'o')
        ax1.plot(xax, Yuncorr, 'o')
        ax1.set_ylabel('X/Y Shifts (Angst)')
        ax1.set_xticks = xax
        ax1.set_xlabel('Particle Number')

        ax2 = fig.add_subplot(1, 2, 2)
        ax2.set_title('Corrected')
        ax2.plot(xax, Xcorr, 'o')
        ax2.plot(xax, Ycorr, 'o')
        ax2.set_ylabel('X/Y Shifts (Angst)')
        ax2.set_xticks = xax
        ax2.set_xlabel('Particle Number')
        
        plt.tight_layout()
        pdfpages.savefig(fig)
        plt.close()



    ###### Seam Checking ######
    def vote_on_seam(self, cutoff, pfnum, rise):
        self._add_stdout('\nMiRP - voting on relative seam position...\n\n', False)
        confidence_data = []
        new_mts = []
        cutoff = float(cutoff)
        pfnum = int(pfnum)
        rise = float(rise)
        
        #for each microtubule, calculate the modal class from 3D seam classification, and use this to correct the seam position relative to the 3D reference
        for ix, microtubule in enumerate(self._data):
            self._add_stdout('Correcting microtubule %i of %i' % (ix, self.mt_tot), True)
            mt_len = self._microtubule_len(microtubule)
            #get modal class for this microtubule, and confidence in class assignment 
            class_assignment = microtubule['rlnClassNumber']
            class_distribution = collections.Counter(class_assignment)
            top_class, top_class_freq = class_distribution.most_common(1)[0]
            confidence = top_class_freq / mt_len * 100

            #remove microtubules with lower confidence than the cutoff 
            if confidence < cutoff:
                confidence_data.append(confidence)
            else:
                # correct microtubules with alpha/beta-tubulin out of register
                if top_class > pfnum:
                    self._shift_along_z(microtubule, 41)  
                    # replace all class assignments with the modal class
                    microtubule['rlnClassNumber'] = [top_class - pfnum for _ in range(mt_len)]
                else:
                    microtubule['rlnClassNumber'] = [top_class for _ in range(mt_len)]
                confidence_data.append(confidence)
                # correct the rot angle based on the modal class
                self._correct_pfregister(pfnum, rise, microtubule)
                new_mts.append(microtubule)
        
        self._plot_confidence(confidence_data, self.job_path)
        self._plot_seam_stats(new_mts)
        self._data = new_mts
        self.outfile = '%sseamCorrected_data.star' % self.job_path
        self._write_microtubules(self.outfile)
        self._add_stdout('\nWrote %s' % self.outfile, False)

    def _plot_seam_stats(self, microtubules):
        #plot the percentage of particles in each relative seam position
        particles = self._microtubules_to_particles(microtubules)
        size = self._get_total_particle_number(microtubules)
        distribution = collections.Counter(particles['rlnClassNumber'])

        stats = []
        for k in distribution:
            p = distribution[k] / size * 100
            stats.append( (k, p) )
        stats.sort(key=operator.itemgetter(0))
        seamclass = [i[0] for i in stats]
        freq = [i[1] for i in stats]

        #write a star file describing the relative seam position distribution
        stats_star = starfileIO.Starfile('seamcorrection_stats.star')
        stats_star.add_datablock('seam_class_distribution', {'seamClassNumber': seamclass, 'percentDistribution': freq})
        stats_star.write_star('%s/seamcorrection_stats.star' % self.job_path)

        plt.bar([i for i in range(1, len(freq)+1)], freq)
        plt.xlabel('Seam Class')
        plt.ylabel('Percent of Data')
        plt.savefig('%s/seamclass_distribution.pdf' % self.job_path)
        plt.close()

    #use seam classification results to calculate seam postion relative to the 3D reference, and correct the Rot angle and X/Y shifts accordingly
    def _correct_pfregister(self, pfnum, rise, mt):
        #convert class number to relaive seam position (e.g. -6 to 7 for 14 protofilament microtubule)
        seampos = convert_pfnum_to_semicircle(mt['rlnClassNumber'][0], pfnum)
        twist = (360 / pfnum)
        #calculate the change in rotation angle needed to correct the seam position
        Drot = seampos * twist
        for i in range(self._microtubule_len(mt)):
            #correct the rotation angle to align the seam of the experimental particles with the seam of the reference
            rot = mt['rlnAngleRot'][i]
            mt['rlnAngleRot'][i] = rot + Drot
            try:
                mt['rlnAngleRotPrior'][i] = rot + Drot
            except KeyError:
                pass
            #correct the X/Y shifts to account for the change in Rot
            self._shift_along_z(mt, seampos*rise)

    #translate particles along the microtubule z-axis using psi angle and desired z-axis translation (hypotenuse) to edit the x/y shifts
    def _shift_along_z(self, mt, shift):
        for i in range(self._microtubule_len(mt)):
            psi = mt['rlnAnglePsi'][i]
            xsh = mt['rlnOriginXAngst'][i]
            ysh = mt['rlnOriginYAngst'][i]
            dx = shift * math.cos(-psi * math.pi / 180)
            dy = shift * math.sin(-psi * math.pi / 180)
            mt['rlnOriginXAngst'][i] = xsh + dx
            mt['rlnOriginYAngst'][i] = ysh + dy
            


    ###### Microtubule operations ######
    #takes any number of relion labels and sets their values to zero
    def reset_eulerxy(self, *rln_labels):
       for label in rln_labels:
            #for Tilt, values a set to 90 degrees
           if label == 'rlnAngleTilt':
               for mt in self._data:
                   mt_len = self._microtubule_len(mt)
                   mt['rlnAngleTilt'] = [90 for _ in range(mt_len)]
           #for Psi, angles are set to their prior (so that after round of protofilament classification, they can be reset to picking angle)
           elif label == 'rlnAnglePsi':
               for mt in self._data:
                   mt_len = self._microtubule_len(mt)
                   mt['rlnAnglePsi'] = [mt['rlnAnglePsiPrior'][i] for i in range(mt_len)]
           else:
               for mt in self._data:
                   mt_len = self._microtubule_len(mt)
                   mt[label] = [0 for _ in range(mt_len)]

    #for one data entry, get all the data from all the microtubules and return as a list
    def _get_global_data(self, microtubules, label):
        data = []
        for mt in microtubules:
            data += mt[label]
        return data

    #convert from list of dictionaries (microtubules) to particles (dictionary of lists)
    def _microtubules_to_particles(self, mts):
        data = {k:[] for k in mts[0].keys()}
        for mt in mts:
            for key in mt.keys():
                data[key] += list(mt[key])
        return data

    #return the number of particles in a given microtubule
    def _microtubule_len(self, microtubule):
        return len(microtubule['rlnHelicalTubeID'])

    #return to the total number of particles for the given list of microtubules
    def _get_total_particle_number(self, mts):
        data = self._microtubules_to_particles(mts)
        return len(data['rlnHelicalTubeID'])



    ###### Methods for microtubule angle, shift, and class correction ######
    #for a list of numbers (e.g. class assignment for a microtubule), smoothen the numbers by calculating the mode over a window of 7
    def _mode_smoothen(self, data):
        smoothened = []
        for i, _ in enumerate(data):
            l, h = helper_fns.get_window(i, 3, 4, len(data))
            smoothened.append(int(stats.mode(data[l:h])[0]))
        return smoothened

    #for data from a microtubule, split the microtubule where changes in data value occur
    def _split_mt_on_change(self, microtubule, data):
        changes = []
        for index, (curr, nxt) in enumerate( zip(data[:-1], data[1:]) ):
            if curr - nxt != 0:
                changes.append(index+1)
        return helper_fns.split_dict_of_list(microtubule, changes)

    #find the modal value for a data entry in a microtubule, and change all values in the data entry to this. Return the confidence in the mode
    def _vote_mode(self, microtubule, label):
        mt_len = self._microtubule_len(microtubule)
        mode, freq = collections.Counter(microtubule[label]).most_common(1)[0]
        microtubule[label] = [mode for _ in range(mt_len)]
        confidence = freq / mt_len * 100
        return confidence

    #for a list of y-values, extract the desired values statedd in xax_data_tofit, and perform linear regression on them. Fit and return all values to this equation.
    def _fit_eulerXY(self, ydata, xax_data_tofit):
        yax_data_tofit = [ydata[x] for x in xax_data_tofit]
        slope, yincept = stats.linregress(xax_data_tofit, yax_data_tofit)[0:2]
        return [yincept + x*slope for x in range(1, len(ydata)+1)]

    #finds clusters where 2D data follows many straight lines with shallow slopes (e.g. microutuble Rot angles)
    def _cluster_shallow_slopes(self, angles, cutoff):
        #create distance matrix (residual of each data point against all others)
        #returns index of all datapoint pairs that were within the cutoff
        linkMtrx = [ (i, i2) for ( (i, j), (i2, j2) ) 
                    in itertools.combinations( enumerate(angles), 2 ) 
                    if -cutoff <= float(j) - float(j2) <= cutoff]
        
        #creates non-redundant set of clusters of related datapoints 
        cluster = []
        while linkMtrx:
            node = linkMtrx[-1]
            for idx, pair in reversed( list( enumerate(linkMtrx[:-1]) ) ):
                if any( i == j for i, j in itertools.combinations(node+pair, 2) ):
                    node += pair
                    node = tuple( set(node) )
                    del linkMtrx[idx]
            del linkMtrx[-1]
            cluster.append( sorted(node) )
        
        if not cluster: 
            return None, None
        topclst = cluster.pop( max( enumerate( [len(i) for i in cluster] ), 
                                    key=operator.itemgetter(1) )[0])
        return topclst, cluster
    
    #find the most populated linear region in 2D data
    #by calculating residuals betwee neighbours, and creating clusters where the boundaries are defined by residuals that are outside a cutoff
    def _cluster_breaks(self, data, cutoff):
        diff = [curr - nxt for curr, nxt in zip(data[:-1], data[1:])]
        
        clusters = []
        clust = [0]
        for index, d in enumerate(diff):
            if -cutoff <= d <= cutoff:
                clust.append(index + 1)
            else:
                clusters.append(clust)
                clust = [index + 1]
        clusters.append(clust)
        return clusters

    #for a PF number of seam classification job, take the list of confidence in MiRP class correction for each microtubule, and plot as a histogram
    def _plot_confidence(self, confidence, cutoff):
        self._add_stdout('\nPlotting confidence to %s/confidence.pdf...' % self.job_path, False)
        fq = plt.hist(confidence, bins=10)[0]
        plt.vlines(cutoff, ymin=0.0, ymax=fq[-1], colors='red', label='cutoff')
        plt.xlabel('Confidence (percent in modal class)')    
        plt.xticks=[0,10,20,30,40,50,60,70,80,90,100]
        plt.ylabel('Frequency')    
        plt.savefig('%s/confidence.pdf' % self.job_path)
        plt.close()

    def __repr__(self):
        return 'Microtubule(%s)' % self.starfile_in
    
    def __str__(self):
        return 'Microtubules from %s' % self.starfile_in
    
    def __setitem__(self, key, value):
        self._data[key] = value
    
    def __getitem__(self, key):
        return self._data[key]


#for the current protofilament number (e.g. a number between 1 and 14), return a symmetry operator between e.g. -6 and 7
def convert_pfnum_to_semicircle(current_pfnumber, total_pfnumber):
    half_circle = math.ceil( float(total_pfnumber) / 2 )
    if current_pfnumber <= half_circle:
        return current_pfnumber - 1
    else:
        return -(total_pfnumber % current_pfnumber) - 1 

    
