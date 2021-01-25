#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MiRP - a microtubule RELION-based pipeline for cryo-EM image processing. 
MiRP v2 is designed to function with RELION v3.1, and is not compatible with earlier versions of RELION.
starfileIO.py provides functions for reading, editing, and writing starfiles. It supports starfiles with multiple datablocks,
and is therefore compatible with any version of RELION.
"""

__author__ = 'Alexander D. Cook & Joseph Atherton'
__license__ = 'GPLv3'
__version__ = '2.0'


from collections import OrderedDict
import helper_fns
import operator
import itertools

class Starfile:
    
    def __init__(self, starfile):
        self.starfile = starfile
        self._datablocks = OrderedDict()

    def read_star(self):
        loop = False 
                        
        for fields in helper_fns.readfile(self.starfile):
            
            #identify if at start of new datablock
            if fields[0].startswith('data_'):
                curr_datablock_id = fields[0]
                curr_datablock = OrderedDict()
                self._datablocks[curr_datablock_id] = curr_datablock
                loop = False
            
            #identify if in loop style datablock
            elif fields[0] == 'loop_':
                loop = True
        
            #identify if encountered a data label
            elif fields[0].startswith('_'):
                label = fields[0][1:]
                #if loop data, make dictionary of datalabels as keys
                if loop:
                    curr_datablock[label] = []
                #if key-val data, make dictionary of key-val data
                else:
                    data = helper_fns.literal_eval(fields[1])
                    curr_datablock[label] = data

            #if loop data and not encountered a data label, populate data label dictionary with data    
            elif loop:
                labels = curr_datablock.keys()
                assert len(labels) == len(fields), ('Error: Some data in starfile %s does not match the labels!\n%s\n%s\n' % (self.starfile, labels, data))
                for label, data in zip(labels, fields):
                    data = helper_fns.literal_eval(data)
                    curr_datablock[label].append(data)
            else:
                print('Error: could not understand this line in starfile %s:\n%s\n' % (self.starfile, fields))


    ###### Getting starfile data ######
    #return the data of a specified datablock
    def get_datablock(self, datablock_id):
        return self._datablocks[datablock_id]
    
    #return the data labels only from a specified datablock
    def get_labels(self, datablock_id):
        return self._datablocks[datablock_id].keys()
    
    #return the data values only from a specified datablock
    def get_data(self, datablock_id):
        return self._datablocks[datablock_id].keys()
    
    #return the data from a specified datablock, associated with a specified data label
    def get_entry(self, datablock_id, label):
        return self._datablocks[datablock_id][label]
    
    #return the number of data entries in loop type datablock
    def get_loopdatablock_len(self, datablock_id):
        key = list(self._datablocks[datablock_id].keys())[0]
        return len(self._datablocks[datablock_id][key])
    

    ###### Updating/adding starfile data ######
    # data must be in dictionary format, where the key is the starfile data label, and the value is the starfile data
    #add a new datablock to an empty or existing starfile
    def add_datablock(self, datablock_id, datablock):
        assert isinstance(datablock, dict), 'Datablock must be a dictionary!'
        self._datablocks[datablock_id] = datablock
        
    #append loop data to an existing loop type datablock.     
    def add_loop_data(self, datablock_id, data):
        assert isinstance(data, dict), 'Data must be a dictionary!'
        
        db = self._datablocks[datablock_id]
        assert len(data) == len(db), 'New loop data must have the same number of labels as in the current datablock.'
        for key in (data):
            assert key in db, 'New loop data must have the same labels as in the current datablock.'
            
        for key in data:
            val = data[key]
            if isinstance(val, list):
                db[key] += val
            else:
                db[key].append(val)
                
    #add data to existing key-val type datablock.
    def add_nonloop_data(self, datablock_id, data):
        assert isinstance(data, dict), 'Data must be a dictionary!'
        
        db = self._datablocks[datablock_id]
        for key in data:
            assert key not in db, 'New non-loop data cannot have a label that is already in the datablock. Use update_nonloop_data instead.'
            assert not isinstance (data[key], list), 'Non-loop data cannot be a list.'
            
        for key in data:
            val = data[key]
            db[key] = val
    
    #replace loop data from a specified datablock, at a specified index
    def update_loop_data(self, datablock_id, data, index):
        assert isinstance(data, dict), 'Data must be a dictionary!'
        
        db = self._datablocks[datablock_id]
        assert len(data) == len(db), 'New loop data must have the same number of labels as in the current datablock'
        for key in (data):
            assert key in db, 'New loop data must have the same labels as in the current datablock.'
            assert not isinstance (data[key], list), 'New loop data to update with cannot be a list.'
    
        for key in data:
            db[key][index] = data[key]
            
    #replace an existing entry in a key-val type datablock   
    def update_nonloop_data(self, datablock_id, data):
        assert isinstance(data, dict), 'Data must be a dictionary!'
        
        db = self._datablocks[datablock_id]
        for key in data:
            assert key in db, 'New non-loop data to update with must have an existing label in the datablock. Use add_nonloop_data instead.'
            assert not isinstance (data[key], list), 'Non-loop data cannot be a list.'
            
        for key in data:
            db[key] = data[key]
        
    #order the data in a loop type datablock based on one or more data
    def sort_loop_datablock(self, datablock_id, *labels):
        data = self.get_datablock(datablock_id)
        data = helper_fns.sort_dict_of_list(data, *labels)
        self.add_datablock(datablock_id, data)

    # save the current starfile data 
    def write_star(self, name):
         
        with open(name, 'w') as f:
        
            for key in self._datablocks:
                f.write('\n%s\n\n' % key)
                datablock = self._datablocks[key]
                
                labels = datablock.keys()
                data = list(datablock.values())
                if isinstance(data[0], list) or isinstance(data[0], tuple):
                    f.write('loop_\n')
                    for idx, label in enumerate(labels):
                        f.write('_%s\t#%i\n' % (label, idx+1))
                    for entry in zip(*data):
                        f.write('\t'.join(str(val) for val in entry) + '\n')
                else:
                    for label, entry in zip(labels, data):
                        f.write('_%s\t%s\n' % (label, str(entry)))


    def __repr__(self):
        return 'Starfile(%s)' % self.starfile
    
    def __str__(self):
        return 'Starfile object with filename: %s' % self.starfile
    
    def __setitem__(self, key, value):
        self._datablocks[key] = value

    def __getitem__(self, key):
        return self._datablocks[key]

