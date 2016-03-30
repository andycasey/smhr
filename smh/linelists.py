from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import numpy as np
from astropy.io import ascii,fits
from astropy.table import Table, Column, MaskedColumn
from astropy import table
from utils import element_to_species, species_to_element

class LineList():
    full_colnames = ['wavelength','species','expot','loggf','damp_vdw','dissoc_E','comments',
                     'E_hi','E_lo','lande_hi','lande_lo','damp_stark','damp_rad','references']
    full_dtypes = [np.float,np.float,np.float,np.float,np.float,np.float,str,
                   np.float,np.float,np.float,np.float,np.float,np.float,str]
    moog_colnames = ['wavelength','species','expot','loggf','damp_vdw','dissoc_E','comments','references']
    moog_dtypes = [np.float,np.float,np.float,np.float,np.float,np.float,str,str]

    def __init__(self,data,verbose=False):
        # check required columns in data
        self.data = data
        self.verbose = verbose
        self.default_thresh = 0.01

    def merge(self,new_ll,thresh=None,override_current=True):
        """
        new_ll: 
            new LineList object to merge into this one

        thresh: 
            threshold for wavelength check when matching lines
            Defaults to self.default_thresh (0.01)

        override_current: 
            If True (default), uses new lines whenever duplicate lines are found.
            If False, keep current lines whenever duplicate lines are found.
        """
        if thresh==None: thresh = self.default_thresh

        num_in_list = 0
        lines_to_add = []
        lines_with_multiple_matches = []
        for new_line in new_ll.data:
            index = self.find_match(new_line,thresh)
            if index==-1:
                lines_to_add.append(new_line)
            elif index < -1:
                # use self.pick_best_line to find best line
                lines_with_multiple_matches.append(new_line)
                index = self.pick_best_line(new_line,thresh)
                if override_current:
                    self.data[index] = new_line
            elif index >= 0:
                num_in_list += 1
                if override_current:
                    self.data[index] = new_line
        num_lines_added = len(lines_to_add)
        if len(lines_to_add) > 0:
            new_lines = Table(rows=lines_to_add,names=lines_to_add[0].colnames)
            self.data = table.vstack([self.data,new_lines])
        
        if self.verbose:
            print("Num lines added: {}".format(num_lines_added))
            print("Num lines {}: {}".format('replaced' if override_current else 'ignored', num_in_list))
            num_lines_with_multiple_matches = len(lines_with_multiple_matches)
            print("Num lines with multiple matches: {}".format(num_lines_with_multiple_matches))
        
    def pick_best_line(self,new_line,thresh):
        """
        Given a line and assuming there are multiple matches,
        """
        indices = self.find_match(new_line,thresh=thresh,return_multiples=True)
        matches = self.data[indices]
        if self.verbose:
            print("----{} Matches----".format(len(matches)))
            print(new_line)
            print(matches)
        # Pick the line that is closest in wavelength
        best = np.argmin(np.abs(new_line['wavelength']-matches['wavelength']))
        return indices[best]

    def find_match(self,line,thresh=None,return_multiples=False):
        """
        See if the given line is in this line list
        Conditions: 
        (1) species match
        (2) wavelengeth match to within thresh
        (3) expot match to within 0.01 (hardcoded)

        Returns the index if all conditions are true

        return_multiples: 
            if True, returns list/array of indices
            if False, Returns negative number for no or multiple matches
            -1: no matches
            < -1: that number of matches
        """
        if thresh==None: thresh = self.default_thresh
        ii1 = self.data['species']==line['species']
        ii2 = np.abs(self.data['wavelength']-line['wavelength']) < thresh
        ii3 = np.abs(self.data['expot']-line['expot']) < 0.01
        ii = np.logical_and(np.logical_and(ii1,ii2),ii3)
        num_match = np.sum(ii)
        if num_match==0: return -1
        if num_match==1: return np.where(ii)[0]
        if return_multiples:
            return np.where(ii)[0]
        else:
            if self.verbose:
                print("Error: {} matches!".format(num_match))
                print(line)
            return -1 * num_match

    def find_duplicates(self,thresh=None):
        # The idea here is that you can increase the threshold to see if you were too weak in finding duplicates
        # This is not useful if you have molecular lines (e.g. carbon) because there are too many collisions
        # TODO untested
        if thresh==None: thresh = self.default_thresh
        duplicate_indices = []
        mask = np.ones(len(self.data),dtype=bool)
        # This is really inefficient and can likely be redone
        for i,line in enumerate(self.data):
            mask[i] = False
            tdata = LineList(self.data[mask])
            if tdata.find_match(line) >= 0: duplicate_indices.append(i)
            mask[i] = True
        duplicate_lines = self.data[np.array(duplicate_indices)]
        return duplicate_indices,duplicate_lines

    @staticmethod
    def is_line(line):
        # Check if the given thing is a line
        # Basically it has to be an astropy row with the right columns
        # If yes, return True
        # If not, print an error as to why
        # TODO
        raise NotImplementedError

    @classmethod
    def read(cls,filename):
        for reader in [cls.read_moog, cls.read_GES]:
            try:
                return reader(filename)
            except IOError:
                continue
        raise IOError("Cannot identify linelist format")

    @classmethod
    def read_moog(cls,filename,full_columns=True):
        if full_columns:
            colnames = cls.full_colnames
            dtypes = cls.full_dtypes
        else:
            colnames = cls.moog_colnames
            dtypes = cls.moog_dtypes
    
        with open(filename) as f:
            lines = f.readlines()
        N = len(lines)
        wl = np.zeros(N)
        species = np.zeros(N)
        EP = np.zeros(N)
        loggf = np.zeros(N)
        damping = np.zeros(N) * np.nan
        dissoc = np.zeros(N) * np.nan
        comments = ['' for i in range(N)]
        # wl, transition, EP, loggf, VDW damping C6, dissociation D0, EW, comments
        has_header_line = False
        for i,line in enumerate(lines):
            line = line.strip()
            s = line.split()
            try:
                _wl,_species,_EP,_loggf = map(float,s[:4])
            except:
                if i==0:
                    has_header_line = True
                    continue
                else:
                    raise ValueError("Invalid line: {}".format(line))
            if len(s) > 4:
                try: damping[i] = float(line[40:50])
                except: pass
                
                try: dissoc[i] = float(line[50:60])
                except: pass
                
                comments[i] = line[60:].strip()
            wl[i] = _wl; species[i] = _species; EP[i] = _EP; loggf[i] = _loggf
        if has_header_line:
            wl = wl[1:]; species = species[1:]; EP = EP[1:]; loggf = loggf[1:]
            damping = damping[1:]; dissoc = dissoc[1:]; comments = comments[1:]
        
        # check if gf by assuming there is at least one line with loggf < 0
        if np.all(loggf >= 0): 
            loggf = np.log10(loggf)
            # TODO this is the MOOG default, but it may not be a good idea...
            print("Warning: no lines with loggf < 0 in {}, assuming input is gf".format(filename))
        
        # Cite the filename as the reference for now
        # TODO
        refs = [filename for x in wl]
    
        # Fill required non-MOOG fields with nan
        if full_columns:
            nans = np.zeros_like(wl)*np.nan
            data = [wl,species,EP,loggf,damping,dissoc,comments,
                    nans,nans,nans,nans,nans,nans,refs]
        else:
            data = [wl,species,EP,loggf,damping,dissoc,comments,refs]
    
        return cls(Table(data,names=colnames,dtype=dtypes))

    @classmethod
    def read_GES(cls,filename):
        raise NotImplementedError

    def write_moog(self,filename):
        #TODO untested
        fmt = "{:10.3f}{:10.5f}{:10.3f}{:10.3f}{}{}{}{}"
        space = " "*10
        with open(filename,'w') as f:
            for line in self.data:
                C6 = space if np.ma.is_masked(line['damp_vdw']) else "{:10.3}".format(line['damp_vdw'])
                D0 = space if np.ma.is_masked(line['dissoc_E']) else "{:10.3}".format(line['dissoc_E'])
                comments = '' if np.ma.is_masked(line['comments']) else line['comments']
                f.write(fmt.format(line['wavelength'],line['species'],line['expot'],line['loggf'],C6,D0,space,line['comments'])+"\n")

    def write_latex(self,filename):
        #TODO untested
        #TODO rename columns to something nice?
        #TODO restrict columns?
        self.data.write(filename,format='ascii.latex')

    #TODO make "in" operator that calls self.contains so can do "line in ll"
    #TODO make "[]" operator that accesses the data columns?
    #TODO make iterable so "for line in ll" works
