from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import numpy as np
import astropy.units as u
from astropy.io import ascii,fits
from astropy.table import Table, Column, MaskedColumn
from astropy import table
from utils import element_to_species, species_to_element
import os

class LineList(Table):
    full_colnames = ['wavelength','species','expot','loggf','damp_vdw','dissoc_E','comments',
                     'E_hi','E_lo','lande_hi','lande_lo','damp_stark','damp_rad','references',
                     'element']
    full_dtypes = [np.float,np.float,np.float,np.float,np.float,np.float,str,
                   np.float,np.float,np.float,np.float,np.float,np.float,str,
                   str]
    moog_colnames = ['wavelength','species','expot','loggf','damp_vdw','dissoc_E','comments',
                     'references','element']
    moog_dtypes = [np.float,np.float,np.float,np.float,np.float,np.float,str,
                   str,str]

    def __init__(self,*args,**kwargs):
        # Pull out some default kwargs
        if 'verbose' in kwargs: 
            self.verbose = kwargs.pop('verbose')
        else:
            self.verbose = False
        if 'moog_columns' in kwargs: 
            self.moog_columns = kwargs.pop('moog_columns')
        else:
            self.moog_columns = False
        if 'default_thresh' in kwargs:
            self.default_thresh = kwargs.pop('default_thresh')
        else:
            self.default_thresh = 0.01

        super(LineList, self).__init__(*args,**kwargs)

        self.validate_colnames()

    def validate_colnames(self,error=True):
        ## This has to be included b'c table.vstack() creates an empty list
        if len(self.columns)==0: return False

        if self.moog_columns:
            colnames = self.moog_colnames
        else:
            colnames = self.full_colnames
        badcols = []
        for col in colnames:
            if col not in self.columns: badcols.append(col)
        
        error_msg = "Missing columns: {}".format(badcols)
        if len(badcols) == 0: return True
        if error:
            raise IOError(error_msg)
        else:
            print(error_msg)
        return False

    def merge(self,new_ll,thresh=None,override_current=True,
              in_place=True):
        """
        new_ll: 
            new LineList object to merge into this one

        thresh: 
            threshold for wavelength check when matching lines
            Defaults to self.default_thresh (0.01)

        override_current: 
            If True (default), uses new lines whenever duplicate lines are found.
            If False, keep current lines whenever duplicate lines are found.
        
        in_place:
            If True, merge new lines into this object
            If False, return a new LineList
        """
        if thresh==None: thresh = self.default_thresh

        num_in_list = 0
        lines_to_add = []
        lines_with_multiple_matches = []
        for new_line in new_ll:
            index = self.find_match(new_line,thresh)
            if index==-1:
                lines_to_add.append(new_line)
            elif index < -1:
                # use self.pick_best_line to find best line
                lines_with_multiple_matches.append(new_line)
                index = self.pick_best_line(new_line,thresh)
                # index < 0 is the convention that you should just skip the line rather than overwriting
                if index < 0: continue 
                if override_current:
                    #self.data[index] = new_line
                    self[index] = new_line
            elif index >= 0:
                num_in_list += 1
                if override_current:
                    self[index] = new_line
        num_lines_added = len(lines_to_add)
        if len(lines_to_add) > 0:
            if in_place:
                for line in lines_to_add:
                    self.add_row(line)
            else:
                new_lines = Table(rows=lines_to_add,names=lines_to_add[0].colnames)
                old_lines = self.copy()
                # During the vstack creates an empty LineList and warns
                new_data = table.vstack([old_lines,new_lines])
        
        if self.verbose:
            print("Num lines added: {}".format(num_lines_added))
            print("Num lines {}: {}".format('replaced' if override_current else 'ignored', num_in_list))
            num_lines_with_multiple_matches = len(lines_with_multiple_matches)
            print("Num lines with multiple matches: {}".format(num_lines_with_multiple_matches))

        if not in_place:
            return LineList(new_data)
        else:
            return None
        
    def pick_best_line(self,new_line,thresh):
        """
        Given a line and assuming there are multiple matches, pick the best line.
        By default picks line closest in wavelength.

        The convention is to return -1 if you want to skip the line.
        (This is so if you replace this function with some sort of interactive
        line picking, you can choose to not replace any lines.)
        """
        indices = self.find_match(new_line,thresh=thresh,return_multiples=True)
        if isinstance(indices,int): return -1 #Only one match, skip
        assert len(indices) >= 2
        matches = self[indices]
        # If there is an identical line, return -1 to skip
        for line in matches:
            dwl = np.abs(new_line['wavelength']-line['wavelength'])
            dEP = np.abs(new_line['expot']-line['expot'])
            dgf = np.abs(new_line['loggf']-line['loggf'])
            # TODO does this make sense?
            if dwl < .001 and dEP < .01 and dgf < .001: 
                if self.verbose:
                    print("Found identical match: {:8.3f} {:4.1f} {:5.2f} {:6.3f}".format(new_line['wavelength'],new_line['species'],new_line['expot'],new_line['loggf']))
                return -1

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
        ii1 = self['species']==line['species']
        ii2 = np.abs(self['wavelength']-line['wavelength']) < thresh
        ii3 = np.abs(self['expot']-line['expot']) < 0.01
        ii = np.logical_and(np.logical_and(ii1,ii2),ii3)
        num_match = np.sum(ii)
        if num_match==0: return -1
        if num_match==1: return np.where(ii)[0]
        if return_multiples:
            return np.where(ii)[0]
        else:
            return -1 * num_match

    def find_duplicates(self,thresh=None):
        # TODO test
        # The idea here is that you can increase the threshold to see if you were too weak in finding duplicates
        # This is not useful if you have molecular lines (e.g. carbon) because there are too many collisions
        if thresh==None: thresh = self.default_thresh
        duplicate_indices = []
        mask = np.ones(len(self),dtype=bool)
        # This is really inefficient and can likely be redone
        for i,line in enumerate(self):
            mask[i] = False
            tdata = LineList(self[mask])
            if tdata.find_match(line) >= 0: duplicate_indices.append(i)
            mask[i] = True
        duplicate_lines = self[np.array(duplicate_indices)]
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
    def read(cls,filename,*args,**kwargs):
        if not os.path.exists(filename):
            raise IOError("No such file or directory: {}".format(filename))
        for reader in [cls.read_moog, cls.read_GES]:
            try:
                return reader(filename)
            except IOError as e:
                continue
        #TODO this last part is untested
        try:
            return cls(super(LineList, cls).read((filename,)+args, **kwargs))
        except IOError:
            pass
        raise IOError("Cannot identify linelist format")

    @classmethod
    def read_moog(cls,filename,moog_columns=False):
        if moog_columns:
            colnames = cls.moog_colnames
            dtypes = cls.moog_dtypes
        else:
            colnames = cls.full_colnames
            dtypes = cls.full_dtypes
    
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
                    raise IOError("Invalid line: {}".format(line))
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
        
        # TODO
        # Cite the filename as the reference for now
        refs = [filename for x in wl]
    
        # Element to species
        spec2elem = {}
        for this_species in np.unique(species):
            spec2elem[this_species] = species_to_element(this_species)
        elements = [spec2elem[this_species] for this_species in species]

        # Fill required non-MOOG fields with nan
        if moog_columns:
            data = [wl*u.angstrom,species,EP*u.eV,loggf,damping,dissoc*u.eV,comments,refs,elements]
        else:
            nans = np.zeros_like(wl)*np.nan
            data = [wl,species,EP,loggf,damping,dissoc,comments,
                    nans,nans,nans,nans,nans,nans,refs,elements]
    
        return cls(Table(data,names=colnames,dtype=dtypes),moog_columns=moog_columns)

    @classmethod
    def read_GES(cls,filename):
        raise IOError("Not implemented")

    def write_moog(self,filename):
        #TODO untested
        fmt = "{:10.3f}{:10.5f}{:10.3f}{:10.3f}{}{}{}{}"
        space = " "*10
        with open(filename,'w') as f:
            for line in self:
                C6 = space if np.ma.is_masked(line['damp_vdw']) else "{:10.3}".format(line['damp_vdw'])
                D0 = space if np.ma.is_masked(line['dissoc_E']) else "{:10.3}".format(line['dissoc_E'])
                comments = '' if np.ma.is_masked(line['comments']) else line['comments']
                f.write(fmt.format(line['wavelength'],line['species'],line['expot'],line['loggf'],C6,D0,space,line['comments'])+"\n")

    #def write_latex(self,filename):
    #    #TODO untested
    #    #TODO rename columns to something nice?
    #    #TODO restrict columns?
    #    self..write(filename,format='ascii.latex')

    #TODO make "in" operator that calls self.contains so can do "line in ll"
    #TODO make "[]" operator that accesses the data columns?
    #TODO make iterable so "for line in ll" works
