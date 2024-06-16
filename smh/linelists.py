from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import numpy as np
import astropy.units as u
from astropy.io import ascii,fits,registry
from astropy.table import Table, Column, MaskedColumn
from astropy import table
from .utils import element_to_species, species_to_element
from .utils import elems_isotopes_ion_to_species, species_to_elems_isotopes_ion

import logging
logger = logging.getLogger(__name__)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter(
    "%(asctime)s [%(levelname)-8s] %(message)s"))
logger.addHandler(handler)

import os

import hashlib

class LineListConflict(Exception):
    """Exception raised for merging conflicts"""
    def __init__(self,conflicts1,conflicts2):
        assert len(conflicts1) == len(conflicts2)
        self.conflicts1 = conflicts1
        self.conflicts2 = conflicts2
    def __str__(self):
        return "LineListConflict ({} conflicts)\nConflicts from current list\n{}\nConflicts from new list\n{}".format(len(self.conflicts1), repr(self.conflicts1), repr(self.conflicts2)) 

class LineList(Table):
    full_colnames = ['wavelength','species','expot','loggf','damp_vdw','dissoc_E','comments',
                     'numelems','elem1','isotope1','elem2','isotope2','ion',
                     'E_hi','lande_hi','lande_lo','damp_stark','damp_rad','references','element']
    full_dtypes = [float,float,float,float,float,float,str,
                   int,str,int,str,int,int,
                   float,float,float,float,float,str,str]
    moog_colnames = ['wavelength','species','expot','loggf','damp_vdw','dissoc_E','comments',
                     'numelems','elem1','isotope1','elem2','isotope2','ion',
                     'references','element']
    moog_dtypes = [float,float,float,float,float,float,str,
                   int,str,int,str,int,int,
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
            self.default_thresh = 0.1
        if 'default_loggf_thresh' in kwargs:
            self.default_loggf_thresh = kwargs.pop('default_loggf_thresh')
        else:
            self.default_loggf_thresh = 0.01
        if 'default_expot_thresh' in kwargs:
            self.default_expot_thresh = kwargs.pop('default_expot_thresh')
        else:
            self.default_expot_thresh = 0.01
        if "check_for_duplicates" in kwargs:
            # If you check for duplicates, you do not have duplicates
            # (because a ValueError is thrown otherwise)
            self.has_duplicates = ~kwargs.pop("check_for_duplicates")
        else:
            # By default, do NOT check for duplicates
            self.has_duplicates = True

        super(LineList, self).__init__(*args,**kwargs)

        #self.validate_colnames(False)

    def compute_hashes(self):
        return np.array([self.hash(line) for line in self])
    
    def check_for_duplicates(self):
        """
        Check for exactly duplicated lines. This has to fail because hashes
        are assumed to be unique in a LineList.
        Exactly duplicated lines may occur for real reasons, e.g. if there is
        insufficient precision to distinguish two HFS lines in the line list.
        In these cases, it may be okay to combine the two lines into one 
        total line with a combined loggf.
        """
        hashes = self.compute_hashes()

        if len(self) != len(np.unique(hashes)):
            error_msg = \
                "This LineList contains lines with identical hashes.\n" \
                "The problem is most likely due to completely identical lines\n" \
                "(e.g., because of insufficient precision in HFS).\n" \
                "If that is the case, it may be reasonable to combine the\n" \
                "loggf for the two lines into a single line.\n" \
                "We now print the duplicated lines:\n"
            fmt = "{:.3f} {:.3f} {:.3f} {:5} {}\n"
            total_duplicates = 0
            for i,hash in enumerate(hashes):
                N = np.sum(hashes==hash)
                if N > 1: 
                    line = self[i]
                    total_duplicates += 1
                    error_msg += fmt.format(line['wavelength'],line['expot'],line['loggf'],line['element'],hashes[i])
            raise ValueError(error_msg)
        self.has_duplicates = False
        return None
    def validate_colnames(self,error=False):
        """
        error: if True, raise error when validating.
            This is False by default because many astropy.table operations
            create empty or small tables.
        """
        ## This is included b'c table.vstack() creates an empty table
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
            logger.warn(error_msg)
        return False

    @staticmethod
    def vstack(tables, **kwargs):
        """
        Wraps astropy.table.vstack and returns a LineList
        """
        return LineList(table.vstack(tables, **kwargs))

    @staticmethod
    def identify_conflicts(ll1, ll2, 
                           skip_exactly_equal_lines=False,
                           skip_equal_loggf=False,
                           dwl_thresh=.001,dEP_thresh=.01,dgf_thresh=.001):
        """
        skip_exactly_equal_lines: if True, skips single-line conflicts that are identical hashes
        skip_equal_loggf: if True, skips single-line conflicts that are almost identical
        """
        if len(ll1) > len(ll2): #swap
            swap = True
            _ll = ll1; ll1 = ll2; ll2 = _ll
        else:
            swap = False
        matches = np.zeros((len(ll1),len(ll2)),dtype=bool)
        for i,line in enumerate(ll1):
            ii1 = np.logical_and(np.logical_and(ll2['elem1']==line['elem1'], ll2['elem2']==line['elem2']), ll2['ion']==line['ion'])
            ii2 = np.abs(ll2['wavelength']-line['wavelength']) < .1
            ii3 = np.abs(ll2['expot']-line['expot']) < .01
            matches[i,ii1&ii2&ii3] = True
        #indices of all matches
        conflicts = np.array(np.where(matches)).T
        used = np.zeros(len(conflicts),dtype=bool)
        equivalence_classes = [] #list of lists of indices of conflicts
        while np.sum(used) != len(used):
            next_conflict = np.min(np.where(~used)[0])
            eq_class = [next_conflict]
            used[next_conflict] = True
            while True:
                N_last = len(eq_class)
                for conflict in conflicts[eq_class]:
                    # find indices of all unused conflicts that match this conflict
                    new_ii = np.where(np.logical_and(~used, np.logical_or(conflict[0]==conflicts[:,0], conflict[1]==conflicts[:,1])))[0]
                    if len(new_ii)==0: continue
                    eq_class.extend(list(new_ii))
                    used[new_ii] = True
                if len(eq_class) == N_last: break
            equivalence_classes.append(eq_class)
        assert np.sum([len(eq_class) for eq_class in equivalence_classes]) == len(conflicts)
        equivalence_lines1 = []
        equivalence_lines2 = []
        for eq_class in equivalence_classes:
            this_conflicts = conflicts[eq_class]
            equivalence_lines1.append( ll1[np.unique(this_conflicts[:,0])] )
            equivalence_lines2.append( ll2[np.unique(this_conflicts[:,1])] )
        if skip_exactly_equal_lines or skip_equal_loggf:
            if skip_equal_loggf: 
                equal_fn = lambda x,y: LineList.lines_equal(x,y,dwl_thresh=dwl_thresh,
                                                            dEP_thresh=dEP_thresh,dgf_thresh=dgf_thresh)
            if skip_exactly_equal_lines: #overwrite skip_equal_loggf
                equal_fn = lambda x,y: LineList.lines_exactly_equal(x,y)

            _drop_indices = []
            for i,(x,y) in enumerate(zip(equivalence_lines1,equivalence_lines2)):
                if len(x)==1 and len(y)==1 and equal_fn(x[0],y[0]):
                    _drop_indices.append(i)
                elif len(x)==len(y): # check for identical HFS/molecule blocks
                    for _x, _y in zip(x,y):
                        if not equal_fn(_x, _y): break
                    else: #all equal
                        _drop_indices.append(i)
            equivalence_lines1 = [v for i,v in enumerate(equivalence_lines1) if i not in _drop_indices]
            equivalence_lines2 = [v for i,v in enumerate(equivalence_lines2) if i not in _drop_indices]
        if swap: return equivalence_lines2, equivalence_lines1
        return equivalence_lines1, equivalence_lines2

    def merge(self,new_ll,thresh=None,loggf_thresh=None,raise_exception=True,
              skip_exactly_equal_lines=False,
              skip_equal_loggf=False,
              override_current=False,in_place=True,
              add_new_lines=True,
              ignore_conflicts=False):
        """
        new_ll: 
            new LineList object to merge into this one

        thresh: 
            threshold for wavelength check when matching lines
            Defaults to self.default_thresh (0.1)

        loggf_thresh: 
            threshold for loggf check when finding identical lines
            Defaults to self.default_loggf_thresh (0.01)

        raise_exception:
            If True (default), finds all the conflicts and raises LineListConflict
              Note: if in_place == True, then it merges new lines BEFORE raising the error
            If False, uses self.pick_best_line() to pick a line to overwrite

        skip_exactly_equal_lines:
            If True, skips lines that have equal hashes during the merge
            If False (default), raises exception for duplicate lines

        skip_equal_loggf:
            If True, skips lines that are almost exactly equal during the merge
            If False (default), raises exception for duplicate lines

        override_current: 
            If True, uses new lines whenever duplicate lines are found.
            If False (default), keep current lines whenever duplicate lines are found.
            Ignored if raise_exception == True
        
        in_place:
            If True (default), merge new lines into this object. It will do so BEFORE throwing any LineListConflict exceptions!
            If False, return a new LineList
        
        add_new_lines:
            If True (default), add new lines when merging.
            If False, do not add new lines. This is to replace lines from a list without adding them.

        ignore_conflicts:
            If True, merge the linelists without checking for conflicts
            If False (default), check for conflicts during merge

        """
        if thresh==None: thresh = self.default_thresh
        if loggf_thresh==None: loggf_thresh = self.default_loggf_thresh
        if len(self)==0: 
            if not in_place:
                return new_ll.copy()
            else:
                n_cols = len(new_ll.colnames)
                names = new_ll.colnames
                dtype = [None] * n_cols
                self._init_indices = self._init_indices and new_ll._copy_indices
                self._init_from_table(new_ll, names, dtype, n_cols, True)
                return None

        if ignore_conflicts:
            if self.verbose:
                logger.warn("Ignoring conflicts: adding {} lines".format(len(new_ll)))
            if not in_place:
                return table.vstack([self, new_ll])
            else:
                #combined = table.vstack([self, new_ll])
                #names = combined.colnames
                #dtype = [None] * n_cols
                #self._init_indices = self._init_indices and combined._copy_indices
                #self._init_from_table(combined, names, dtype, n_cols, True)
                #return None
                raise NotImplementedError
        num_in_list = 0
        num_with_multiple_conflicts = 0
        lines_to_add = []

        for j,new_line in enumerate(new_ll):
            index = self.find_match(new_line,thresh)
            if index==-1: # New Line
                lines_to_add.append(new_line)
            elif raise_exception: # Record all conflicts later
                pass
            else: # use self.pick_best_line to find best line
                if index < -1:
                    num_with_multiple_conflicts += 1
                    index = self.pick_best_line(new_line,thresh)
                    # index < 0 is the convention that you should just skip the line rather than overwriting
                    if index < 0: continue 
                    if override_current:
                        self[index] = new_line
                elif index >= 0:
                    num_in_list += 1
                    if override_current:
                        self[index] = new_line
        num_lines_added = len(lines_to_add)
        if add_new_lines and len(lines_to_add) > 0:
            if in_place:
                for line in lines_to_add:
                    self.add_row(line)
            else:
                new_lines = Table(rows=lines_to_add,names=lines_to_add[0].colnames)
                old_lines = self.copy()
                # During the vstack creates an empty LineList and warns
                new_data = table.vstack([old_lines,new_lines])
        else:
            if not in_place:
                new_data = self.copy()
        
        # Note: if in_place == True, then it merges new lines BEFORE raising the exception
        if raise_exception:
            conflicts1,conflicts2 = self.identify_conflicts(self,new_ll,
                                                            skip_exactly_equal_lines=skip_exactly_equal_lines,
                                                            skip_equal_loggf=skip_equal_loggf,
                                                            dwl_thresh=thresh, dgf_thresh=loggf_thresh)
            if len(conflicts1) > 0:
                raise LineListConflict(conflicts1, conflicts2)
        if self.verbose:
            logger.info("Num lines added: {}".format(num_lines_added))
            logger.info("Num lines {}: {}".format('replaced' if override_current else 'ignored', num_in_list))
            logger.info("Num lines with multiple matches: {}".format(num_with_multiple_conflicts))

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
            if self.lines_equal(new_line,line):
                if self.verbose:
                    logger.info("Found identical match: {:8.3f} {:4.1f} {:5.2f} {:6.3f}".format(new_line['wavelength'],new_line['species'],new_line['expot'],new_line['loggf']))
                return -1

        if self.verbose:
            logger.info("----{} Matches----".format(len(matches)))
            logger.info(new_line)
            logger.info(matches)
        
        # Pick the line that is closest in wavelength
        best = np.argmin(np.abs(new_line['wavelength']-matches['wavelength']))
        return indices[best]


    @property
    def unique_elements(self):
        """ Return the unique elements that are within this line list. """

        elements = list(self["elem1"]) + list(self["elem2"])
        return list(set(elements).difference([""]))


    def find_match(self,line,thresh=None,return_multiples=False):
        """
        (This method is terrible and will probably be removed)
        See if the given line is in this line list
        Conditions: 
        (1) elem1, elem2, ion match (but isotopes do not have to!)
        (2) wavelengeth match to within thresh
        (3) expot match to within 0.01 (hardcoded)

        Returns the index if all conditions are true

        thresh:
            Wavelength tolerance to be considered identical lines
            (defaults to self.default_thresh, which is 0.1 by default)

        return_multiples: 
            if True, returns array of indices
            if False, Returns negative number for no or multiple matches
            -1: no matches
            < -1: that number of matches
        """
        if thresh==None: thresh = self.default_thresh
        ii1 = np.logical_and(np.logical_and(self['elem1']==line['elem1'], self['elem2']==line['elem2']), self['ion']==line['ion'])
        ii2 = np.abs(self['wavelength']-line['wavelength']) < thresh
        ii3 = np.abs(self['expot']-line['expot']) < self.default_expot_thresh
        ii = np.logical_and(np.logical_and(ii1,ii2),ii3)
        num_match = np.sum(ii)
        if num_match==0: return -1
        if num_match==1: return np.where(ii)[0]
        if return_multiples:
            return np.where(ii)[0]
        else:
            return -1 * num_match

    def find_duplicates(self,thresh=None):
        # The idea here is that you can increase the threshold to see if you were too weak in finding duplicates
        # This is not useful if you have molecular lines (e.g. carbon) because there are too many collisions
        if thresh==None: thresh = self.default_thresh
        duplicate_indices = []
        mask = np.ones(len(self),dtype=bool)
        # This is really inefficient and can likely be redone
        for i,line in enumerate(self):
            matches = self.find_match(line,thresh=thresh,return_multiples=True)
            if len(matches) > 1: duplicate_indices.append(i)
        duplicate_lines = self[np.array(duplicate_indices)]
        return duplicate_indices,duplicate_lines

    def remove_exact_duplicates(self, in_place=False):
        """
        Returns a linelist with only unique hashes from this linelist
        """
        #n_cols = len(new_ll.colnames)
        #names = new_ll.colnames
        #dtype = [None] * n_cols
        #self._init_indices = self._init_indices and new_ll._copy_indices
        #self._init_from_table(new_ll, names, dtype, n_cols, True)
        #return None
        if in_place: raise NotImplementedError

        uniq, ix = np.unique(self.compute_hashes(), return_index=True)
        return self[ix]


    @staticmethod
    def hash(line):
        s = "{:.3f}_{:.3f}_{:.3f}_{}_{}_{}_{}_{}".format(line['wavelength'],line['expot'],line['loggf'],
                                                         line['elem1'],line['elem2'],line['ion'],line['isotope1'],line['isotope2'])
        #return md5.new(s).hexdigest()
        return hashlib.md5(s.encode("utf-8")).hexdigest()

    @staticmethod
    def lines_equal(l1,l2,dwl_thresh=.001,dEP_thresh=.01,dgf_thresh=.001):
        dwl = np.abs(l1['wavelength']-l2['wavelength'])
        dEP = np.abs(l1['expot']-l2['expot'])
        dgf = np.abs(l1['loggf']-l2['loggf'])
        return dwl < dwl_thresh and dEP < dEP_thresh and dgf < dgf_thresh

    @staticmethod
    def lines_exactly_equal(l1,l2):
        #return LineList.lines_equal(l1,l2,dwl_thresh=1e-4,dEP_thresh=1e-4,dgf_thresh=1e-4)
        return LineList.hash(l1)==LineList.hash(l2)

    @classmethod
    def create_basic_linelist(cls,wavelength,species,expot,loggf, **kwargs):
        """
        Create a minimum LineList for the common case of having only the 
        4 minimum linelist information:
        wavelength, species, expot, loggf
        """
        assert len(wavelength) == len(species) == len(expot) == len(loggf)
        wavelength = np.array(wavelength)
        species = np.array(species)
        expot = np.array(expot)
        loggf = np.array(loggf)
        N = len(wavelength)
        nans = np.zeros(N) + np.nan
        empty = np.array(["" for x in range(N)])

        spec2element = {}
        spec2elem1= {}
        spec2elem2= {}
        spec2iso1 = {}
        spec2iso2 = {}
        spec2ion  = {}
        for this_species in np.unique(species):
            spec2element[this_species] = species_to_element(this_species)
            _e1, _e2, _i1, _i2, _ion = species_to_elems_isotopes_ion(this_species)
            spec2elem1[this_species] = _e1
            spec2elem2[this_species] = _e2
            spec2iso1[this_species] = _i1
            spec2iso2[this_species] = _i2
            spec2ion[this_species] = _ion
        numelems = np.array([2 if x >= 100 else 1 for x in species])
        elements = [spec2element[this_species] for this_species in species]
        elem1 = [spec2elem1[this_species] for this_species in species]
        elem2 = [spec2elem2[this_species] for this_species in species]
        isotope1 = [spec2iso1[this_species] for this_species in species]
        isotope2 = [spec2iso2[this_species] for this_species in species]
        ion  = [spec2ion[this_species] for this_species in species]

        # Fill required non-MOOG fields with nan
        data = [wavelength, species, expot, loggf, nans, nans, empty,
                numelems, elem1, isotope1, elem2, isotope2, ion, empty, elements, nans]
        columns = cls.moog_colnames + ['equivalent_width']
        dtypes = cls.moog_dtypes + [float]

        return cls(Table(data, names=columns, dtype=dtypes), moog_columns=True, **kwargs)

    @classmethod
    def read(cls,filename,*args,**kwargs):
        """
        filename: name of the file
        To use the default Table reader, must specify 'format' keyword.
        Otherwise, tries to read moog and then GES fits
        """
        
        if 'format' in kwargs: 
            return cls(super(LineList, cls).read(*((filename,)+args), **kwargs))

        if not os.path.exists(filename):
            raise IOError("No such file or directory: {}".format(filename))
        for reader in [cls.read_moog, cls.read_GES]:
            try:
                return reader(filename,**kwargs)
            except (IOError, KeyError, UnicodeDecodeError) as e:
                # KeyError: Issue #87
                # UnicodeDecodeError: read_moog fails this way for fits
                pass
        raise IOError("Cannot identify linelist format (specify format if possible)")

    @classmethod
    def read_moog(cls,filename,moog_columns=False,**kwargs):
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
        ew = np.zeros(N) * np.nan
        damping = np.zeros(N) * np.nan
        dissoc = np.zeros(N) * np.nan
        comments = ['' for i in range(N)]
        # wl, transition, EP, loggf, VDW damping C6, dissociation D0, EW, comments
        has_header_line = False
        for i,line in enumerate(lines):
            s = line.split()
            try:
                _wl,_species,_EP,_loggf = list(map(float,s[:4]))
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
                
                try: 
                    _ew = float(line[60:70])
                    # It seems some linelists have -1 in the EW location as a placeholder
                    if _ew <= 0: 
                        raise ValueError("EW <= 0: {}".format(_ew))
                    ew[i] = _ew
                    comments[i] = line[70:].strip()
                except:
                    comments[i] = line[60:].strip()
            wl[i] = _wl; species[i] = _species; EP[i] = _EP; loggf[i] = _loggf
        if has_header_line:
            wl = wl[1:]; species = species[1:]; EP = EP[1:]; loggf = loggf[1:]
            damping = damping[1:]; dissoc = dissoc[1:]; comments = comments[1:]
            ew = ew[1:]
        
        # check if gf by assuming there is at least one line with loggf < 0
        ## TURNED THIS OFF IT WAS CAUSING HUGE PROBLEMS
        #if np.all(loggf >= 0): 
        #    loggf = np.log10(loggf)
        #    # TODO this is the MOOG default, but it may not be a good idea...
        #    logger.warn("MOOG: no lines with loggf < 0 in {}, assuming input is gf".format(filename))
        
        # TODO
        # Cite the filename as the reference for now
        refs = [filename for x in wl]
    
        # Species to element
        spec2element = {}
        spec2elem1= {}
        spec2elem2= {}
        spec2iso1 = {}
        spec2iso2 = {}
        spec2ion  = {}
        for this_species in np.unique(species):
            spec2element[this_species] = species_to_element(this_species)
            _e1, _e2, _i1, _i2, _ion = species_to_elems_isotopes_ion(this_species)
            spec2elem1[this_species] = _e1
            spec2elem2[this_species] = _e2
            spec2iso1[this_species] = _i1
            spec2iso2[this_species] = _i2
            spec2ion[this_species] = _ion
        numelems = np.array([2 if x >= 100 else 1 for x in species])
        elements = [spec2element[this_species] for this_species in species]
        elem1 = [spec2elem1[this_species] for this_species in species]
        elem2 = [spec2elem2[this_species] for this_species in species]
        isotope1 = [spec2iso1[this_species] for this_species in species]
        isotope2 = [spec2iso2[this_species] for this_species in species]
        ion  = [spec2ion[this_species] for this_species in species]

        # Fill required non-MOOG fields with nan
        if moog_columns:
            data = [wl,species,EP,loggf,damping,dissoc,comments,
                    numelems,elem1,isotope1,elem2,isotope2,ion,
                    refs,elements]
        else:
            nans = np.zeros_like(wl)*np.nan
            E_hi = EP + 12398.42/wl #hc = 12398.42 eV AA
            data = [wl,species,EP,loggf,damping,dissoc,comments,
                    numelems,elem1,isotope1,elem2,isotope2,ion,
                    E_hi,nans,nans,nans,nans,refs,elements]
        # add EW if needed
        #if not np.all(np.isnan(ew)):
        #    print("Read {} EWs out of {} lines".format(np.sum(~np.isnan(ew)),len(ew)))
        colnames = colnames + ['equivalent_width']
        dtypes = dtypes + [float]
        data = data + [ew]
        
        return cls(Table(data,names=colnames,dtype=dtypes),moog_columns=moog_columns,**kwargs)

    @classmethod
    def read_GES(cls,filename,moog_columns=False,**kwargs):
        if moog_columns:
            colnames = cls.moog_colnames
            dtypes = cls.moog_dtypes
        else:
            colnames = cls.full_colnames
            dtypes = cls.full_dtypes

        tab = Table.read(filename)
        wl = tab['LAMBDA']

        elem1 = [x.strip() for x in tab['NAME'][:,0]]
        elem2 = [x.strip() for x in tab['NAME'][:,1]]
        numelems = [2 if x=='' else 1 for x in elem2]
        isotope1 = tab['ISOTOPE'][:,0]
        isotope2 = tab['ISOTOPE'][:,1]
        ion = tab['ION']

        # it turns out to be super slow to loop through all the rows; memoizing helps a bit
        memo = {}
        def _get_key(row): 
            return "{}_{}_{}_{}_{}".format(row['NAME'][0],row['NAME'][1],row['ISOTOPE'][0],row['ISOTOPE'][1],row['ION'])
        def _get_species(row):
            key = _get_key(row)
            if key in memo: return memo[key]
            species = elems_isotopes_ion_to_species(row['NAME'][0],row['NAME'][1],row['ISOTOPE'][0],row['ISOTOPE'][1],row['ION'])
            memo[key] = species
            return species
        import time
        start = time.time()
        species = list(map(_get_species, tab))
        print('{:.1f}s to compute species'.format(time.time()-start))

        memo = {}
        def _get_element(species):
            if species in memo: return memo[species]
            element = species_to_element(species)
            memo[species] = element
            return element
        elements = list(map(_get_element, species))

        expot = tab['E_LOW']
        loggf = tab['LOG_GF']
        damp_vdw = tab['VDW_DAMP']
        dissoc_E = np.zeros_like(wl)*np.nan #TODO
        E_hi = tab['E_UP']
        lande_hi = tab['LANDE_UP']
        lande_lo = tab['LANDE_LOW']
        damp_stark = tab['STARK_DAMP']
        damp_rad = tab['RAD_DAMP']
        
        # Takes too long to do this loop; will fill in later
        #refstr = "WL:{},GF:{},EL:{},EU:{},LA:{},RD:{},SD:{},VD:{}"
        #ref_concatenator = lambda row: refstr.format(row['LAMBDA_REF'],row['LOG_GF_REF'],row['E_LOW_REF'],row['E_UP_REF'],
        #                                             row['LANDE_REF'],row['RAD_DAMP_REF'],row['STARK_DAMP_REF'],row['VDW_DAMP_REF'])
        #comments = map(ref_concatenator, tab)
        comments = ['' for row in tab]

        refs = tab['LOG_GF_REF']

        if moog_columns:
            data = [wl,species,expot,loggf,damp_vdw,dissoc_E,comments,
                    numelems,elem1,isotope1,elem2,isotope2,ion,
                    refs,elements]
        else:
            data = [wl,species,expot,loggf,damp_vdw,dissoc_E,comments,
                    numelems,elem1,isotope1,elem2,isotope2,ion,
                    E_hi,lande_hi,lande_lo,damp_stark,damp_rad,refs,elements]
        print('Constructing line list')
        return cls(Table(data,names=colnames,dtype=dtypes),moog_columns=moog_columns,**kwargs)

    def write_moog(self,filename):
        fmt = "{:10.3f}{:10.5f}{:10.3f}{:10.3f}{}{}{}{}"
        space = " "*10
        with open(filename,'w') as f:
            f.write("\n")
            for line in self:
                C6 = space if np.ma.is_masked(line['damp_vdw']) or np.isnan(line['damp_vdw']) else "{:10.3f}".format(line['damp_vdw'])
                D0 = space if np.ma.is_masked(line['dissoc_E']) or np.isnan(line['dissoc_E']) else "{:10.3f}".format(line['dissoc_E'])
                comments = '' if np.ma.is_masked(line['comments']) else line['comments']
                if 'equivalent_width' in line.colnames:
                    EW = space if np.ma.is_masked(line['equivalent_width']) or np.isnan(line['equivalent_width']) else "{:10.3f}".format(line['equivalent_width'])
                else:
                    EW = space
                f.write(fmt.format(line['wavelength'],line['species'],line['expot'],line['loggf'],C6,D0,EW,line['comments'])+"\n")

    def write_latex(self,filename,sortby=['species','wavelength'],
                    write_cols = ['wavelength','element','expot','loggf']):
        new_table = self.copy()
        new_table.sort(sortby)
        new_table = new_table[write_cols]
        new_table.write(filename,format='ascii.aastex')

## Add to astropy.io registry
def _moog_identifier(*args, **kwargs):
    return isinstance(args[0], basestring) and args[0].lower().endswith(".moog")
registry.register_writer("moog", LineList, LineList.write_moog)
registry.register_reader("moog", LineList, LineList.read_moog)
registry.register_identifier("moog", LineList, _moog_identifier)

