from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import numpy as np
from astropy.io import ascii,fits
from astropy.table import Table, Column, MaskedColumn

def line_in_list(line,ll,thresh=.01,verbose=False):
    ii1 = ll['species']==line['species']
    ii2 = np.abs(ll['wavelength']-line['wavelength']) < thresh
    ii = np.logical_and(ii1,ii2)
    num_match = np.sum(ii)
    if num_match==0: return -1
    if num_match==1: return np.where(ii)[0]
    if verbose:
        print("Error: {} matches!".format(num_match))
        print(line)
    return -1 * num_match

def read_moog_linelist(filename,full_columns=True):
    if full_columns:
        colnames = ['wavelength','species','expot','loggf','damp_vdw','dissoc_E','comments',
                    'E_hi','E_lo','lande_hi','lande_lo','damp_stark','damp_rad','references']
        dtypes = [np.float,np.float,np.float,np.float,np.float,np.float,str,
                  np.float,np.float,np.float,np.float,np.float,np.float,str]
    else:
        colnames = ['wavelength','species','expot','loggf','damp_vdw','dissoc_E','comments','references']
        dtypes = [np.float,np.float,np.float,np.float,np.float,np.float,str,str]

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

    return Table(data,names=colnames,dtype=dtypes)

def write_moog_linelist(ll,filename):
    fmt = "{:10.3f}{:10.5f}{:10.3f}{:10.3f}{}{}{}{}"
    space = " "*10
    with open(filename,'w') as f:
        for line in ll:
            C6 = space if np.ma.is_masked(line['damp_vdw']) else "{:10.3}".format(line['damp_vdw'])
            D0 = space if np.ma.is_masked(line['dissoc_E']) else "{:10.3}".format(line['dissoc_E'])
            comments = '' if np.ma.is_masked(line['comments']) else line['comments']
            f.write(fmt.format(line['wavelength'],line['species'],line['expot'],line['loggf'],C6,D0,space,line['comments'])+"\n")
