from smh import linelists
import numpy as np

import os
datadir = os.path.dirname(os.path.abspath(__file__))+'/test_data'

def test_conflicts():
    ll1 = linelists.LineList.read(datadir+'/linelists/complete.list')
    ll2 = linelists.LineList.read(datadir+'/linelists/tiII.moog')
    ll3 = linelists.LineList.read(datadir+'/linelists/lin4554new')
    
    # Simple case: 1-1 conflicts
    try:
        ll = ll1.merge(ll2,in_place=False)
    except linelists.LineListConflict as e:
        c1 = e.conflicts1
        c2 = e.conflicts2
        for x,y in zip(c1,c2):
            assert len(x)==len(y)==1
    else:
        raise RuntimeError
    
    # 1-many conflicts
    try:
        ll = ll3.merge(ll1,in_place=False)
    except linelists.LineListConflict as e:
        c1 = e.conflicts1
        c2 = e.conflicts2
        for x in c1:
            if np.all(x['element']=='Ba II'):
                assert len(x)==15
    else:
        raise RuntimeError
    
    # Merging multiple conflicts
    try:
        ll = ll1.merge(ll3,in_place=False,thresh=.1)
    except linelists.LineListConflict as e:
        c1 = e.conflicts1
        c2 = e.conflicts2
        for y in c2:
            if np.all(y['element']=='Ba II'):
                assert len(y)==15
    else:
        raise RuntimeError

if __name__=="__main__":
    ll1 = linelists.LineList.read('test_data/linelists/complete.list')
    ll2 = linelists.LineList.read('test_data/linelists/tiII.moog')
    ll3 = linelists.LineList.read('test_data/linelists/lin4554new')
    
    #c1, c2 = linelists.identify_conflicts(ll1,ll2)
    #c1, c2 = linelists.identify_conflicts(ll1,ll3)
    #c1, c2 = linelists.identify_conflicts(ll3,ll1)
    test_conflicts()
    
    ll2_copy = linelists.LineList.read('test_data/linelists/tiII.moog')
    c1, c2 = linelists.LineList.identify_conflicts(ll2, ll2_copy)
    assert len(c1) == len(ll2)
    
    c1, c2 = linelists.LineList.identify_conflicts(ll2, ll2_copy, skip_equal_loggf=True)
    assert len(c1) == 0
