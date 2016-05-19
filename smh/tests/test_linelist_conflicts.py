from smh import linelists
import numpy as np

def test_conflicts():
    ll1 = linelists.LineList.read('test_data/linelists/complete.list')
    ll2 = linelists.LineList.read('test_data/linelists/tiII.moog')
    ll3 = linelists.LineList.read('test_data/linelists/lin4554new')
    
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
    
if __name__=="__main__":
    test_conflicts()
