from smh.session import Session

"""
1. Create a new session from multi-spec files.
2. Measure the radial velocity from an un-normalsied order.
3. Apply that radial velocity measurement to the session.
"""

a = Session([
    "/Users/arc/codes/smh/hd44007red_multi.fits",
    "/Users/arc/codes/smh/hd44007blue_multi.fits",
    ])

a.rv_measure("hd140283.fits")