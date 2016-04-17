#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" A common application interface for radiative transfer codes. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)


"""
Standardized object for radiative transfer:

def synthesize(photosphere, transitions, **kwargs)
def abundance_cog(photosphere, transitions, **kwargs)

Maybe??:
with smh.rt.instance() as rt:
    rt.synthesize()

"""


from .moog import *