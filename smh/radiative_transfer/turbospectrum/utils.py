

import logging
import os
from smh.utils import mkdtemp

logger = logging.getLogger(__name__)

def twd_path(twd=None,**kwargs):
    """
    Create a temporary working directory and return a function that will format
    basenames from that temporary working directory.
    """

    if twd is None:
        kwds = {}
        kwds["dir"] = kwargs.get("dir", "/tmp/")
        kwds["prefix"] = kwargs.get("prefix", "smh-")
        kwds["suffix"] = kwargs.get("suffix", "")
        #kwds = kwargs.copy()
        #kwds.setdefault("dir", "/tmp/")
        #kwds.setdefault("prefix", "smh-")
        twd = mkdtemp(**kwds)
    return lambda filename: os.path.join(twd, filename)