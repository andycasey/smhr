
import numpy as np

from smh.radiative_transfer import (moog, turbospectrum as ts)
from smh.photospheres.marcs import Interpolator
from smh.linelists import LineList

from time import time


marcs = Interpolator()
solar_photosphere = marcs(5777, 4.4, 0)
linelist = LineList.read("atomic.list")

# Just restrict to a few lines.
linelist = linelist[500:505]

t = [time()]
ts_result = ts.synthesize(solar_photosphere, linelist)

t.append(time())
moog_result = moog.synthesize(solar_photosphere, linelist)
t.append(time())

t_ts, t_moog = np.diff(t)

ts_dispersion, ts_flux = ts_result[0][:2]
moog_dispersion, moog_flux = moog_result[0][:2]


import matplotlib.pyplot as plt
fig, ax = plt.subplots()

ax.plot(ts_dispersion, ts_flux)
ax.plot(moog_dispersion, moog_flux)

linelist["equivalent_width"] = np.array([
    60, 0.1, 260, 90, 150])

result = ts.abundance_cog(solar_photosphere, linelist)