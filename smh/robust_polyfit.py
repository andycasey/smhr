# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

"""
Small collection of robust statistical estimators based on functions from
Henry Freudenriech (Hughes STX) statistics library (called ROBLIB) that have
been incorporated into the AstroIDL User's Library.  Function included are:
  * biweightMean - biweighted mean estimator
  * mean - robust estimator of the mean of a data set
  * mode - robust estimate of the mode of a data set using the half-sample
    method
  * std - robust estimator of the standard deviation of a data set
  * checkfit - return the standard deviation and biweights for a fit in order 
    to determine its quality
  * linefit - outlier resistant fit of a line to data
  * polyfit - outlier resistant fit of a polynomial to data

For the fitting routines, the coefficients are returned in the same order as
numpy.polyfit, i.e., with the coefficient of the highest power listed first.

For additional information about the original IDL routines, see:
  http://idlastro.gsfc.nasa.gov/contents.html#C17

This code was taken by Alex Ji from Terese Hansen, who got it from Josh Simon.
"""

import math
import numpy

__version__ = '0.4'
__revision__ = '$Rev$'
__all__ = ['biweightMean', 'mean', 'mode', 'std', 'checkfit', 'linefit', 'polyfit', '__version__', '__revision__', '__all__']

__iterMax = 25
__delta = 5.0e-7
__epsilon = 1.0e-20

#print("Note that for the outlier rejection, BisquareLimit=3.0 is used")

def biweightMean(inputData, axis=None, dtype=None):
	"""
	Calculate the mean of a data set using bisquare weighting.  
	
	Based on the biweight_mean routine from the AstroIDL User's 
	Library.
	
	.. versionchanged:: 1.0.3
		Added the 'axis' and 'dtype' keywords to make this function more
		compatible with numpy.mean()
	"""
	
	if axis is not None:
		fnc = lambda x: biweightMean(x, dtype=dtype)
		y0 = numpy.apply_along_axis(fnc, axis, inputData)
	else:
		y = inputData.ravel()
		if type(y).__name__ == "MaskedArray":
			y = y.compressed()
		if dtype is not None:
			y = y.astype(dtype)
			
		n = len(y)
		closeEnough = 0.03*numpy.sqrt(0.5/(n-1))
		
		diff = 1.0e30
		nIter = 0
		
		y0 = numpy.median(y)
		deviation = y - y0
		sigma = std(deviation)
		
		if sigma < __epsilon:
			diff = 0
		while diff > closeEnough:
			nIter = nIter + 1
			if nIter > __iterMax:
				break
			uu = ((y-y0)/(6.0*sigma))**2.0
			uu = numpy.where(uu > 1.0, 1.0, uu)
			weights = (1.0-uu)**2.0
			weights /= weights.sum()
			y0 = (weights*y).sum()
			deviation = y - y0
			prevSigma = sigma
			sigma = std(deviation, Zero=True)
			if sigma > __epsilon:
				diff = numpy.abs(prevSigma - sigma) / prevSigma
			else:
				diff = 0.0
				
	return y0


def mean(inputData, Cut=3.0, axis=None, dtype=None):
	"""
	Robust estimator of the mean of a data set.  Based on the 
	resistant_mean function from the AstroIDL User's Library.
	
	.. versionchanged:: 1.0.3
		Added the 'axis' and 'dtype' keywords to make this function more
		compatible with numpy.mean()
	"""
	
	if axis is not None:
		fnc = lambda x: mean(x, dtype=dtype)
		dataMean = numpy.apply_along_axis(fnc, axis, inputData)
	else:
		data = inputData.ravel()
		if type(data).__name__ == "MaskedArray":
			data = data.compressed()
		if dtype is not None:
			data = data.astype(dtype)
			
		data0 = numpy.median(data)
		maxAbsDev = numpy.median(numpy.abs(data-data0)) / 0.6745
		if maxAbsDev < __epsilon:
			maxAbsDev = (numpy.abs(data-data0)).mean() / 0.8000
			
		cutOff = Cut*maxAbsDev
		good = numpy.where( numpy.abs(data-data0) <= cutOff )
		good = good[0]
		dataMean = data[good].mean()
		dataSigma = math.sqrt( ((data[good]-dataMean)**2.0).sum() / len(good) )

		if Cut > 1.0:
			sigmaCut = Cut
		else:
			sigmaCut = 1.0
		if sigmaCut <= 4.5:
			dataSigma = dataSigma / (-0.15405 + 0.90723*sigmaCut - 0.23584*sigmaCut**2.0 + 0.020142*sigmaCut**3.0)
			
		cutOff = Cut*dataSigma
		good = numpy.where(  numpy.abs(data-data0) <= cutOff )
		good = good[0]
		dataMean = data[good].mean()
		if len(good) > 3:
			dataSigma = math.sqrt( ((data[good]-dataMean)**2.0).sum() / len(good) )
			
		if Cut > 1.0:
			sigmaCut = Cut
		else:
			sigmaCut = 1.0
		if sigmaCut <= 4.5:
			dataSigma = dataSigma / (-0.15405 + 0.90723*sigmaCut - 0.23584*sigmaCut**2.0 + 0.020142*sigmaCut**3.0)
			
		dataSigma = dataSigma / math.sqrt(len(good)-1)
		
	return dataMean


def mode(inputData, axis=None, dtype=None):
	"""
	Robust estimator of the mode of a data set using the half-sample mode.
	
	.. versionadded: 1.0.3
	"""
	
	if axis is not None:
		fnc = lambda x: mode(x, dtype=dtype)
		dataMode = numpy.apply_along_axis(fnc, axis, inputData)
	else:
		# Create the function that we can use for the half-sample mode
		def _hsm(data):
			if data.size == 1:
				return data[0]
			elif data.size == 2:
				return data.mean()
			elif data.size == 3:
				i1 = data[1] - data[0]
				i2 = data[2] - data[1]
				if i1 < i2:
					return data[:2].mean()
				elif i2 > i1:
					return data[1:].mean()
				else:
					return data[1]
			else:
				wMin = data[-1] - data[0]
				N = data.size/2 + data.size%2 
				for i in xrange(0, N):
					w = data[i+N-1] - data[i] 
					if w < wMin:
						wMin = w
						j = i
				return _hsm(data[j:j+N])
				
		data = inputData.ravel()
		if type(data).__name__ == "MaskedArray":
			data = data.compressed()
		if dtype is not None:
			data = data.astype(dtype)
			
		# The data need to be sorted for this to work
		data = numpy.sort(data)
		
		# Find the mode
		dataMode = _hsm(data)
		
	return dataMode


def std(inputData, Zero=False, axis=None, dtype=None):
	"""
	Robust estimator of the standard deviation of a data set.  
	
	Based on the robust_sigma function from the AstroIDL User's Library.
	
	.. versionchanged:: 1.0.3
		Added the 'axis' and 'dtype' keywords to make this function more
		compatible with numpy.std()
	"""
	
	if axis is not None:
		fnc = lambda x: std(x, dtype=dtype)
		sigma = numpy.apply_along_axis(fnc, axis, inputData)
	else:
		data = inputData.ravel()
		if type(data).__name__ == "MaskedArray":
			data = data.compressed()
		if dtype is not None:
			data = data.astype(dtype)
			
		if Zero:
			data0 = 0.0
		else:
			data0 = numpy.median(data)
		maxAbsDev = numpy.median(numpy.abs(data-data0)) / 0.6745
		if maxAbsDev < __epsilon:
			maxAbsDev = (numpy.abs(data-data0)).mean() / 0.8000
		if maxAbsDev < __epsilon:
			sigma = 0.0
			return sigma
			
		u = (data-data0) / 6.0 / maxAbsDev
		u2 = u**2.0
		good = numpy.where( u2 <= 1.0 )
		good = good[0]
		if len(good) < 3:
			print("WARNING:  Distribution is too strange to compute standard deviation")
			sigma = -1.0
			return sigma
			
		numerator = ((data[good]-data0)**2.0 * (1.0-u2[good])**2.0).sum()
		nElements = (data.ravel()).shape[0]
		denominator = ((1.0-u2[good])*(1.0-5.0*u2[good])).sum()
		sigma = nElements*numerator / (denominator*(denominator-1.0))
		if sigma > 0:
			sigma = math.sqrt(sigma)
		else:
			sigma = 0.0
			
	return sigma


def checkfit(inputData, inputFit, epsilon, delta, BisquareLimit=3.0):
	"""
	Determine the quality of a fit and biweights.  Returns a tuple
	with elements:
	  0. Robust standard deviation analog
	  1. Fractional median absolute deviation of the residuals
	  2. Number of input points given non-zero weight in the calculation
	  3. Bisquare weights of the input points
	  4. Residual values scaled by sigma
	
	This function is based on the rob_checkfit routine from the AstroIDL 
	User's Library.
	"""
	
	data = inputData.ravel()
	fit = inputFit.ravel()
	if type(data).__name__ == "MaskedArray":
		data = data.compressed()
	if type(fit).__name__ == "MaskedArray":
		fit = fit.compressed()

	deviation = data - fit
	sigma = std(deviation, Zero=True)
	if sigma < epsilon:
		return (sigma, 0.0, 0, 0.0, 0.0)
	
	toUse = (numpy.where( numpy.abs(fit) > epsilon ))[0]
	if len(toUse) < 3:
		fracDev = 0.0
	else:
		fracDev = numpy.median(numpy.abs(deviation[toUse]/fit[toUse]))
	if fracDev < delta:
		return (sigma, fracDev, 0, 0.0, 0.0)
		
	biweights = numpy.abs(deviation)/(BisquareLimit*sigma)
	toUse = (numpy.where(biweights > 1))[0]
	if len(toUse) > 0:
		biweights[toUse] = 1.0
	nGood = len(data) - len(toUse)
	
	scaledResids = (1.0 - biweights**2.0)
	scaledResids = scaledResids / scaledResids.sum()
	
	return (sigma, fracDev, nGood, biweights, scaledResids)


def linefit(inputX, inputY, iterMax=25, Bisector=False, BisquareLimit=6.0, CloseFactor=0.03):
	"""
	Outlier resistance two-variable linear regression function.
	
	Based on the robust_linefit routine in the AstroIDL User's Library.
	"""
	
	xIn = inputX.ravel()
	yIn = inputY.ravel()
	if type(yIn).__name__ == "MaskedArray":
		xIn = xIn.compress(numpy.logical_not(yIn.mask))
		yIn = yIn.compressed()
	n = len(xIn)
	
	x0 = xIn.sum() / n
	y0 = yIn.sum() / n
	x = xIn - x0
	y = yIn - y0
	
	cc = numpy.zeros(2)
	ss = numpy.zeros(2)
	sigma = 0.0
	yFit = yIn
	badFit = 0
	nGood = n
	
	lsq = 0.0
	yp = y
	if n > 5:
		s = numpy.argsort(x)
		u = x[s]
		v = y[s]
		nHalf = n/2 -1
		x1 = numpy.median(u[0:nHalf])
		x2 = numpy.median(u[nHalf:])
		y1 = numpy.median(v[0:nHalf])
		y2 = numpy.median(v[nHalf:])
		if numpy.abs(x2-x1) < __epsilon:
			x1 = u[0]
			x2 = u[-1]
			y1 = v[0]
			y2 = v[-1]
		cc[1] = (y2-y1)/(x2-x1)
		cc[0] = y1 - cc[1]*x1
		yFit = cc[0] + cc[1]*x
		sigma, fracDev, nGood, biweights, scaledResids = checkfit(yp, yFit, __epsilon, __delta)
		if nGood < 2:
			lsq = 1.0
			
	if lsq == 1 or n < 6:
		sx = x.sum()
		sy = y.sum()
		sxy = (x*y).sum()
		sxx = (x*x).sum()
		d = sxx - sx*sx
		if numpy.abs(d) < __epsilon:
			return (0.0, 0.0)
		ySlope = (sxy - sx*sy) / d
		yYInt = (sxx*sy - sx*sxy) / d
		
		if Bisector:
			syy = (y*y).sum()
			d = syy - sy*sy
			if numpy.abs(d) < __epsilon:
				return (0.0, 0.0)
			tSlope = (sxy - sy*sx) / d
			tYInt = (syy*sx - sy*sxy) / d
			if numpy.abs(tSlope) < __epsilon:
				return (0.0, 0.0)
			xSlope = 1.0/tSlope
			xYInt = -tYInt / tSlope
			if ySlope > xSlope:
				a1 = yYInt
				b1 = ySlope
				r1 = numpy.sqrt(1.0+ySlope**2.0)
				a2 = xYInt
				b2 = xSlope
				r2 = numpy.sqrt(1.0+xSlope**2.0)
			else:
				a2 = yYInt
				b2 = ySlope
				r2 = numpy.sqrt(1.0+ySlope**2.0)
				a1 = xYInt
				b1 = xSlope
				r1 = numpy.sqrt(1.0+xSlope**2.0)
			yInt = (r1*a2 + r2*a1) / (r1 + r2)
			slope = (r1*b2 + r2*b1) / (r1 + r2)
			r = numpy.sqrt(1.0+slope**2.0)
			if yInt > 0:
				r = -r
			u1 = slope / r
			u2 = -1.0/r
			u3 = yInt / r
			yp = u1*x + u2*y + u3
			yFit = y*0.0
			ss = yp
		else:
			slope = ySlope
			yInt = yYInt
			yFit = yInt + slope*x
		cc[0] = yInt
		cc[1] = slope
		sigma, fracDev, nGood, biweights, scaledResids = checkfit(yp, yFit, __epsilon, __delta)
		
	if nGood < 2:
		cc[0] = cc[0] + y0 - cc[1]*x0
		return cc[::-1]
		
	sigma1 = (100.0*sigma)
	closeEnough = CloseFactor * numpy.sqrt(0.5/(n-1))
	if closeEnough < __delta:
		closeEnough = __delta
	diff = 1.0e20
	nIter = 0
	while diff > closeEnough:
		nIter = nIter + 1
		if nIter > iterMax:
			break
		sigma2 = sigma1
		sigma1 = sigma
		sx = (biweights*x).sum()
		sy = (biweights*y).sum()
		sxy = (biweights*x*y).sum()
		sxx = (biweights*x*x).sum()
		d = sxx - sx*sx
		if numpy.abs(d) < __epsilon:
			return (0.0, 0.0)
		ySlope = (sxy - sx*sy) / d
		yYInt = (sxx*sy - sx*sxy) / d
		slope = ySlope
		yInt = yYInt
		
		if Bisector:
			syy = (biweights*y*y).sum()
			d = syy - sy*sy
			if numpy.abs(d) < __epsilon:
				return (0.0, 0.0)
			tSlope = (sxy - sy*sx) / d
			tYInt = (syy*sx - sy*sxy) / d
			if numpy.abs(tSlope) < __epsilon:
				return (0.0, 0.0)
			xSlope = 1.0/tSlope
			xYInt = -tYInt / tSlope
			if ySlope > xSlope:
				a1 = yYInt
				b1 = ySlope
				r1 = numpy.sqrt(1.0+ySlope**2.0)
				a2 = xYInt
				b2 = xSlope
				r2 = numpy.sqrt(1.0+xSlope**2.0)
			else:
				a2 = yYInt
				b2 = ySlope
				r2 = numpy.sqrt(1.0+ySlope**2.0)
				a1 = xYInt
				b1 = xSlope
				r1 = numpy.sqrt(1.0+xSlope**2.0)
			yInt = (r1*a2 + r2*a1) / (r1 + r2)
			slope = (r1*b2 + r2*b1) / (r1 + r2)
			r = numpy.sqrt(1.0+slope**2.0)
			if yInt > 0:
				r = -r
			u1 = slope / r
			u2 = -1.0/r
			u3 = yInt / r
			yp = u1*x + u2*y + u3
			yFit = y*0.0
			ss = yp
		else:
			yFit = yInt + slope*x
		cc[0] = yInt
		cc[1] = slope
		sigma, fracDev, nGood, biweights, scaledResids = checkfit(yp, yFit, __epsilon, __delta)
		
		if nGood < 2:
			badFit = 1
			break
		diff1 = numpy.abs(sigma1 - sigma)/sigma
		diff2 = numpy.abs(sigma2 - sigma)/sigma
		if diff1 < diff2:
			diff = diff1
		else:
			diff = diff2
			
	cc[0] = cc[0] + y0 - cc[1]*x0
	return cc[::-1]


def polyfit(inputX, inputY, order, iterMax=25):
	"""
	Outlier resistance two-variable polynomial function fitter.
	
	Based on the robust_poly_fit routine in the AstroIDL User's 
	Library.
	
	Unlike robust_poly_fit, two different polynomial fitters are used
	because numpy.polyfit does not support non-uniform weighting of the
	data.  For the weighted fitting, the SciPy Orthogonal Distance
	Regression module (scipy.odr) is used.
	"""
	
	from scipy import odr
	
	def polyFunc(B, x, order=order):
		out = x*0.0
		for i in range(order+1):
			out = out + B[i]*x**i
	
	model = odr.Model(polyFunc)
	
	x = inputX.ravel()
	y = inputY.ravel()
	if type(y).__name__ == "MaskedArray":
		x = x.compress(numpy.logical_not(y.mask))
		y = y.compressed()
	n = len(x)
	
	x0 = x.sum() / n
	y0 = y.sum() / n
	u = x
	v = y
	
	nSeg = order + 2
	if (nSeg//2)*2 == nSeg:
		nSeg = nSeg + 1
	minPts = nSeg*3
	if n < 1000:
		lsqFit = 1
		cc = numpy.polyfit(u, v, order)
		yFit = numpy.polyval(cc, u)
	else:
		lsqfit = 0
		q = numpy.argsort(u)
		u = u[q]
		v = v[q]
		nPerSeg = numpy.zeros(nSeg, dtype=int) + n//nSeg
		nLeft = n - nPerSeg[0]*nSeg
		nPerSeg[nSeg//2] = nPerSeg[nSeg//2] + nLeft
		r = numpy.zeros(nSeg)
		s = numpy.zeros(nSeg)
		r[0] = numpy.median(u[0:nPerSeg[0]])
		s[0] = numpy.median(v[0:nPerSeg[0]])
		i2 = nPerSeg[0]-1
		for i in range(1,nSeg):
			i1 = i2
			i2 = i1 + nPerSeg[i]
			r[i] = numpy.median(u[i1:i2])
			s[i] = numpy.median(v[i1:i2])
		cc = numpy.polyfit(r, s, order)
		yFit = numpy.polyval(cc, u)
		
	sigma, fracDev, nGood, biweights, scaledResids = checkfit(v, yFit, __epsilon, __delta)
	if nGood == 0:
		return cc, np.nan
	if nGood < minPts:
		if lsqFit == 0:
			cc = numpy.polyfit(u, v, order)
			yFit = numpy.polyval(cc, u)
			sigma, fracDev, nGood, biweights, scaledResids = checkfit(yp, yFit, __epsilon, __delta)
			if nGood == 0:
				return __processPoly(x0, y0, order, cc)
			nGood = n - nGood
		if nGood < minPts:
			return 0, np.nan
			
	closeEnough = 0.03*numpy.sqrt(0.5/(n-1))
	if closeEnough < __delta:
		closeEnough = __delta
	diff = 1.0e10
	sigma1 = 100.0*sigma
	nIter = 0
	while diff > closeEnough:
		nIter = nIter + 1
		if nIter > iterMax:
			break
		sigma2 = sigma1
		sigma1 = sigma
		g = (numpy.where(biweights < 1))[0]
		if len(g) < len(biweights):
			u = u[g]
			v = v[g]
			biweights = biweights[g]
		try:
			## Try the fancy method...
			data = odr.RealData(u, v, sy=1.0/biweights)
			fit = odr.ODR(data, model, beta0=cc[::-1])
			out = fit.run()
			cc = out.beta[::-1]
		except:
			## And then give up when it doesn't work
			cc = numpy.polyfit(u, v, order)
		yFit = numpy.polyval(cc, u)
		sigma, fracDev, nGood, biweights, scaledResids = checkfit(v, yFit, __epsilon, __delta)
		if nGood < minPts:
			return cc, np.nan
		diff1 = numpy.abs(sigma1 - sigma)/sigma
		diff2 = numpy.abs(sigma2 - sigma)/sigma
		if diff1 < diff2:
			diff = diff1
		else:
			diff = diff2
	return cc, sigma

import numpy as np
from scipy import optimize
def gaussian(x,A,x0,err,B):
	return A * np.exp(-(x-x0)**2/(2.*err**2)) + B
def fit_gaussian(x,y,p0=None,yerr=None, **kwargs):
	assert np.all(np.isfinite(x)) & np.all(np.isfinite(y))
	if p0 is None:
		p0 = [np.max(y), (np.max(x)-np.min(x))/2., np.median(x), np.min(y)]
		popt, pcov = optimize.curve_fit(gaussian, x, y, p0=p0,
						bounds=([0,           np.min(x), 0,                       0],
							[2*np.max(y), np.max(x), 3*(np.max(x)-np.min(x)), np.max(y)]),
						sigma=yerr,
						**kwargs)
	return popt, pcov

def gfunc3(x, *theta):
	z = (x-theta[1])/theta[2]
	return theta[0] * np.exp(-z**2/2.)
def gfunc4(x, *theta):
	z = (x-theta[1])/theta[2]
	return theta[0] * np.exp(-z**2/2.) + theta[3]
def gfunc5(x, *theta):
	z = (x-theta[1])/theta[2]
	return theta[0] * np.exp(-z**2/2.) + theta[3] + theta[4]*x
def gfunc6(x, *theta):
	z = (x-theta[1])/theta[2]
	return theta[0] * np.exp(-z**2/2.) + theta[3] + theta[4]*x + theta[5]*x**2
def gaussfit(xdata, ydata, p0, **kwargs):
	"""
	p0 = (amplitude, mean, sigma) (bias; linear; quadratic)
	"""
	NTERMS = len(p0)
	if NTERMS == 3:
		func = gfunc3
	elif NTERMS == 4:
		func = gfunc4
	elif NTERMS == 5:
		func = gfunc5
	elif NTERMS == 6:
		func = gfunc6
	else:
		raise ValueError("p0 must be 3-6 terms long, {}".format(p0))
        
	popt, pcov = optimize.curve_fit(func, xdata, ydata, p0, **kwargs)
	return popt
