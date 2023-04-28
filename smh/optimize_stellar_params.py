#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import numpy as np
from scipy.optimize import fsolve
from . import (photospheres, radiative_transfer, utils)
from functools import wraps
import sys, os, time
import atexit
from shutil import rmtree
import warnings
from copy import deepcopy

#from astropy.stats.biweight import biweight_scale
from .photospheres.abundances import asplund_2009 as solar_composition
from .utils import mkdtemp

import logging
logger = logging.getLogger(__name__)


def check_satisfies_tolerance(final_parameters_result, total_tolerance, individual_tolerances):
    acquired_total_tolerance = np.sum(np.array(final_parameters_result)**2)
    return acquired_total_tolerance < total_tolerance and \
        (individual_tolerances is None or np.all(np.less_equal(np.abs(final_parameters_result), individual_tolerances)))
    
class OptimizationSuccess(BaseException):
    pass


def stellar_optimization(func):
    """ A decorator to wrap the stellar optimization function such
    that it finishes with the parameter tolerance we require and
    doesn't make unncessary (i.e., repeated) calls to MOOG """

    def _decorator(request, *args, **kwargs):

        previously_sampled_points, total_tolerance, individual_tolerances, use_nlte_grid = args
        
        previously_sampled_points = np.array(previously_sampled_points)

        # Checking if this point is sampled already
        if len(previously_sampled_points) > 0:

            this_point = np.array(request)
            #import pdb; pdb.set_trace()

            try:
                sampled_before = (previously_sampled_points[:, :4] == this_point).all(axis=1)
            except:
                sampled_before = (previously_sampled_points[:, :4] == this_point).all()
                #np.arange(previously_sampled_points[:, :4].ndim - this_point.ndim, previously_sampled_points[:, :4].ndim))
            #print(sampled_before)

            if np.any(sampled_before):
                index = np.where(sampled_before)[0][0]
                logger.debug("Sampled these stellar parameters already, so we're just returning those values.")
                logger.debug(previously_sampled_points[index])

                return previously_sampled_points[index, 4:]

        # Perform the optimization
        response = func(request, *args, **kwargs)

        # Check for accepted tolerances, both total tolerance and individual tolerances if they are specified
        acquired_total_tolerance = np.sum(pow(response, 2))
        if total_tolerance >= acquired_total_tolerance and \
            (individual_tolerances is None or np.all(np.less_equal(np.abs(response), individual_tolerances))):

            message = "Optimization complete. Total tolerance of <{0:.1e} met: {1:.1e}".format(
                total_tolerance, acquired_total_tolerance)

            if individual_tolerances is not None:
                message += ", and individual tolerances of <[{0:.1e}, {1:.1e}, {2:.1e}, {3:.1e}] met: [{4:.1e}, {5:.1e}, {6:.1e}, {7:.1e}]".format(
                    *tuple(list(response) + list(individual_tolerances)))

            else:
                message += ", and no individual tolerances specified."

            raise OptimizationSuccess(message)


        return response
    return wraps(func)(_decorator)


# E. Holmbeck just copied the above decorator, but changed 4-->2.
# It might be wrong...
def feh_optimization(func):
    """ A decorator to wrap the stellar optimization function such
    that it finishes with the parameter tolerance we require and
    doesn't make unncessary (i.e., repeated) calls to MOOG """

    def _decorator(request, *args, **kwargs):

        params_to_optimize, previously_sampled_points, total_tolerance, individual_tolerances, use_nlte_grid = args

        previously_sampled_points = np.array(previously_sampled_points)
        
        # Checking if this point is sampled already
        if len(previously_sampled_points) > 0:

            this_point = np.array(request)
            #import pdb; pdb.set_trace()

            try:
                sampled_before = (previously_sampled_points[:, :4] == this_point).all(axis=1)
            except:
                sampled_before = (previously_sampled_points[:, :4] == this_point)
                #np.arange(previously_sampled_points[:, :4].ndim - this_point.ndim, previously_sampled_points[:, :4].ndim))
            #print sampled_before

            if np.any(sampled_before):
                index = np.where(sampled_before)[0][0]
                logger.debug("Sampled these stellar parameters already, so we're just returning those values.")
                logger.debug(previously_sampled_points[index])

                return previously_sampled_points[index, 4:]

        # Perform the optimization
        response = func(request, *args, **kwargs)

        # Check for accepted tolerances, both total tolerance and individual tolerances if they are specified
        acquired_total_tolerance = np.sum(pow(response, 2))
        if total_tolerance >= acquired_total_tolerance and \
            (individual_tolerances is None or np.all(np.less_equal(np.abs(response), individual_tolerances))):

            message = "Optimization complete. Total tolerance of <{0:.1e} met: {1:.1e}".format(
                total_tolerance, acquired_total_tolerance)

            if individual_tolerances is not None:
                message += ", and individual tolerances of <[{0:.1e}, {1:.1e}, {2:.1e}, {3:.1e}] met: [{4:.1e}, {5:.1e}, {6:.1e}, {7:.1e}]".format(
                    *tuple(list(response) + list(individual_tolerances)))

            else:
                message += ", and no individual tolerances specified."

            raise OptimizationSuccess(message)


        return response
    return wraps(func)(_decorator)
    


def optimize_stellar_parameters(initial_guess, transitions, EWs=None, 
                                jacobian=utils.approximate_sun_hermes_jacobian,
                                rt=radiative_transfer.moog,
                                max_attempts=5, total_tolerance=1e-4, 
                                individual_tolerances=None, 
                                maxfev=30, use_nlte_grid=None,
                                alpha=0.4):
    """
    Assumes these are all transitions you want to use for stellar parameters
    Assumes you only want to balance neutral ions against Expot and REW
    
    If specify EWs, all EWs are in same order as transitions

    initial_guess order : [teff, vt, logg, feh]
    """
    
    if EWs is None:
        EWs = transitions["equivalent_width"]
        REWs = np.log10(1e-3 * EWs / transitions['wavelength'])
    else:    
        REWs = np.log10(1e-3 * EWs / transitions['wavelength'])
        transitions["equivalent_width"] = EWs
    transitions["reduced_equivalent_width"] = REWs
    idx_I  = transitions["ion"] == 1
    idx_II = transitions["ion"] == 2

    photosphere_interpolator = photospheres.interpolator()
    parameter_ranges = {
        "teff": (3500, 7000),
        "vt": (0.0, 4.0),
        "logg": (0, 5),
        "[Fe/H]": (-5, 0.5)
        }
    
    # Create a mutable copy for the initial guess
    solver_guess = []
    solver_guess.extend(initial_guess)

    # Create a new temporary working directory for cog
    twd = mkdtemp()
    atexit.register(rmtree, twd)

    @stellar_optimization
    def minimisation_function(stellar_parameters, *args):
        """The function we want to minimise (e.g., calculates the quadrature
        sum of all minimizable variables.)
        
        stellar_parameters : [teff, vt, logg, feh]
        """
        ## "Global" vars: transitions, idx_I, idx_II, photosphere_interpolator

        teff, vt, logg, feh = stellar_parameters
        #if teff < parameter_ranges["teff"][0] or teff > parameter_ranges["teff"][1] or \
        #        vt < parameter_ranges["vt"][0] or vt > parameter_ranges["vt"][1] or \
        #        logg < parameter_ranges["logg"][0] or logg > parameter_ranges["logg"][1] or \
        #        feh < parameter_ranges["[Fe/H]"][0] or feh > parameter_ranges["[Fe/H]"][1]:
        #    return np.array([np.nan, np.nan, np.nan, np.nan])

        all_sampled_points, total_tolerance, individual_tolerances, use_nlte_grid = args

        #if not (5 > vt > 0):
        #    return np.array([np.nan, np.nan, np.nan, np.nan])
        
        photosphere = photosphere_interpolator(teff, logg, feh, alpha)
        photosphere.meta["stellar_parameters"]["microturbulence"] = vt
        
        ## TODO: ADJUST ABUNDANCES TO ASPLUND?
        abundances = rt.abundance_cog(photosphere,transitions,twd=twd)
        transitions["abundance"] = abundances
        
        ## Calculate slopes and differences that are being minimized
        out = utils.equilibrium_state(transitions[idx_I],
                                      ("expot", "reduced_equivalent_width"))
        dAdchi = out[26.0]['expot'][0]
        dAdREW = out[26.0]['reduced_equivalent_width'][0]
        dFe = np.mean(abundances[idx_I]) - np.mean(abundances[idx_II])
        dM  = np.mean(abundances[idx_I]) - (feh + solar_composition("Fe"))
        results = np.array([dAdchi, dAdREW, 0.1 * dFe, 0.1 * dM])
        acquired_total_tolerance = np.sum(results**2)

        point = [teff, vt, logg, feh] + list(results)
        all_sampled_points.append(point)

        logger.info("Atmosphere with Teff = {0:.0f} K, vt = {1:.2f} km/s, logg = {2:.2f}, [Fe/H] = {3:.2f}, [alpha/Fe] = {4:.2f}"
                    " yields sum {5:.1e}:\n\t\t\t[{6:.1e}, {7:.1e}, {8:.1e}, {9:.1e}]".format(teff, vt, logg, feh, alpha,
                                                                                              acquired_total_tolerance, *results))
        return results

    all_sampled_points = []
        
    start = time.time()
    for i in range(1, 1 + max_attempts):
        sampled_points = []
        args = (sampled_points, total_tolerance, individual_tolerances, 
                use_nlte_grid)

        #results = fsolve(minimisation_function, solver_guess, args=args, fprime=utils.approximate_stellar_jacobian,
        #                 col_deriv=1, epsfcn=0, xtol=1e-10, full_output=1, maxfev=maxfev)
        try:
            results = fsolve(minimisation_function, solver_guess, args=args, fprime=utils.approximate_stellar_jacobian,
                             col_deriv=1, epsfcn=0, xtol=1e-10, full_output=1, maxfev=maxfev)
        except OptimizationSuccess as e:#: #Exception as e:# OptimizationSuccess as e:
            print(e)
            
            # Optimization is complete and tolerances have been reached
            t_elapsed = time.time() - start
            #logger.info("Successful after {0:.0f} seconds".format(t_elapsed))
    
            point_results = np.sum(np.array(sampled_points)[:, 4:]**2, axis=1)
            min_index = np.argmin(point_results)

            final_parameters = sampled_points[min_index][:4]
            final_parameters[0] = int(np.round(final_parameters[0])) # Effective temperature

            final_parameters_result = sampled_points[min_index][4:]
            num_moog_iterations = len(sampled_points)
            all_sampled_points.extend(sampled_points)

            acquired_total_tolerance = sum(pow(np.array(final_parameters_result), 2))
            return (True, initial_guess, num_moog_iterations, i, t_elapsed, 
                    final_parameters, final_parameters_result, 
                    np.array(all_sampled_points))
        else:
            t_elapsed = time.time() - start
            ier, mesg = results[-2:]
            #if ier != 1:
            #    logger.info("Optimizer has returned with the message: {0}".format(mesg))
            
            final_parameters = results[0]
            final_parameters[0] = int(np.round(final_parameters[0])) # Effective temperature
            final_parameters_result = results[1]["fvec"]
            num_moog_iterations = len(sampled_points)
            all_sampled_points.extend(sampled_points)
            acquired_total_tolerance = sum(pow(np.array(final_parameters_result), 2))

            # Have we achieved tolerance?
            if total_tolerance >= acquired_total_tolerance and \
                    (individual_tolerances is None or np.all(np.less_equal(np.abs(final_parameters_result), individual_tolerances))):
                
                return (True, initial_guess, num_moog_iterations, i, t_elapsed, final_parameters, final_parameters_result, np.array(all_sampled_points))
            else:
                # This solution will fail. Try a new starting point
                solver_guess = [np.random.uniform(*parameter_ranges[parameter]) for parameter in ["teff", "vt", "logg", "[Fe/H]"]]
                
        t_elapsed = time.time() - start

        if total_tolerance >= acquired_total_tolerance and \
            (individual_tolerances is None or np.all(np.less_equal(np.abs(final_parameters_result), individual_tolerances))):
            tolerance_achieved = True
            
        else:
            tolerance_achieved = False

        return (tolerance_achieved, initial_guess, num_moog_iterations, i, t_elapsed, final_parameters, final_parameters_result, np.array(all_sampled_points))

    pass

def optimize_stellar_parameters_2(initial_guess, transitions, EWs=None,
                                  alphafe=0.4,
                                  expot_balance_species=26.0,
                                  ionization_balance_species_1=26.0,
                                  ionization_balance_species_2=26.1,
                                  rew_balance_species=26.0,
                                  mh_species=26.0,
                                  min_eqw=0.01, max_eqw=9999.9,
                                  min_lines=3,
                                  sigma_clip=5,
                                  max_attempts=9, total_tolerance=1e-4, 
                                  individual_tolerances=None, 
                                  maxfev=30,
                                  rt=radiative_transfer.moog,
                                  use_nlte_grid=None,
                                  teff_range=(3500,7000),logg_range=(0.0,5.0),
                                  vt_range=(0.0,4.0),mh_range=(-5.0,0.5),
                                  twd=None):
    """
    An updated and more flexible way to optimize stellar parameters.
    
    initial_guess order : [teff, logg, vt, MH] (DIFFERENT!)
    transitions: the list of lines to use for optimization
    EWs: if not in transitions, the equivalent widths of those lines
    alphafe = 0.4: the fixed [alpha/Fe] for the model atmosphere
    
    expot_balance_species = 26.0: the species used to optimize for Teff. Usually 26.0 or 26.1
    ionization_balance_species_1/2 = 26.0/26.1: the species used to optimize for logg. Usually 26.0/26.1
    rew_balance_species = 26.0: the species used to optimize for vt. Usually 26.0 or 26.1
    mh_species = 26.0: the species used to determine [M/H]
    min_eqw, max_eqw = 0.01, 9999.9: outside of this range, assume the eqw is invalid
    
    min_lines = 3: minimum number of lines for any balance
    sigma_clip = 5: remove outlier abundances, unless you would go below min_lines
    
    max_attempts = 5: maximum number of iterations for optimization
    total_tolerance = 1e-4: the abs squared sum of dA/dchi, dA/dREW, 0.1*(FeI-FeII), 0.1*dMH
    individual_tolerances = None: 
    maxfev = 30:
    
    jacobian = utils.approximate_sun_hermes_jacobian: the 
    rt = radiative_transfer.moog: which radiative transfer code to use (only MOOG right now)
    use_nlte_grid = None: not implemented
    
    teff_range = (3500,7000), logg_range = (0.0, 5.0), vt_range = (0.0, 4.0), mh_range = (-5.0, 0.5)
      maximum limits for optimization
    twd = None: if you want to specify a working directory for the MOOG files, you can do that here
    
    output: 
    tolerance_achieved, initial_guess, num_moog_iterations, num_iters, t_elapsed,
    final_parameters, final_parameters_result, all_sampled_points, valid
    
    """
    
    ## Validate inputs
    if EWs is None:
        EWs = transitions["equivalent_width"]
        REWs = np.log10(1e-3 * EWs / transitions['wavelength'])
    else:    
        REWs = np.log10(1e-3 * EWs / transitions['wavelength'])
        transitions["equivalent_width"] = EWs
    transitions["reduced_equivalent_width"] = REWs
    
    finite = (transitions["equivalent_width"] > min_eqw) & (transitions["equivalent_width"] < max_eqw)
    if np.sum(finite) != len(finite): warnings.warn("Only {}/{} eqw are finite or within {} < eqw < {}".format(
            finite.sum(), len(finite), min_eqw, max_eqw))
    
    idx_expot = (transitions["species"] == expot_balance_species) & finite
    if np.sum(idx_expot) < min_lines:
        raise RuntimeError("Fewer than {} lines for expot species {:.1f}".format(min_lines, expot_balance_species))
    idx_rew = (transitions["species"] == rew_balance_species) & finite
    if np.sum(idx_rew) < min_lines:
        raise RuntimeError("Fewer than {} lines for rew species {:.1f}".format(min_lines, rew_balance_species))
    idx_I = (transitions["species"] == ionization_balance_species_1) & finite
    if np.sum(idx_I) < min_lines:
        raise RuntimeError("Fewer than {} lines for Ion I species {:.1f}".format(min_lines, ionization_balance_species_1))
    idx_II = (transitions["species"] == ionization_balance_species_2) & finite
    if np.sum(idx_II) < min_lines:
        raise RuntimeError("Fewer than {} lines for Ion II species {:.1f}".format(min_lines, ionization_balance_species_2))
    idx_MH = (transitions["species"] == mh_species) & finite
    
    unique_species = np.unique([expot_balance_species, rew_balance_species, ionization_balance_species_1,
                                ionization_balance_species_2, mh_species])
    
    valid = idx_expot | idx_rew | idx_I | idx_II | idx_MH
    if valid.sum() != len(transitions):
        warnings.warn("Not all lines are being used!")
    
    ## Set up needed variables
    photosphere_interpolator = photospheres.interpolator()
    parameter_ranges = deepcopy({
        "teff": teff_range, "logg": logg_range,
        "vt": vt_range, "[Fe/H]": mh_range
    })
    
    # Create a mutable copy for the initial guess
    solver_guess = []
    solver_guess.extend(initial_guess)
    
    # Create a new temporary working directory for cog
    if twd is None:
        twd = mkdtemp()
        atexit.register(rmtree, twd)
    else:
        assert os.path.exists(twd), twd
    

    ## Careful: caching does NOT work if we use sigma clipping,
    ## have to reset the caches after each clip
    @stellar_optimization
    def minimisation_function(stellar_parameters, *args):
        """The function we want to minimise (e.g., calculates the quadrature
        sum of all minimizable variables.)
        
        stellar_parameters : [teff, logg, vt, feh]
        """
        
        teff, logg, vt, feh = stellar_parameters
        sampled_points, total_tolerance, individual_tolerances, use_nlte_grid = args
        
        if teff < parameter_ranges["teff"][0] or teff > parameter_ranges["teff"][1] or \
           logg < parameter_ranges["logg"][0] or logg > parameter_ranges["logg"][1] or \
           vt < parameter_ranges["vt"][0] or vt > parameter_ranges["vt"][1] or \
           feh < parameter_ranges["[Fe/H]"][0] or feh > parameter_ranges["[Fe/H]"][1] or \
           np.isnan(teff) or np.isnan(logg) or np.isnan(vt) or np.isnan(feh):
            #return np.array([np.nan, np.nan, np.nan, np.nan])
            return np.array([np.inf, np.inf, np.inf, np.inf])
        
        photosphere = photosphere_interpolator(teff, logg, feh, alphafe)
        photosphere.meta["stellar_parameters"]["microturbulence"] = vt
        
        abundances = rt.abundance_cog(photosphere,transitions,twd=twd)
        transitions["abundance"] = abundances
        
        ## Calculate slopes and differences that are being minimized
        x, y = transitions[idx_expot]["expot"], abundances[idx_expot]
        dAdchi, b_expot, medy, stdy, stdm_expot, N = utils.fit_line(x, y, None)
        
        x, y = transitions[idx_rew]["reduced_equivalent_width"], abundances[idx_rew]
        dAdREW, b_rew, medy, stdy, stdm_rew, N = utils.fit_line(x, y, None)
        
        dFe = np.mean(abundances[idx_I]) - np.mean(abundances[idx_II])
        dM  = np.mean(abundances[idx_MH]) - (feh + solar_composition("Fe"))
        
        # Some metric has been imposed here. Could make it depend on the slope standard error, but that seems not stable.
        results = np.array([dAdchi, dAdREW, 0.1 * dFe, 0.1 * dM])
        acquired_total_tolerance = np.sum(results**2)
        
        point = [teff, logg, vt, feh] + list(results)
        sampled_points.append(point)

        acquired_total_tolerance = np.sum(results**2)
        logger.debug("Atmosphere with Teff = {0:.0f} K, logg = {2:.2f}, vt = {1:.2f} km/s, [Fe/H] = {3:.2f}, [alpha/Fe] = {4:.2f}"
                     " yields sum {5:.1e}:\n\t\t\t[{6:.1e}, {7:.1e}, {8:.1e}, {9:.1e}]".format(teff, logg, vt, feh, alphafe,
                                                                                               acquired_total_tolerance, *results))
        return results
    
    all_sampled_points = []

    start = time.time()
    
    for i in range(1, 1 + max_attempts):
        sampled_points = []
        args = (sampled_points, total_tolerance, individual_tolerances, 
                use_nlte_grid)
        
        ## Hack to get the output of the minimization function
        #b_expot = np.nan
        #b_rew = np.nan
        
        ## Minimize function
        try:
            results = fsolve(minimisation_function, solver_guess, args=args, fprime=utils.approximate_stellar_jacobian_2,
                             col_deriv=1, epsfcn=0, xtol=1e-10, full_output=1, maxfev=maxfev)
        except OptimizationSuccess as e:
            logger.info(e)
            t_elapsed = time.time() - start
            point_results = np.sum(np.array(sampled_points)[:, 4:]**2, axis=1)
            min_index = np.nanargmin(point_results)
        else:
            min_index = -1
        
        ## Take the result
        final_parameters = sampled_points[min_index][:4]
        final_parameters[0] = int(np.round(final_parameters[0])) # Effective temperature
        final_parameters_result = sampled_points[min_index][4:]
        num_moog_iterations = len(sampled_points)
        all_sampled_points.extend(sampled_points)
        
        ## Check tolerance
        tolerance_achieved = check_satisfies_tolerance(
            final_parameters_result, total_tolerance, individual_tolerances
        )
        
        ## Sigma Clip Outlier Lines
        abundances = transitions["abundance"]
        to_clip = np.zeros(len(abundances), dtype=bool)
        for species in unique_species:
            ii = (transitions["species"]==species) & valid
            med = np.median(abundances[ii])
            #scale = biweight_scale(abundances[ii])
            raise NotImplementedError("Import errors for others")
            this_clip = ii & (np.abs(abundances-med) > sigma_clip*scale)
            to_clip = to_clip | this_clip
        logger.info("Iteration {}/{}: clipping {} lines".format(i,max_attempts,to_clip.sum()))
        logger.info("Iteration {}/{}: {}".format(i,max_attempts,np.where(to_clip)[0]))
        valid[valid] = np.logical_not(to_clip[valid])
        logger.info("Iteration {}/{}: now {} lines".format(i,max_attempts,valid.sum()))
        # Propagate to idx
        idx_expot = idx_expot & valid
        idx_rew = idx_rew & valid
        idx_I = idx_I & valid
        idx_II = idx_II & valid
        idx_MH = idx_MH & valid
        
        # Success if no more clipping and (OptimizationSuccess or tolerance achieved)
        if to_clip.sum() == 0 and ((min_index != -1) or tolerance_achieved):
            # We did it! This is the main exit point for the optimization loop
            break
        
        if to_clip.sum() == 0 and not tolerance_achieved:
            # This solution will fail. Try a new starting point
            solver_guess = [np.random.uniform(*parameter_ranges[parameter]) for parameter in ["teff", "logg", "vt", "[Fe/H]"]]
            logger.info("No lines to clip and tolerance not achieved: trying random restart at {}".format(solver_guess))
        
        ## If to_clip.sum() > 0, then go through another sigma clipping iteration.
        
        pass
    else:
        warnings.warn("Did not converge in {} attempts".format(max_attempts))
    
    t_elapsed = time.time() - start
    logger.info("Took {:.1f}s, used {}/{} lines; Tol. achieved? {}".format(t_elapsed, valid.sum(), len(valid), tolerance_achieved))
    logger.info("{} -> {}".format(initial_guess, final_parameters))
    return (tolerance_achieved, initial_guess, num_moog_iterations, i, t_elapsed,
            final_parameters, final_parameters_result, np.array(all_sampled_points), valid)

# E. Holmbeck copied the above function and basically changed teff and logg to constant.
def optimize_feh(initial_guess, transitions, params_to_optimize, EWs=None, 
                                jacobian=utils.approximate_feh_jacobian,
                                rt=radiative_transfer.moog,
                                max_attempts=5, total_tolerance=1e-4, 
                                individual_tolerances=None, 
                                maxfev=30, use_nlte_grid=None,
                                ):
    """
    Assumes these are all transitions you want to use for stellar parameters
    Assumes you only want to balance neutral ions against Expot and REW
    
    If specify EWs, all EWs are in same order as transitions

    initial_guess order : [teff, vt, logg, feh]
    """
    
    if EWs is None:
        EWs = transitions["equivalent_width"]
        REWs = np.log10(1e-3 * EWs / transitions['wavelength'])
    else:    
        REWs = np.log10(1e-3 * EWs / transitions['wavelength'])
        transitions["equivalent_width"] = EWs
    transitions["reduced_equivalent_width"] = REWs
    idx_I  = transitions["ion"] == 1
    idx_II = transitions["ion"] == 2

    photosphere_interpolator = photospheres.interpolator()
    
    '''
    parameter_ranges = {
        "teff": (3500, 7000),
        "vt": (0.0, 4.0),
        "logg": (0, 5),
        "[Fe/H]": (-5, 0.5)
        }
	'''

    plist = ['Teff','vt','logg','[Fe/H]']
    nconst = int(4-sum(params_to_optimize))
	
    if nconst == 1:
        logger.info('Holding 1 parameter constant: '+', '.join([plist[p]+'='+str(initial_guess[p]) for p in range(len(plist)) if not params_to_optimize[p]]))
    elif nconst > 1:
        logger.info('Holding %s parameters constant: '%nconst +', '.join([plist[p]+'='+str(initial_guess[p]) for p in range(len(plist)) if not params_to_optimize[p]]))

    parameter_ranges = {}
    if params_to_optimize[0]:
        parameter_ranges["teff"] = (3500, 7000)
    else:
        parameter_ranges["teff"] = (initial_guess[0], initial_guess[0])
	
    if params_to_optimize[1]:
        parameter_ranges["vt"] = (0.0, 4.0)
    else:
        parameter_ranges["vt"] = (initial_guess[1], initial_guess[1])

    if params_to_optimize[2]:
        parameter_ranges["logg"] = (0, 5),
    else:
        parameter_ranges["logg"] = (initial_guess[2], initial_guess[2])

    if params_to_optimize[3]:
        parameter_ranges["[Fe/H]"] = (-5, 0.5)
    else:
        parameter_ranges["[Fe/H]"] = (initial_guess[3], initial_guess[3])
	
    
    # Create a mutable copy for the initial guess
    #solver_guess = []
    #solver_guess.extend(initial_guess)
    solver_guess = initial_guess[params_to_optimize]

    # Create a new temporary working directory for cog
    twd = mkdtemp()
    atexit.register(rmtree, twd)
	
    @feh_optimization
    def minimisation_function(stellar_parameters, *args):
        """The function we want to minimise (e.g., calculates the quadrature
        sum of all minimizable variables.)
        
        stellar_parameters : [teff, vt, logg, feh]
        """
        
        params_to_optimize, all_sampled_points, total_tolerance, individual_tolerances, use_nlte_grid = args

        # Old way:
        #teff, vt, logg, feh = [initial_guess[0], stellar_parameters[0], initial_guess[2], stellar_parameters[1]]

        # First set all to initial_guess, then replace the ones with a "True" value	
        params_list = initial_guess
        next_param = 0
        for pi,p in enumerate(params_to_optimize):
        	if p==True:
        		params_list[pi] = stellar_parameters[next_param]
        		next_param += 1
        
        teff, vt, logg, feh = params_list

        photosphere = photosphere_interpolator(teff, logg, feh)
        photosphere.meta["stellar_parameters"]["microturbulence"] = vt
        
        ## TODO: ADJUST ABUNDANCES TO ASPLUND?
        abundances = rt.abundance_cog(photosphere,transitions,twd=twd)
        transitions["abundance"] = abundances
        
        out = utils.equilibrium_state(transitions[idx_I],
                                      ("expot", "reduced_equivalent_width"))
        dAdchi = out[26.0]['expot'][0]
        dAdREW = out[26.0]['reduced_equivalent_width'][0]
        dFe = np.mean(abundances[idx_I]) - np.mean(abundances[idx_II])
        # E. Holmbeck changed dM to be w.r.t. FeII abundances.
        dM  = np.mean(abundances[idx_II]) - (feh + solar_composition("Fe"))
		
        results = np.array([dAdchi, dAdREW, 0.1 * dFe, 0.1 * dM])
        acquired_total_tolerance = np.sum(results**2)

        point = list(np.array([teff, vt, logg, feh])*params_to_optimize) + list(results*params_to_optimize)
        all_sampled_points.append(point)

        logger.info("Atmosphere with Teff = {0:.0f} K, vt = {1:.2f} km/s, logg = {2:.2f}, [Fe/H] = {3:.2f}, [alpha/Fe] = {4:.2f}"
                    " yields sum {5:.1e}:\n\t\t\t[{6:.1e}, {7:.1e}, {8:.1e}, {9:.1e}]".format(teff, vt, logg, feh, 0.4,
                                                                                              acquired_total_tolerance, *results))
        return results[params_to_optimize]

    all_sampled_points = []
    
    start = time.time()
    for i in range(1, 1 + max_attempts):
        sampled_points = []
        args = (params_to_optimize, sampled_points, total_tolerance, individual_tolerances, 
                use_nlte_grid)
	
        try:
            results = fsolve(minimisation_function, solver_guess, args=args, fprime=utils.approximate_feh_jacobian,
                             col_deriv=1, epsfcn=0, xtol=1e-10, full_output=1, maxfev=maxfev)
        except OptimizationSuccess as e:#: #Exception as e:# OptimizationSuccess as e:
            print(e)
            
            t_elapsed = time.time() - start
            
            sampled_points_tols = np.array(sampled_points)[:,4:]*params_to_optimize

            point_results = np.sum(sampled_points_tols**2, axis=1)
            min_index = np.argmin(point_results)

            final_parameters = (sampled_points[min_index][:4]*params_to_optimize) + (initial_guess*(-1*params_to_optimize + 1))
            final_parameters_result = (sampled_points[min_index][4:]*params_to_optimize) + (initial_guess*(-1*params_to_optimize + 1))

            num_moog_iterations = len(sampled_points)
            all_sampled_points.extend(sampled_points)
            
            acquired_total_tolerance = sum(pow(np.array(final_parameters_result), 2))
            return (True, initial_guess, num_moog_iterations, i, t_elapsed, 
                    final_parameters, final_parameters_result,
                    np.array(all_sampled_points))
        else:
            t_elapsed = time.time() - start
            #ier, mesg = results[-len(sampled_points)/2:]
            ier, mesg = results[-2:]
            
            final_parameters = (results[0]*params_to_optimize) + (initial_guess*(-1*params_to_optimize + 1))
            final_parameters[0] = int(np.round(final_parameters[0])) # Effective temperature
            final_parameters_result = results[1]["fvec"]
            num_moog_iterations = len(sampled_points)
            all_sampled_points.extend(sampled_points)
            acquired_total_tolerance = sum(pow(np.array(final_parameters_result), 2))

            # Have we achieved tolerance?
            if total_tolerance >= acquired_total_tolerance and \
                    (individual_tolerances is None or np.all(np.less_equal(np.abs(final_parameters_result), individual_tolerances))):
                
                return (True, initial_guess, num_moog_iterations, i, t_elapsed, final_parameters, final_parameters_result, np.array(all_sampled_points))
            else:
                # This solution will fail. Try a new starting point
                # It should select a new point where teff and logg are constant though...
                #solver_guess = [teff, np.random.uniform(*parameter_ranges["vt"]), logg, np.random.uniform(*parameter_ranges["[Fe/H]"])]
                solver_guess = [np.random.uniform(*parameter_ranges[parameter]) for parameter in ['teff', 'vt', 'logg', '[Fe/H]']]
                
        t_elapsed = time.time() - start

        if total_tolerance >= acquired_total_tolerance and \
            (individual_tolerances is None or np.all(np.less_equal(np.abs(final_parameters_result), individual_tolerances))):
            tolerance_achieved = True
            
        else:
            tolerance_achieved = False

        return (tolerance_achieved, initial_guess, num_moog_iterations, i, t_elapsed, final_parameters, final_parameters_result, np.array(all_sampled_points))

    pass
