import numpy as np
from scipy.optimize import fsolve
from . import (photospheres, radiative_transfer, utils)
from functools import wraps
import sys, os, time

from smh.photospheres.abundances import asplund_2009 as solar_composition

import logging
logger = logging.getLogger(__name__)


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
            sampled_before = (previously_sampled_points[:, :4] == this_point).all(
                np.arange(previously_sampled_points[:, :4].ndim - this_point.ndim, previously_sampled_points[:, :4].ndim))

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
                                maxfev=30, use_nlte_grid=None):
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
        
        photosphere = photosphere_interpolator(teff, logg, feh)
        photosphere.meta["stellar_parameters"]["microturbulence"] = vt
        
        ## TODO: ADJUST ABUNDANCES TO ASPLUND?
        abundances = rt.abundance_cog(photosphere,transitions)
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
                    " yields sum {5:.1e}:\n\t\t\t[{6:.1e}, {7:.1e}, {8:.1e}, {9:.1e}]".format(teff, vt, logg, feh, 0.4,
                                                                                              acquired_total_tolerance, *results))
        return results

    all_sampled_points = []
        
    start = time.time()
    for i in xrange(1, 1 + max_attempts):
        sampled_points = []
        args = (sampled_points, total_tolerance, individual_tolerances, 
                use_nlte_grid)

        #results = fsolve(minimisation_function, solver_guess, args=args, fprime=utils.approximate_stellar_jacobian,
        #                 col_deriv=1, epsfcn=0, xtol=1e-10, full_output=1, maxfev=maxfev)
        try:
            results = fsolve(minimisation_function, solver_guess, args=args, fprime=utils.approximate_stellar_jacobian,
                             col_deriv=1, epsfcn=0, xtol=1e-10, full_output=1, maxfev=maxfev)
        except: #Exception as e:# OptimizationSuccess as e:
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

