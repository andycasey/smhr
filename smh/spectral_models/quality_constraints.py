#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Quality constraints for spectral models. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["constraint", "constraints"]

import logging
import numpy as np

logger = logging.getLogger(__name__)

# Convenience function for accessing spectral model metadata.
_meta = lambda model: model.metadata.get("fitted_result", [{}])[-1]

_common_constraints = {
    "wavelength": lambda model: model.wavelength,   

    # Note: Abundance and abundance uncertainty constraints are applied
    #       to *all* abundances in that spectral model (e.g., spectral
    #       synthesis models with multiple elements to fit)

    "abundance": lambda model: _meta(model)["abundances"],
    "abundance_uncertainty": lambda model: \
        np.abs(_meta(model)["abundance_uncertainties"]),

}

_profile_constraints = {
    "equivalent_width": lambda model: 1e3 * _meta(model)["equivalent_width"][0],
    "equivalent_width_uncertainty": lambda model: \
        1e3 * np.nanmax(np.abs(_meta(model)["equivalent_width"][1:])),
    "equivalent_width_percentage_uncertainty": lambda model: \
        100 * (np.nanmax(np.abs(_meta(model)["equivalent_width"][1:])) \
            / _meta(model)["equivalent_width"][0]),
    "reduced_equivalent_width": lambda model: \
        _meta(model)["reduced_equivalent_width"][0],

    # Constraints related to the transition properties.
    "loggf": lambda model: model.transitions["loggf"][0],
    "excitation_potential": lambda model: model.transitions["expot"][0],
}

_synthesis_constraints = {}

# Join them all together.
_get_data_from_model = {}
[_get_data_from_model.update(_constraint) for _constraint in \
    (_common_constraints, _profile_constraints, _synthesis_constraints)]


def constraint(spectral_model, constraint_name, lower, upper, full_output=False):
    """
    Returns a boolean value indicating whether the supplied spectral model meets
    a constraint.

    :param spectral_model:
        A spectral model.

    :param constraint_name:
        A valid constraint name specified in the 
        `_get_data_from_model` function above.

    :param lower:
        A (inclusive) lower bound value.

    :param upper:
        An (inclusive) upper bound value.

    :param full_output: [optional]
        Optionally also return the model value for this constraint.

    :returns:
        `True` if the model meets the constraint, `False` if it doesn't. If
        `full_output` is used, a two-length tuple containing 
        `(is_ok, model_value)` will be returned.
    """

    lower, upper = (lower or -np.inf, upper or +np.inf)

    if not np.any(np.isfinite([lower, upper])):
        # Nothing to check.
        return True if not full_output else (True, None)

    try:
        func = _get_data_from_model[constraint_name]

    except:
        raise KeyError(
            "unrecognized constraint name '{}' (available: {})".format(
                constraint_name, ", ".join(get_model_data.keys())))

    try:
        value = func(spectral_model)

    except KeyError:
        # That entry of data does not exist in this model.
        return True if not full_output else (True, None)

    # If we get back an array of values from get_model_data, we will
    # apply the constraint to *every* value. All values must fit within
    # the constraints.
    qualifier = np.all if isinstance(value, np.ndarray) else np.any    
    is_ok = qualifier((upper >= value) * (value >= lower))
    return is_ok if not full_output else (is_ok, value)


def constraints(spectral_model, quality_constraints):
    """
    Returns a boolean flag as to whether the supplied spectral model meets a
    set of quality constraints.

    :param spectral_model:
        A spectral model.

    :param quality_constraints:
        A dictionary containing constraint names as keys and a 2-length tuple
        with the (lower, upper) bounds specified as values.

    :returns:
        `True` if the model meets the constraints, `False` if it doesn't.
    """

    for constraint_name, (lower, upper) in quality_constraints.items():
        is_ok, value = constraint(
            spectral_model, constraint_name, lower, upper, full_output=True)
        if not is_ok:
            logger.debug(
                "Spectral model {} does not meet quality constraint {}: "\
                "({} >= {} >= {}) not met".format(
                    spectral_model, constraint_name, lower, value, upper))
            return False

    return True