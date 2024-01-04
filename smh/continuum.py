
import numpy as np
import gzip
import warnings
import pickle
import os
from functools import cache
from typing import Optional, Tuple
from scipy import optimize as op
from sklearn.decomposition._nmf import _fit_coordinate_descent
from sklearn.exceptions import ConvergenceWarning
from scipy.signal.windows import tukey
from scipy.optimize import leastsq
from scipy import interpolate


SMALL = 1e-12

def expand_path(path):
    return os.path.abspath(os.path.expanduser(path))

@cache
def load_basis_vectors(path, P, pad=0):    
    full_path = expand_path(path)
    if full_path.lower().endswith(".gz"):
        with gzip.open(full_path, "rb") as fp:
            masked_basis_vectors = pickle.load(fp)
    else:        
        with open(full_path, "rb") as fp:
            masked_basis_vectors = pickle.load(fp)
    if pad > 0:
        basis_vectors = np.zeros((masked_basis_vectors.shape[0], P + 2 * pad))
        basis_vectors[:, pad:-pad] = masked_basis_vectors
    else:
        basis_vectors = masked_basis_vectors
    return basis_vectors


def fit_sines_and_cosines(wavelength, flux, ivar, deg: Optional[int] = 3, L: Optional[float] = None, A = None):
    L = L or np.ptp(wavelength)
    if A is None:
        A = design_matrix(wavelength, deg)
    MTM = A @ (ivar[:, None] * A.T)
    MTy = A @ (ivar * flux)
    try:
        theta = np.linalg.solve(MTM, MTy)
    except np.linalg.LinAlgError:
        if np.any(ivar > 0):
            raise
    continuum = A @ theta
    return (continuum, theta, A)


def template_logwl_resample(
    wavelength, 
    flux,
    template_wavelength,
    template_flux,
    wblue=None, 
    wred=None,
    delta_log_wavelength=None,
):
    """
    Resample a spectrum and template onto a common log-spaced spectral grid.
    """

    if wblue:
        w0 = np.log10(wblue)
    else:
        ws0 = np.log10(wavelength[0])
        wt0 = np.log10(template_wavelength[0])
        w0 = min(ws0, wt0)

    if wred:
        w1 = np.log10(wred)
    else:
        ws1 = np.log10(wavelength[-1])
        wt1 = np.log10(template_wavelength[-1])
        w1 = max(ws1, wt1)

    if delta_log_wavelength is None:
        ds = np.log10(wavelength[1:]) - np.log10(wavelength[:-1])
        dw = ds[np.argmin(ds)]
    else:
        dw = delta_log_wavelength

    nsamples = int((w1 - w0) / dw)

    log_wave_array = np.ones(nsamples) * w0
    for i in range(nsamples):
        log_wave_array[i] += dw * i

    # Build the corresponding wavelength array
    wave_array = np.power(10., log_wave_array)

    # Resample spectrum and template into wavelength array so built
    resampled_flux = np.interp(wave_array, wavelength, flux, left=0, right=0)
    resampled_template_flux = np.interp(wave_array, template_wavelength, template_flux, left=0, right=0)
    
    return (wave_array, resampled_flux, resampled_template_flux)



def cross_correlate(wavelength, flux, template_wavelength, template_flux, apodization_window=0.2, resample=True):
    if resample:
        if resample is True:
            resample_kwargs = dict()  # use defaults
        else:
            resample_kwargs = resample
        resampled_wavelength, resampled_flux, resampled_template_flux = template_logwl_resample(
            wavelength, flux, template_wavelength, template_flux, **resample_kwargs)

    else:
        resampled_wavelength = wavelength
        resampled_flux = flux
        resampled_template_flux = template_flux

    resampled_flux, resampled_template_flux = _apodize(resampled_flux, resampled_template_flux, apodization_window)

    resampled_flux -= np.nanmean(resampled_flux)
    resampled_template_flux -= np.nanmean(resampled_template_flux)

    corr = np.correlate(resampled_flux, resampled_template_flux, mode='full')

    delta_log_wave = np.log10(resampled_wavelength[1]) - np.log10(resampled_wavelength[0])
    deltas = (np.array(range(len(corr))) - len(corr)/2 + 0.5) * delta_log_wave
    lags = np.power(10., deltas) - 1.
    return (corr, lags * 299792.458)



def measure_relative_velocity(wavelength, flux, template_wavelength, template_flux, v_lim=(-250, 250), deg=1, N=10, **kwargs):
    
    corr, lags = cross_correlate(wavelength, flux, template_wavelength, template_flux, **kwargs)
    
    si, ei = lags.searchsorted(v_lim)
    pi = si + np.argmax(corr[si:ei])
    v_rel = lags[pi]
    
    # Get initial guess of peak.
    p0 = np.array([v_rel, np.max(corr[si:ei]), 10] + [0] * deg)

    gaussian = lambda p, x: p[1] * np.exp(-(x - p[0])**2 / (2.0 * p[2]**2)) + np.polyval(p[3:], x)
    errfunc = lambda p, x, y: y - gaussian(p, x)

    # only consider +/- N pixels
    ps, pe = (pi - N, pi + N + 1)
    try:
        p1, ier = leastsq(errfunc, p0.copy(), args=(lags[ps:pe], corr[ps:pe]))
    except:
        raise 
    
    v_rel, e_v_rel = (p1[0], p1[2]) 
    chi2 = -2 * np.max(corr[si:ei])
    return (v_rel, e_v_rel, chi2, corr, lags)
    


def _apodize(spectrum, template, apodization_window):
    # Apodization. Must be performed after resampling.
    if apodization_window is None:
        clean_spectrum = spectrum
        clean_template = template
    else:
        if callable(apodization_window):
            window = apodization_window
        else:
            def window(wlen):
                return tukey(wlen, alpha=apodization_window)
        clean_spectrum = spectrum * window(len(spectrum))
        clean_template = template * window(len(template))

    return clean_spectrum, clean_template



class BaseContinuumModel(object):

    def __init__(
        self,
        wavelength: np.array,
        basis_vectors: np.ndarray,
        deg: Optional[int] = 3,
    ):       
        self.wavelength = wavelength
        _check_wavelength_basis_vectors_shape(wavelength, basis_vectors)
        self.basis_vectors = basis_vectors
        self.deg = deg
        return None

    def fit_sines_and_cosines(self, wavelength, flux, ivar, A=None):
        if A is None:
            A = design_matrix(wavelength, self.deg)
        MTM = A @ (ivar[:, None] * A.T)
        MTy = A @ (ivar * flux)
        try:
            theta = np.linalg.solve(MTM, MTy)
        except np.linalg.LinAlgError:
            if np.any(ivar > 0):
                raise
        continuum = theta @ A
        return (continuum, theta, A)

    def _W_step(self, mean_rectified_flux, W, **kwargs):
        absorption = 1 - mean_rectified_flux
        use = np.ones(mean_rectified_flux.size, dtype=bool)
        use *= (
            np.isfinite(absorption) 
        &   (absorption >= 0) 
        &   (mean_rectified_flux > 0)
        )
        W_next, _, n_iter = _fit_coordinate_descent(
            absorption[use].reshape((1, -1)),
            W,
            self.basis_vectors[:, use],
            update_H=False,
            verbose=False,
            shuffle=True
        )        
        rectified_model_flux = 1 - (W_next @ self.basis_vectors)[0]
        return (W_next, rectified_model_flux, np.sum(use), n_iter)


    def _predict(self, theta, A_slice, C, P):
        return (1 - theta[:C] @ self.basis_vectors) * (A_slice @ theta[C:]).reshape((-1, P))

    
    def full_design_matrix(self, N, R=1):        
        C, P = self.basis_vectors.shape
        
        K = R * (2 * self.deg + 1)
        A = np.zeros((N * P, C + N * K), dtype=float)
        for i in range(N):
            A[i*P:(i+1)*P, :C] = self.basis_vectors.T
            A[i*P:(i+1)*P, C + i*K:C + (i+1)*K] = design_matrix(self.wavelength, self.deg).T
        return A

    def get_mask(self, ivar):
        N, P = np.atleast_2d(ivar).shape        
        use = np.zeros((N, P), dtype=bool)
        #use[:, np.hstack(self.region_masks)] = True
        use *= (ivar > 0)
        return ~use

    

    def _get_initial_guess_with_small_W(self, si, ei, flux, ivar, initial_theta=None, small=1e-12):
        with warnings.catch_warnings():        
            warnings.filterwarnings("ignore", category=RuntimeWarning)
            N = 1
            C, P = self.basis_vectors.shape
            R = 1
            P = ei - si
            
            K = R * (2 * self.deg + 1)
            A = np.zeros((N * P, C + N * K), dtype=float)
            for i in range(N):
                A[i*P:(i+1)*P, :C] = self.basis_vectors[:, si:ei].T
                A[i*P:(i+1)*P, C + i*K:C + (i+1)*K] = design_matrix(self.wavelength[si:ei], self.deg).T


            Y = flux.flatten()
            use = ~self.get_mask(ivar).flatten()
            result = op.lsq_linear(
                A[use],
                Y[use],
                bounds=self.get_bounds(A.shape[1], [-np.inf, 0]),
            )            
            C, P = self.basis_vectors.shape       
            theta = result.x[C:]     
            if initial_theta is not None and np.allclose(theta, np.zeros_like(theta)):
                theta = initial_theta
            return np.hstack([small * np.ones(C), theta])
                
    
    def _get_initial_guess(self, si, ei, flux, continuum, initial_theta):
                
        absorption = np.clip(1 - flux / continuum, 0, 1)
        use = np.ones(absorption.size, dtype=bool)
        use *= (
            np.isfinite(absorption) 
        &   (absorption >= 0) 
        )
        W = SMALL * np.ones((1, self.basis_vectors.shape[0]))
        
        W_next, _, n_iter = _fit_coordinate_descent(
            absorption[use].reshape((1, -1)),
            W,
            self.basis_vectors[:, si:ei][:, use],
            update_H=False,
            verbose=False,
            shuffle=True
        )               
        return np.hstack([W_next.flatten(), initial_theta])




    def fit(self, wavelength, flux, ivar, v_rel=0.0):
        """
        Fit the continuum in every order simultaneously.
        """
        
        # interpolate the basis vectors onto the observed wavelength arrays
        
        # create design matrix for each order
        raise NotImplementedError
        
        


    def fit_relative_velocity(self, wavelength, flux, ivar, x0=None, **kwargs):
        
        # TODO: check that wavelength, flux, etc are 1D arrays only
    
        A_obs = design_matrix(wavelength, self.deg)
        (initial_continuum, initial_theta, _) = self.fit_sines_and_cosines(wavelength, flux, ivar, A=A_obs)

        si, ei = np.searchsorted(self.wavelength, wavelength[[0, -1]])
        
        C, P = self.basis_vectors.shape
        v_rels, chi2s = (np.nan * np.ones(C), np.nan * np.ones(C))
        for i, basis_vector in enumerate(self.basis_vectors[:, si:ei]):
            v_rel, e_v_rel, chi2, corr, lags = measure_relative_velocity(
                wavelength,
                flux / initial_continuum,
                self.wavelength[si:ei],
                1 - basis_vector/np.max(basis_vector)
            )
            chi2s[i] = chi2
            v_rels[i] = v_rel
        
        initial_v_rel = v_rels[np.argmin(chi2s)]
        
        # interpolate data to model (this is bad practice, but all of this is approximate)
        si, ei = self.wavelength.searchsorted(overlap(self.wavelength, wavelength))
        interp_wavelength = self.wavelength[si:ei]
        interp_flux = np.interp(interp_wavelength, wavelength * (1 - initial_v_rel/3e5), flux, left=0, right=0)
        interp_ivar = np.interp(interp_wavelength, wavelength * (1 - initial_v_rel/3e5), ivar, left=0, right=0)
        interp_continuum = np.interp(interp_wavelength, wavelength * (1 - initial_v_rel/3e5), initial_continuum, left=0, right=0)
        
        sigma = interp_ivar**-0.5
        
        x0 = self._get_initial_guess(si, ei, interp_flux, interp_continuum, initial_theta)
        
        A = design_matrix(interp_wavelength, self.deg).T

        basis_vectors_slice = self.basis_vectors[:, si:ei]
        
        def f(_, *theta):            
            return (1 - theta[:C] @ basis_vectors_slice) * (A @ theta[C:]).flatten()
        
        p_opt, cov = op.curve_fit(
            f,
            None,
            interp_flux,
            p0=x0,
            sigma=sigma,
            bounds=self.get_bounds(C + 2 * self.deg + 1)
        )

        '''
        fig, ax = plt.subplots()
        ax.plot(interp_wavelength, interp_flux, c='k')
        ax.plot(interp_wavelength, f(None, *x0), c="tab:red")
        ax.plot(interp_wavelength, f(None, *x1), c="tab:blue")
        ax.plot(interp_wavelength, f(None, *p_opt), c="tab:green")
        '''

        rectified_model_flux = 1 - p_opt[:C] @ basis_vectors_slice
        # compute the continuum using the original design matrix
        # TODO: this is not quite right... change it all to use design matrices in one reference frame!!
        continuum = np.interp(wavelength * (1 - initial_v_rel/3e5), interp_wavelength, p_opt[C:] @ A.T, left=0, right=0)
        
        
        # measure the radial velocity again
        v_rel, e_v_rel, chi2, corr, lags = measure_relative_velocity(
            wavelength,
            flux / continuum,
            interp_wavelength,
            rectified_model_flux
        )
        
        return (v_rel, continuum, interp_wavelength, rectified_model_flux, p_opt)


    def get_bounds(self, N, component_bounds=(0, +np.inf)):
        C, P = self.basis_vectors.shape          

        return np.vstack([
            np.tile(component_bounds, C).reshape((C, 2)),
            np.tile([-np.inf, +np.inf], N - C).reshape((-1, 2))
        ]).T            


def overlap(A, B):
    # get range of values in A and B that overlap
    return (np.max([A[0], B[0]]), np.min([A[-1], B[-1]]))






def _check_wavelength_basis_vectors_shape(wavelength, basis_vectors):
    P = wavelength.size
    assert wavelength.ndim == 1, "wavelength must be a one-dimensional array." 
    C, P2 = basis_vectors.shape
    assert P == P2, "`basis_vectors` should have shape (C, P) where P is the size of `wavelength`"


def design_matrix(wavelength: np.array, deg: int) -> np.array:
    L = np.ptp(wavelength)
    scale = 2 * (np.pi / L)
    return np.vstack(
        [
            np.ones_like(wavelength).reshape((1, -1)),
            np.array(
                [
                    [np.cos(o * scale * wavelength), np.sin(o * scale * wavelength)]
                    for o in range(1, deg + 1)
                ]
            ).reshape((2 * deg, wavelength.size)),
        ]
    )




'''
class ContinuumModel(BaseContinuumModel):
    
    def __init__(
        self,
    ):
        super(BaseContinuumModel, self).__init__(
            wavelength=wavelength,
            basis_vectors=basis_vectors,
        )
        return None    
'''

if __name__ == "__main__":
    
    
    from specutils import Spectrum1D
    
    #wavelength, flux, ivar, meta = Spectrum1D.read_fits_multispec("j174239-133332blue_multi.fits")
    wavelength, flux, ivar, meta = Spectrum1D.read_fits_multispec(expand_path("~/Downloads/hd122563red_multi.fits"), flux_ext=6, ivar_ext=3)
    
    index = 14
    model_wavelength = 10 * (10**(2.57671464 + np.arange(167283) * 2.25855074e-06))
    basis_vectors = load_basis_vectors("H.pkl.gz", wavelength.size)
    
    from time import time
    model = BaseContinuumModel(model_wavelength, basis_vectors, deg=20)

    t_init = time()
    (v_rel, continuum, template_wavelength, template_flux, p_opt) = model.fit_relative_velocity(wavelength[index], flux[index], ivar[index])
    t_fit = time() - t_init
    print(v_rel, t_fit)
    
    fig, ax = plt.subplots()
    ax.plot(wavelength[index], flux[index] / continuum, c='#cccccc')
    ax.plot(wavelength[index] * (1 - v_rel/3e5), flux[index] / continuum, c='k')
    ax.plot(template_wavelength, template_flux, c="tab:red")
    ax.set_ylim(0, 1.2)
    
    raise a

    continuum, theta, A = model.fit_sines_and_cosines(wavelength[index], flux[index], ivar[index])
    
    si, ei = np.searchsorted(model.wavelength, wavelength[index][[0, -1]])
    
    wl = model.wavelength[si:ei]
    
    vals = []
    for ti in range(32):
        template_flux = model.basis_vectors[ti, si:ei]
        template_flux = 1 - template_flux/np.max(template_flux)
        
        v_rel, e_v_rel, chi2, corr, lags = measure_relative_velocity(wavelength[index], flux[index] / continuum, model.wavelength[si:ei], template_flux)
        print(ti, v_rel, chi2)
        vals.append([v_rel, chi2])
    
    vals = np.array(vals)
    ti = np.argmin(vals[:, 1])
    v_rel, chi2 = vals[ti]

    # 


    fig, ax = plt.subplots()
    ax.plot(lags, corr)
    
    fig, ax = plt.subplots()
    ax.plot(wavelength[index], flux[index] / continuum, c='#666666')
    ax.plot(wavelength[index] * (1 - v_rel/3e5), flux[index] / continuum, c="k")
    template_flux = model.basis_vectors[ti, si:ei]
    template_flux = 1 - template_flux/np.max(template_flux)
    
    ax.plot(model.wavelength[si:ei], template_flux, c="tab:red")
    ax.set_ylim(0, 5)
    
    # Now fit with the continuum model.
    
