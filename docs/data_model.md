# Data model description for a session

Overview
--------
* New SMH files will be gzipped tar balls containing many files.
* The primary file will be `session.pkl` (Pickle protocol 2 for Legacy Python) which contains all of the descriptive information about the session, and associated measurements.
* The tarball will also contain the input spectra used to initiate the session: these should never be directly edited.
* There may be ancillary files that will also be included in the tarball (e.g., the spectrum used for cross-correlation)


`session.pkl` description
-------------------------

Top-level dictionary containing the following keys:

- `input_spectra_paths`: A list of basenames which should be in the session tar ball.
- `discarded_orders`: indices of any orders (read from `input_spectrum_filenames`) that were ignored/deleted.



Other stuff that will need to be in there, but for which nomenclature and data structure needs to be determined:

- `rv`: Radial velocity
  - **Measurements and related information that was not explicitly inputted**
  - `rv_applied`: The radial velocity that was actually applied
  - `rv_measured`: The last radial velocity that was measured (all other items like the `CCF`, etc will refer to this measurement)
  - `rv_uncertainty`: The uncertainty in the measured radial velocity.
  - `order_index`: The zero-th indexed order used for RV determination by CCF
  - `normalized_order`: Normalized echelle order used to measure radial velocity.
  - `ccf`: The CCF.
  - `heliocentric_correction`: The heliocentric correction calculated based on header information.
  - `barycentric_correction`: The barycentric correction calculated based on header information.
  - **Input settings**
  - `template_spectrum`: The `Spectrum1D` object used to perform the cross-correlation.
  - `wavelength_region`: The (start, end) region used for CCF.
  - `resample`: Whether the 'template' or 'observed' spectrum was resampled.
  - `apodize`: The apodization fraction employed.
  - `normalization`:
    - `knot_spacing`: The knot spacing
    - `sigma_clip`: The low- and high-sigma clipping limits used.
    - `max_iterations`: The maximum number of iterations employed for normalization.
    - `order`: The order of the function used for normalization.
    - `exclude`: Regions that were excluded from the normalization.
    - `include`: Regions that should *always* be included (e.g., additional points), even if they exist within a larger excluded region.
    - `function`: The function used for continuum normalization.
    - `scale`: The scaling factor applied to additional points.


- Continuum (for each order in the input spectra -- do not overlook discarded orders! they are just ignored, not deleted):
  - pixel-by-pixel continuum calculated for each rest-frame (wavelength-shifted, no rebinning yet!) order, including the settings used to determine the continuum (e.g., masks, method, additional points + weights etc)
  - NB: Because the RV correction happens first, there may be default continuum masks in SMH. In this case we store only the mask relevant to the wavelength extents of each order. 

- Stacked spectrum:
  - The rebinned, stacked normalised rest-frame spectrum. Include any rebinning options (e.g., bin size, linear/log binning)
 
- Input physics:
  - model atmospheres used, which version, the ADS reference, and a hash of the entire model grid in case there are unexpected updates to check against.
  - radiative transfer code used, version, and default input settings for each driver
  - spectral_lines:
    - a full list of **unique** spectral lines. One line may actually link to a line list where there is a single (or multiple!) element to be measured
    - the line list should include everything (damping coefficients, lower and upper levels, multiplet description) because different RT codes use different things
    - When each line is added, a md5 hash is made of its atomic data. That will make the line unique, and any measurements (in another section) will refer to the hash.
  

- spectral line measurements:
  - These can be measured by EW, Synth, or in principle some other way. Here are the things to store:
  - the hash (see above) of the spectral line that this measurement was made for
  - the best-fitting model spectrum (be it synth or a Gaussian profile)
  - the abundances of all (relevant) elements employed for the fit. If this is an EW measurement then these will just be NaN.
  - Any nuisance parameters solved at the time of fitting: residual RV (synth or EW), smoothing kernel (synth), continuum params, etc
  - Any mask used for the input 
  - method (Ew/Synth/other) and relevant input parameters (including C12/C13 isotope ratios for synth, etc)
  - whether this line is 'acceptable' or not
  - comments on each measurement (these can be exported as footnotes later in LaTeX tables)
  - whether the measurement should be considered as a limit or not
  - associated uncertainties (formal ones from the fitting procedure + ones calculated upstream at stellar parameter stage) 

NB: This allows for multiple line fits for a single spectral line (e.g., one EW and one Synth). We will just specify default behaviour under-the-hood in SMH so that any new EW fit for a spectral line will just overwrite a previous EW fit for that same spectral line.

- stellar parameters:
  - we should allow for these to be determined in different ways (e.g., excitation/ionisation, Teff from H-lines, whatever)
  - store the method (e.g., excitation/ionisation)
  - inputs (e.g., spectral line hashes)
  - TODO: Need to consider this more in terms of how multiple stellar parameter methods will interplay with each other. Eg if we have teff from excitation/ionisation and H-line fitting, which should be used? Which should be shown as default?


