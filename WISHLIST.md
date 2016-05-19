

**WISHLIST**




*Long-term issues which are high priority (e.g., improve gradually over time):*
- [X] Travis to be faster (e.g., temporarily remove Py3 dependency)
- [X] Coveralls running
- [ ] Unit test improvements
- [ ] Document string improvements
- [ ] Some documentation on the data model.


*Hard*:
- [X] Widget for when an exception occurs, asking to submit GH issue. Catch all GUI-related exceptions with this. #gui
- [ ] Fix matplotlib visualization bug for figures on the RHS. #gui
- [ ] Edit masks GUI for normalization tab. #gui
- [ ] I/O: Save/load new session format.
- [ ] I/O: Create session .from_filename.
- [ ] Line list manager widget --> transitions/spectral_models/assigned to #category.
- [ ] Can the spectral_models __init__ arguments be abstracted more from the session?
- [ ] Include interactive fitting options in the spectral_models management GUI. #gui
- [ ] Functionality to export all intermediate-step data to maximize reproducibility.
- [ ] Script to format old SMH session files into the new format.
- [ ] Sketch up flexible GUI to be able to use stellar parameter information from other sources (e.g., teff from H-alpha, logg from isochrones). Similarly, ensure this is in context of *where* the H-alpha, isochrone stuff can be inputted.
- [ ] Do spectral_model fitting ordered by setup complexity.
- [ ] Summarise chemical abundances from all spectral_models, such that the abundances used for RT are updated *and* this can be reduced again and exported to FITS tables.
- [ ] Parse radiative transfer options from the session file and update that directly in the session.
- [ ] Optimization of stellar parameters given spectral models, including covariance matrix + jacobian.
- [ ] Uncertainties in stellar parameters given uncertainties in linear fits of excitation balance.
- [ ] Propagate uncertainties from stellar parameters to all "acceptable" spectral models.

*Medium*:
- [ ] An edit list widget for wavelength regions in the RV tab. #gui
- [ ] Allow click-to-mask + additional point functionality in RV tab. #gui
- [ ] Implement literature comparison stuff based on data description in the defaults.yaml file.
- [ ] Show 'heliocentric correction for UTDATE' info in RV tab if it could be calculated.
- [ ] Ensure inverse variance arrays are being properly used in continuum fitting, and that the fits are correct at the edges.
- [ ] 'Open recent' list in the File menu.
- [ ] Fix white-text/blue-button GUI bug in OSX. #gui
- [ ] Quality control metrics filter (and GUI) for spectral models. #gui
- [ ] Refactor method names in tabs to make them more consistent (e.g., _populate_widgets, redraw_*, __init_ui__)
- [ ] Least-recently-used cacher for the stellar photosphere.
- [ ] A widget to select element(s) from a list that will be part of a synthesis.


*Easy*:
- [ ] Configurable summary plot.
- [ ] Fix spacing between top and middle figure in RV tab. #gui
- [ ] Consistent color scheme between all tabs, which is read from the session defaults. #gui
- [ ] Go through all 'TODO', 'HACK', or 'MAGIC' entries in the code and create
      GitHub issues. Assign labels and people.
- [X] Implement 'Query Simbad..' button.
- [ ] Change font style for name and position data in Summary tab. #gui
- [ ] Compile line lists together based on current "best" knowledge.
- [ ] Integrate spectral_models widget (sans any models) into the chemical abundances tab. #gui
