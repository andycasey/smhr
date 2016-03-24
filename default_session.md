Default SMH settings file
-------------------------

This file will be shipped with SMH, but can be edited by the user.

# The `rv` settings that end up getting used will go into the data model for that `Session` object.
- `rv`
  - `wavelength_regions`: A list giving (start, end) regions preferred for the cross-correlation
    - [8450, 8700] # Ca II triplet
    - [5100, 5200] # Mg I triplet
    - [6510, 6610] # H-alpha
    - [4810, 4910] # H-beta
    - [4290, 4390] # H-gamma
  - `template_path`: A SMH-relative path containing a normalized spectrum to use for comparison
  - `resample`: `template` # or `observed`
  - `apodize`: 0.10 # TODO: Currently not implemented in SMH but maybe not needed?
  - `normalization`:
    - `knot_spacing`: 200
    - `sigma_clip`: (1.0, 0.2)
    - `max_iterations`: 3
    - `order`: 3
    - `exclude`: None
    - `include`: None
    - `function`: 'spline'
    - `scale`: 1.0