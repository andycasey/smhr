from smh.session import Session
from smh.linelists import LineList
from smh import spectral_models
from smh import specutils

star = Session([
    "j1808-5104blue_multi.fits",
    "j1808-5104red_multi.fits",
])

# Measure and correct the radial velocity.
rv, rv_uncertainty = star.rv_measure("hd140283.fits")
#star.rv_correct(rv)

# Fake a normalized spectrum.
star.normalized_spectrum = specutils.Spectrum1D.read("hd140283.fits")

# Create a few spectral models.
transitions = LineList.read("/Users/arc/research/ges/linelist/vilnius.ew")
star.metadata["line_list"] = transitions

star.metadata["spectral_models"] = [
    spectral_models.ProfileFittingModel(star, transitions["hash"][[5]]),
    spectral_models.ProfileFittingModel(star, transitions["hash"][[6]]),
    spectral_models.ProfileFittingModel(star, transitions["hash"][[7]]),
    spectral_models.SpectralSynthesisModel(star, transitions["hash"][[15]], ("Si", ))
]


