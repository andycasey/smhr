

import smh
from smh import linelists

session = smh.Session([
    "/Users/arc/codes/smh/hd44007red_multi.fits",
    "/Users/arc/codes/smh/hd44007blue_multi.fits",
])

spectrum = smh.specutils.Spectrum1D.read("../smh/hd44007-rest.fits")


# Check against best GES node using MOOG.
transitions = smh.linelists.LineList.read("/Users/arc/research/ges/linelist/vilnius.ew")

# Generate photosphere.
import smh.photospheres
interpolator = smh.photospheres.interpolator(kind="marcs")
photosphere = interpolator(5737, 4.395, 0.03)
photosphere.meta["stellar_parameters"]["microturbulence"] = 0.87

import smh.radiative_transfer as rt
#foo = rt.moog.abundance_cog(photosphere, transitions)


y2 = smh.linelists.LineList.read("../smh/smh/data/linelists/yII_5320.txt")
li = smh.linelists.LineList.read("../smh/smh/data/linelists/lin_li6707")
#rt.moog.synthesize(photosphere, y2)

single_trans = transitions[[30]]

#rt.moog.synthesize(photosphere, transitions)
import smh.spectral_models as sm

single_line = sm.SpectralSynthesisModel(transitions[[172]], session, "Fe")
li_region = sm.SpectralSynthesisModel(li, session, "Li")

y_region = sm.SpectralSynthesisModel(y2, session, ("Y", "Ca"))


raise a

session.metadata["stellar_parameters"] = {
    "effective_temperature": 5522,
    "metallicity": -2.36,
    "surface_gravity": 3.58,
    "microturbulence": 1.07
}


