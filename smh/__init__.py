
import logging
import os
from shutil import copyfile

__version__ = "0.1"

# Set up logging.
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter(
    "%(asctime)s [%(levelname)-8s] %(message)s"))

#handler.setFormatter(logging.Formatter(
#    "%(asctime)s [%(levelname)-8s] (%(name)s/%(lineno)d): %(message)s"))

logger.addHandler(handler)


# Get the location for the Session defaults file.
from .session import Session
from . import (photospheres, radiative_transfer, spectral_models)

# If there isn't a local copy of the default Session settings file, create one.
if not os.path.exists(Session._default_settings_path):
    copyfile(
        os.path.join(
            os.path.dirname(os.path.join(__file__)),
            "../default_session.yaml"
        ),
        Session._default_settings_path)
