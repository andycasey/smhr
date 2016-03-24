

import os
from shutil import copyfile

# Get the location for the Session defaults file.
from session import Session

# If there isn't a local copy of the default Session settings file, create one.
if not os.path.exists(Session._default_settings_path):
    copyfile(
        os.path.join(
            os.path.dirname(os.path.join(__file__)),
            "../default_session.yaml"
        ),
        Session._default_settings_path)
