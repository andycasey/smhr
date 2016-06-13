


from base import RadiativeTransfer



class Turbospectrum(RadiativeTransfer):

    name = "Turbospectrum"
    _executables = ("bsyn_lu", )

    def __init__(self, executable_path=None, **kwargs):
        super(Turbospectrum, self).__init__(executable_path, **kwargs)
        self._version = self._get_version()
        return None

    def _get_version(self):
        """ Parse the version from the output of the bsyn_lu executable. """
        code, stdout, stderr = self()
        return stdout.split("\n")[2].strip(" *").split()[-1]


    def synthesize(self, atmosphere, transitions, **kwargs):



        raise a



if __name__ == "__main__":
    rt = Turbospectrum(
        "/Users/arc/codes/turbospectrum/EXPORT-15.1/exec-v15.1/bsyn_lu")
