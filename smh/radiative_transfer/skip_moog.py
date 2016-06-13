


from base import RadiativeTransfer



class MOOG(RadiativeTransfer):

    # Executables to search for if no executable_path is given when initiated.
    name = "MOOGSILENT"
    _executables = ("moogsilent", "MOOGSILENT")

    def __init__(self, executable_path=None, **kwargs):
        super(MOOG, self).__init__(executable_path, **kwargs)
        self._version = self._get_version()
        return None


    def _get_version(self):
        """ Parse the version from the output of the MOOG executable. """
        code, stdout, stderr = self()
        index = stdout.find("VERSION")
        return (" ".join(stdout[index:].split(" ")[1:3])).strip("()")


    def __call__(self, input_filename=None, **kwargs):
        """
        Execute MOOG.

        :param input_filename: [optional]
            An input file to execute.
        """

        kwds = kwargs.copy()
        if input_filename is not None:
            kwds.update({
                "cwd": os.path.dirname(input_filename),
                "stdin": "{}\n".join(os.path.basename(input_filename))
            })
        return super(MOOG, self).__call__(**kwds)


