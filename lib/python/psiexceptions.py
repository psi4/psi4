"""Module with non-generic exceptions classes."""
import PsiMod


class PsiException(Exception):
    """Error class for Psi."""
    pass


class ValidationError(PsiException):
    """Error called for problems with the input file. Prints
    error message *msg* to standard output stream and output file.

    """
    def __init__(self, msg):
        PsiException.__init__(self, msg)
        self.msg = msg
        PsiMod.print_out('\nPsiException: %s\n\n' % (msg))
