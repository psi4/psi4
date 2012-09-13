"""Module with non-generic exceptions classes."""


class QcdbException(Exception):
    """Error class for QCDB."""
    pass


class FeatureNotImplemented(QcdbException):
    """Error called for functions defined but not yet implemented.
    Also for functions defined that will never be implemented.

    """
    def __init__(self, msg):
        QcdbException.__init__(self, msg)
        self.msg = msg
        print('\nQcdbException: Feature %s is not yet implemented.\n\n' % (msg))


class ValidationError(QcdbException):
    """Error called for problems with syntax input file. Prints
    error message *msg* to standard output stream.

    """
    def __init__(self, msg):
        QcdbException.__init__(self, msg)
        self.msg = msg
        print('\nQcdbException: %s\n\n' % (msg))


class IncompleteAtomError(QcdbException):
    """Error raised when not all variables in an atom specification
    have been defined at compute time. May be a temporary situation
    so message not printed but appears as traceback when error persists.

    """
    def __init__(self, msg):
        QcdbException.__init__(self, msg)
        self.msg = msg
