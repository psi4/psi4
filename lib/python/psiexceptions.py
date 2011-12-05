import PsiMod

class PsiException(Exception): pass
class ValidationError(PsiException):
    def __init__(self, msg):
        PsiException.__init__(self, msg)
        self.msg = msg
        PsiMod.print_out('\nPsiException: %s\n\n' % (msg))

