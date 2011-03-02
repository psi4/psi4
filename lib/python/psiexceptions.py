import PsiMod

class PsiException: pass
class ValueNotSet (PsiException): pass
class RowAlignmentError(PsiException):

    def __init__(self, l1name, l1, l2name, l2):
        msg = "Rows %s and %s not aligned. Length %d != %d" % (l1name, l2name, l1, l2)
        PsiException.__init__(self, msg)

