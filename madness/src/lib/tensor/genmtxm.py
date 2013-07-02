def doit(transpose):
    t = ""
    if transpose: t = "T"
    for i in range(2,28,2):
        if i == 12:
            print "#ifdef X86_64"

        print "MTXM_ENTRY(%smTxm%d)" % (t,i)

        jloop = ".%sJLOOP%d" % (t,i)
        print "%s:" % jloop

        for j in range(0,i/2):
            print "ZERO(C%02.2d)" % j

        kloop = ".%sKLOOP%d" % (t,i)
        print "%s:" % kloop
        print "LOADBKJ"

        for j in range(0,i/2):
            print "ABC(%d,C%02.2d)" % (j*16,j)

        print ("INCA")

        for j in range(0,i/2):
            print "ABC2(%d,C%02.2d)" % (j*16,j)

        print "INCA"
        print "sub     $2, NK"
        print "jnz     %s" % kloop

        print "mov     C, B"
        for j in range(0,i/2):
            if transpose:
                print "STORET(C%02.2d)" % j
            else:
                print "STORE(C%02.2d)" % j

        if transpose:
            print "INCCT"
        else:
            print "INCC"
        print "NEXTJ"
        print "jnz     %s" % jloop
        print "RETURN"
    print "#endif"

doit(0)
doit(1)
