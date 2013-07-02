#!/usr/bin/python
import sys, re, os, os.path
try:
    sys.argv[1]
except:
    print "Which file? "
    sys.exit("tellme")
inputFile = sys.argv[1]
if not os.path.isfile(inputFile):
    print inputFile, "doesn't exist"
    sys.exit()
print inputFile
#Parse projPsi's output 
writeEX  = False
writeION = False
writeY   = False
f = open(inputFile, 'r')
lines = f.readlines()
#only opening files for which we have data
for line in lines:
    exState = re.match( "^[1-9]", line)
    kState  = re.match( "^0", line)
    Yl      = re.match( "^Y(\d)0", line)
    if exState:
        writeEX = True
    if kState:
        writeION = True
    if Yl:
        writeY = True
        print "matchY"
        break
if writeEX:
    exFile = open("ex.dat", 'w')
if writeION:
    ionFile = open("ion.dat", 'w')
if writeY:
    RlFile = open("Rl.dat", 'w')
while 1:
    if(lines):
        line = lines.pop(0)
        exState = re.match( "^[1-9]", line)
        kState = re.match( "^0", line)
        Yl = re.match( "^Y(\d)0", line)
        if exState:
            words = line.split()
            exFile.write(' '.join(words) + '\n')
        if kState:
            words = line.split()
            words.pop(0)
            for word in words:
                if( word == 'took'):
                    break
                ionFile.write(word + "\t")
            ionFile.write("\n")
        if Yl:
            words  = line.split()
            l  =  Yl.group(1)
            Pl =  lines.pop(0)[0:-1]
            RlFile.write( l + "\t" + Pl + "\t" + " ".join(words[1:]) + "\n" )
    else:
        break
f.close()
if writeEX:
    exFile.close()
if writeION:
    ionFile.close()
if writeY:
    RlFile.close()

#Sort STEP, atomic time, and walltime
f = open("time.dat", 'r')
lines = f.readlines()
legacyTime = 0
lastTime  = 0
lastStep  = 0
fout = open("t.dat", "w")
for line in lines:
    time = line.split()
    step  = float(time[0]) # time step number
    qTime = time[1]        # quantum time
    wTime = float(time[2]) # wall time
    if step <= lastStep:   # Don't allow repeat time steps
        continue
    if wTime < lastTime:   # add to legacy time
        legacyTime += lastTime
    totalTime = wTime + legacyTime
    fout.write(time[0] + "\t" + time[1] + "\t" + str(totalTime) + "\n")
    lastStep = step
    lastTime = wTime
f.close()
tFile = open("tMAX.dat", 'w')
tFile.write( qTime )
tFile.close()

if os.path.isfile('input'):
    os.system('cp input input.dat')
if os.path.isfile('input2'):
    os.system("cp input2 input2.dat")
#make dr.dat
if( writeY ):
    thisDIR = os.getcwd()
    inputFile = thisDIR  + '/input'
    input2File = thisDIR + '/input2'
    if os.path.isfile(inputFile):
        f = open("input", 'r')
        lines = f.readlines()
        for line in lines:
            if line:
                word = line.split()
                if word[0] == 'L':
                    L = float(word[1])
        f.close()
        if os.path.isfile(input2File):
            f = open("input2", 'r')
            lines = f.readlines()
            for line in lines:
                if line:
                    word = line.split()
                    if word[0] == 'nGrid':
                        n = int(word[1])
            f.close()
            f = open("dr.dat", 'w')
            f.write( str(L/n) )
            f.close()
