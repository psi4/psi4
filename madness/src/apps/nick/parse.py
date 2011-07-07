#!/usr/bin/env python
# Write the name of the data file in the argument list
# Known issues: multiple file arguments get included 
# too many times
import sys, string, re, array, os
try:
        sys.argv[1]
except:
        print "Enter the file as an argument"
        sys.exit()
number  = []
time    = []
field   = []
energy  = []
norm    = []
overlap = []
xdip    = []
ydip    = []
zdip    = []
accel   = []
Wtime   = []
fileList = sys.argv[1:]

for file in fileList:
        print file
        f = open(file)
        lines = f.readlines()
        #Find the starting pattern
        while lines:
                line = lines.pop(0)
                m = re.match( "\s*step",line)
                if m:
                        break
        #get the data
        for line in lines:
                if re.match("\A\s*Active",line):
                        print "Found the End"
                        break
                if re.match("\A\s+\d+",line):
                        word = line.split()
                        number.append( word[-1])
                        time.append(   word[1])
                        field.append(  word[2])
                        energy.append( word[3])
                        norm.append(   word[4])
                        overlap.append(word[5])
                        xdip.append(   word[6])
                        ydip.append(   word[7])
                        zdip.append(   word[8])
                        accel.append(  word[9])
                        Wtime.append(  word[10])
                if not line:
                        break
timeFile   = open("time.dat",      'w')
energyFile = open("energy.dat",    'w')
fieldFile  = open("field.dat",     'w')
normFile   = open("norm.dat",      'w')
overlapFile= open("overlap.dat",   'w')
xdipFile   = open("xdip.dat",      'w')
ydipFile   = open("ydip.dat",      'w')
zdipFile   = open("zdip.dat",      'w')
accelFile  = open("accel.dat",     'w')
WtimeFile  = open("Wtime.dat",     'w')
print "writing files"
for i in range(0,len(number)):
        timeFile.write(number[i]  + "\t" + time[i]   + "\t" + Wtime[i] + "\n")
        energyFile.write(time[i]  + "\t" + energy[i] + "\n")
        fieldFile.write(time[i]   + "\t" + field[i]  + "\n")
        normFile.write(time[i]    + "\t" + norm[i]   + "\n")
        overlapFile.write(time[i] + "\t" + overlap[i]+ "\n")
        xdipFile.write(time[i]    + "\t" + xdip[i]   + "\n")
        ydipFile.write(time[i]    + "\t" + ydip[i]   + "\n")
        zdipFile.write(time[i]    + "\t" + zdip[i]   + "\n")
        accelFile.write(time[i]   + "\t" + accel[i]  + "\n")
        WtimeFile.write(Wtime[i]  + "\n")
os.system("cp input input.dat")
#os.system("gnuplot -persist ~/bin/makePlots.gp")
