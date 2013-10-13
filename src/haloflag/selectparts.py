#!/usr/bin/python
import sys
maxmass = 1.0e6
minmass = 0.0

infile = sys.argv[1];
oufile = sys.argv[2];
print "Input File = ", infile
print "Output File = ", oufile

f = open(infile)
maxmass = 1.0e6
minmass = 0.0

s = f.readline()
haloids = [];
for line in f:
    mass = float(line.split()[3])
    hid = int(line.split()[0])
    if(mass <= maxmass and mass >= minmass):
        haloids.append(hid)
                                        
f.close()
ofile = open(oufile, 'w')
ofile.write("%d\n"%len(haloids));
for i in haloids:
    ofile.write("%d\n"%i)
ofile.close
