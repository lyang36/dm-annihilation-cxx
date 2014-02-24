import os
import numpy as np

def getVLData(filename, xc, yc, zc, boxSize, basename="", mask=""):
    cmd = "./VLShow ";
    if(basename != ""):
        cmd = cmd + " -b " + basename;
    cmd = cmd + " -f " + filename;
    cmd = cmd + " -box %10.10f %10.10f %10.10f %10.10f"%(xc, yc, zc, boxSize)
    if(mask != ""):
        cmd = cmd + " -m " + mask;

#print cmd;
    stream = os.popen(cmd);
    data = [];
    for r in stream.readlines():
        floats = [float(x) for x in r.strip().split()]
        data.append(floats)
    return np.array(data)
