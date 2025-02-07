#!/usr/bin/env python
import sys
"""
./12_merge_feature.py dFeatures
"""

def extractFeature(infile):
    print infile
    inf = open(infile,"r")
    lines = inf.readlines()
    inf.close()

    fDict = {}
    for line in lines:
        cells = line.strip().split()
        num = int(cells[0])
        fDict[num] = cells[1:len(cells)]
    return fDict

            
if __name__== '__main__':

    #the dir of all the features derived from different antigenic pairs
    featureDir = "/gluster/home/liquan/antigenicTransition/myProgram/H1N1/myProgram/features/"

    file1 = featureDir+"dPIMA"
    file2 = featureDir+"dDssp"
    file3 = featureDir+"dSEF"
    file4 = featureDir+"dHSEUp"
    file5 = featureDir+"dAmberVanDerWaals"
    file6 = featureDir+"dAmberNonPolarSolv"

    f1 = extractFeature(file1)
    f2 = extractFeature(file2)
    f3 = extractFeature(file3)
    f4 = extractFeature(file4)
    f5 = extractFeature(file5)
    f6 = extractFeature(file6)

    outf = open(sys.argv[1],'w')
    for i in f1.keys():
        line = str(i)+" "
        for j in range(len(f1[i])):
            line += f1[i][j]+" "+f2[i][j]+" "+f3[i][j]+" "+f4[i][j]+" "+f5[i][j]+" "+f6[i][j]+" "
        outf.write(line+"\n")
