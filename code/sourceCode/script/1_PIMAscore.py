#!/usr/bin/env python
"""
 ./script/1_PIMAscore.py ./dData ./features/dPIMA
"""

import sys

PIMA = [
        [1,6,6,6,6,6,6,3,6,6,6,6,6,6,4,4,4,6,6,6],
        [6,1,5,5,6,4,4,6,6,6,6,3,6,6,6,5,5,6,6,6],
        [6,5,1,3,6,5,4,6,6,6,6,5,6,6,6,5,5,6,6,6],
        [6,5,3,1,6,5,3,6,6,6,6,5,6,6,6,5,5,6,6,6],
        [6,6,6,6,1,6,6,6,6,5,5,6,5,5,6,6,6,5,5,5],
        [6,4,5,5,6,1,3,6,6,6,6,4,6,6,6,5,5,6,6,6],
        [6,4,4,3,6,3,1,6,6,6,6,4,6,6,6,5,5,6,6,6],
        [3,6,6,6,6,6,6,1,6,6,6,6,6,6,4,4,4,6,6,6],
        [6,6,6,6,6,6,6,6,1,6,6,6,6,4,6,6,6,4,4,6],
        [6,6,6,6,5,6,6,6,6,1,4,6,4,5,6,6,6,5,5,3],
        [6,6,6,6,5,6,6,6,6,4,1,6,3,5,6,6,6,5,5,4],
        [6,3,5,5,6,4,4,6,6,6,6,1,6,6,6,5,5,6,6,6],
        [6,6,6,6,5,6,6,6,6,4,3,6,1,5,6,6,6,5,5,4],
        [6,6,6,6,5,6,6,6,4,5,5,6,5,1,6,6,6,3,3,5],
        [4,6,6,6,6,6,6,4,6,6,6,6,6,6,1,4,4,6,6,6],
        [4,5,5,5,6,5,5,4,6,6,6,5,6,6,4,1,3,6,6,6],
        [4,5,5,5,6,5,5,4,6,6,6,5,6,6,4,3,1,6,6,6],
        [6,6,6,6,5,6,6,6,4,5,5,6,5,3,6,6,6,1,3,5],
        [6,6,6,6,5,6,6,6,4,5,5,6,5,3,6,6,6,3,1,5],
        [6,6,6,6,5,6,6,6,6,3,4,6,4,5,6,6,6,5,5,1]
        ]

AA = {
        'A':0,
        'R':1,
        'N':2,
        'D':3,
        'C':4,
        'Q':5,
        'E':6,
        'G':7,
        'H':8,
        'I':9,
        'L':10,
        'K':11,
        'M':12,
        'F':13,
        'P':14,
        'S':15,
        'T':16,
        'W':17,
        'Y':18,
        'V':19
        }

featureList = {}
fin = open(sys.argv[1],"r")
for line in fin.readlines():
    iv1 = line.strip().split()[1]
    iv2 = line.strip().split()[2]
    pairNum = int(line.strip().split()[0])

    seqin = open("/gluster/home/liquan/antigenicTransition/myProgram/database/seq/"+iv1,'r')
    seq1 = seqin.readlines()[1].strip()
    seq1 = list(seq1)
    seqin.close()

    seqin = open("/gluster/home/liquan/antigenicTransition/myProgram/database/seq/"+iv2,'r')
    seq2 = seqin.readlines()[1].strip()
    seq2 = list(seq2)
    seqin.close()

    temp = []
    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            temp.append(0)
        else:
            print i, iv1, seq1[i], iv2, seq2[i]
            print AA[seq1[i]], AA[seq2[i]]
            temp.append(PIMA[AA[seq1[i]]][AA[seq2[i]]])
    featureList[pairNum] = temp
fin.close()


outfile = open(sys.argv[2],'w')
for key in featureList:
    outfile.write(str(key)+"\t"+" ".join(map(str,featureList[key]))+"\n")
outfile.close()



