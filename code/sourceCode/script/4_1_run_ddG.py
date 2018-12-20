#!/usr/bin/env python


import sys
sys.path.insert(1,'/gluster/home/liquan/tools/python')
import scSef,preparation

./4_1_run_ddG.py dData dSEF

featureList = {}
fin = open(sys.argv[1],"r")
outDir = "/gluster/home/liquan/antigenicTransition/myProgram/database/features/SEF/" + sys.argv[2] 

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

    model1 = "/gluster/home/liquan/antigenicTransition/myProgram/database/Modeler/"+iv1+".pdb"
    model2 = "/gluster/home/liquan/antigenicTransition/myProgram/database/Modeler/"+iv2+".pdb"
    sefOut = preparation.check_slash(outDir)+str(pairNum)
    preparation.mkdir(sefOut)
    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            continue
        else:
            pos = str(i+1) 
            sefFile1 = scSef.run_sef(model1, seq1[i]+pos+seq2[i],sefOut)
            sefFile2 = scSef.run_sef(model2, seq2[i]+pos+seq1[i],sefOut)
fin.close()

