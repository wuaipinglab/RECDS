#!/usr/bin/python
import sys,os,re
import numpy as np
sys.path.insert(1,'/gluster/home/liquan/tools/python')
import scSef,preparation



"""
./4_2_get_dssp.py ../dData ../features/dDssp
"""

ddGDir = "/gluster/home/liquan/antigenicTransition/myProgram/database/features/SEF"

to1letter = {
        'ALA':'A',
        'ARG':'R',
        'ASN':'N',
        'ASP':'D',
        'CYS':'C',
        'GLU':'E',
        'GLN':'Q',
        'GLY':'G',
        'HIS':'H',
        'ILE':'I',
        'LEU':'L',
        'LYS':'K',
        'MET':'M',
        'PHE':'F',
        'PRO':'P',
        'SER':'S',
        'THR':'T',
        'TRP':'W',
        'TYR':'Y',
        'VAL':'V'}

def ddGSef(infile,aa2):
    #ddG
    if not os.path.isfile(infile):
        print 'Error: '+infile+' is defective.'   
        return None, None
    else:
        f = open(infile, "r")
        lines = f.read().strip().split('\n')
        f.close()
        if len(lines) <20:
            print "Warning: mutation"+' can not be found in '+ infile+'.'
            return None, None
    ddG = {}
    for i in range(len(lines)):
        cells = lines[i].strip().split()
        ddG[to1letter[cells[0]]]=float(cells[1])
    return ddG[aa2]         

def Extrating(iv1, iv2, ddG):
    seqDir = "/gluster/home/liquan/antigenicTransition/myProgram/database/seq"
    ##########dssp
    seqin = open(preparation.check_slash(seqDir)+iv1,'r')
    seq1 = seqin.readlines()[1].strip()
    seq1 = list(seq1)
    seqin.close()

    seqin = open(preparation.check_slash(seqDir)+iv2,'r')
    seq2 = seqin.readlines()[1].strip()
    seq2 = list(seq2)
    seqin.close()

    ddG1List=[]
    ddG2List=[]
    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            ddG1List.append(0)
            ddG2List.append(0)
        else:
            pos = str(i+1)
            ddG1 = ddGSef(preparation.check_slash(ddG)+seq1[i]+pos+seq2[i]+".ddG",seq2[i])
            ddG2 = ddGSef(preparation.check_slash(ddG)+seq2[i]+pos+seq1[i]+".ddG",seq1[i])
            ddG1List.append(ddG1)
            ddG2List.append(ddG2)
    return np.array(ddG1List), np.array(ddG2List)
	
if __name__ == '__main__':
    pdbNameFile = open(sys.argv[1],'r')
    fin = open(sys.argv[1],"r")
    featureddG1List = {}
    featureddG2List = {}
    featureddGList = {}
    for line in fin.readlines():
        iv1 = line.strip().split()[1] 
        iv2 = line.strip().split()[2]
        pairNum = int(line.strip().split()[0])

        ddG = preparation.check_slash(ddGDir)+sys.argv[3]+"/"+str(pairNum)
        featureddG1List[pairNum],featureddG2List[pairNum] = Extrating(iv1,iv2,ddG)	
        featureddGList[pairNum] = (np.fabs(featureddG1List[pairNum])+np.fabs(featureddG2List[pairNum]))/2	
    fin.close()

    outfile = open(sys.argv[2],'w')
    for key in featureddGList:
        outfile.write(str(key)+"\t"+" ".join(map(str,featureddGList[key]))+"\n")
    outfile.close()
