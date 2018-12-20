#!/usr/bin/python
import sys,os,re
import numpy as np


"""
./3_2_get_dssp.py ../dData ../features/dDssp
"""

dsspfile = "/gluster/home/liquan/antigenicTransition/myProgram/database/features/dssp_modeler"

extendedASA = {
	'A':110.2,
	'C':140.4,
	'D':144.1,
	'E':174.7,
	'F':200.7,
	'G':78.7,
	'H':181.9,
	'I':185.0,
	'K':205.7,
	'L':183.1,
	'M':200.1,
	'N':146.4,
	'P':141.9,
	'Q':178.6,
	'R':229.0,
	'S':117.2,
	'T':138.7,
	'V':153.7,
	'W':240.5,
	'Y':213.7};


def Dssp(ivfile,pos_start, pos_end):
	if not os.path.isfile(ivfile):
            print ivfile
	f = open(ivfile, "r")
	content = f.read()
	lines = content.strip().split("\n")
	f.close()

	start = -1
	for i in range(0,len(lines)):
		if re.search('#  RESIDUE AA STRUCTURE',lines[i]):
			start = i+1
			break	
	if start == -1:
		return False, False, False, False, False, False
        ss = []
        acc = []
        acc2 = []
        #kappa = []
        #phi = []
        #psi = []
	for i in range(start,len(lines)):
		residue = lines[i][5:10].strip()
		if residue:
			residue = int(residue)
		else:
			continue
		if residue in range(pos_start, pos_end+1):
			#tmp_ss = lines[i][16:17]
			#if tmp_ss == ' ':
			#	tmp_ss = "C"
                        #ss.append(tmp_ss)
                        tmp_acc = float(lines[i][35:38].strip()) 
			acc.append(tmp_acc)
                        tmp_acc2 = tmp_acc/extendedASA[str(lines[i][13:14])]
			if tmp_acc2 > 1.0:
				#print tmp_acc2
				tmp_acc2 =1.0
                        acc2.append(tmp_acc2)
	return np.array(acc2)



def Extrating(iv1, iv2):
	##########dssp
	iv1_acc2 = Dssp(dsspfile+"/"+iv1+".dssp",1,329)
        iv2_acc2 = Dssp(dsspfile+"/"+iv2+".dssp",1,329)
        features = np.fabs(iv1_acc2 - iv2_acc2)
	return features
	
if __name__ == '__main__':
	pdbNameFile = open(sys.argv[1],'r')
        fin = open(sys.argv[1],"r")
        featureList = {}
        for line in fin.readlines():
            iv1 = line.strip().split()[1] 
            iv2 = line.strip().split()[2]
            pairNum = int(line.strip().split()[0])

            featureList[pairNum] = Extrating(iv1,iv2)	
        fin.close()

        outfile = open(sys.argv[2],'w')
        for key in featureList:
            outfile.write(str(key)+"\t"+" ".join(map(str,featureList[key]))+"\n")
        outfile.close()
