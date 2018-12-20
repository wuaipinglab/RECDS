#!/usr/bin/python
import sys,os,re
import numpy as np


"""
./5_2_get_HSE.py ../dData ../features/dHSE
"""

hsefile = "/gluster/home/liquan/antigenicTransition/myProgram/database/features/HSE"

def HSE(ivfile,pos_start, pos_end):
	if not os.path.isfile(ivfile):
            print ivfile
            return None, None
	f = open(ivfile, "r")
	content = f.read()
	lines = content.strip().split("\n")
	f.close()

        up = []
        down = []
	for i in range(0,len(lines)):
                cells = lines[i].strip().split()
		residue = cells[0]
		if residue:
			residue = int(residue)
		else:
			continue
		if residue in range(pos_start, pos_end+1):
                        tmp_up = float(cells[1]) 
                        tmp_down = float(cells[2]) 
			up.append(tmp_up)
			down.append(tmp_down)
	return np.array(up),np.array(down)



def Extrating(iv1, iv2):
	##########dssp
	iv1_up, iv1_down = HSE(hsefile+"/"+iv1+"/hse_cb.txt",1,329)
        iv2_up, iv2_down = HSE(hsefile+"/"+iv2+"/hse_cb.txt",1,329)
        upD = np.fabs(iv1_up - iv2_up)
        downD = np.fabs(iv1_down - iv2_down)
	return upD, downD
	
if __name__ == '__main__':
        fin = open(sys.argv[1],"r")
        featureUpDList = {}
        featureDownDList = {}
        for line in fin.readlines():
            iv1 = line.strip().split()[1] 
            iv2 = line.strip().split()[2]
            pairNum = int(line.strip().split()[0])

            featureUpDList[pairNum],featureDownDList[pairNum] = Extrating(iv1,iv2)	
        fin.close()

        outfile = open(sys.argv[2]+"Up",'w')
        for key in featureUpDList:
            outfile.write(str(key)+"\t"+" ".join(map(str,featureUpDList[key]))+"\n")
        outfile.close()

        outfile = open(sys.argv[2]+"Down",'w')
        for key in featureDownDList:
            outfile.write(str(key)+"\t"+" ".join(map(str,featureDownDList[key]))+"\n")
        outfile.close()
