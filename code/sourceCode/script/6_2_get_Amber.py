#!/usr/bin/python
import sys,os,re
import numpy as np


"""
illustrating the code for receiving the inpu
./6_2_get_Amber.py ../dData ../features/dAmber
./6_2_get_Amber.py ../dData ../features/dAmber
"""

amberfile = "/gluster/home/liquan/antigenicTransition/myProgram/H1N1/myProgram/features/amber"

def amber(ivfile,pos_start, pos_end):
	if not os.path.isfile(ivfile):
            print ivfile
            return None, None
	f = open(ivfile, "r")
	lines = f.readlines()
	f.close()

        internal = []
        vanDerWaals = []
        electrostatic = []
        polarSolvation = []
        nonPolarSolv = []
        total = []
        if len(lines)<9:
            print ivfile
            return None, None
	for i in range(8,len(lines)):
                cells = lines[i].strip().split()
                if len(cells)<26:
                    continue
		residue = cells[1]
		if residue:
			residue = int(residue)
		else:
			continue
		if residue in range(pos_start, pos_end+1):
                        internal.append(float(cells[3]))
                        vanDerWaals.append(float(cells[7]))
                        electrostatic.append(float(cells[11]))
                        polarSolvation.append(float(cells[15]))
                        nonPolarSolv.append(float(cells[19]))
                        total.append(float(cells[23]))

	return np.array(internal),np.array(vanDerWaals),np.array(electrostatic),np.array(polarSolvation),np.array(nonPolarSolv),np.array(total)



def Extrating(iv1, iv2):
	internal1, vanDerWaals1, electrostatic1, polarSolvation1, nonPolarSolv1, total1 = amber(amberfile+"/"+iv1+"/FINAL_DECOMP_MMPBSA.dat",1,329)
        internal2, vanDerWaals2, electrostatic2, polarSolvation2, nonPolarSolv2, total2 = amber(amberfile+"/"+iv2+"/FINAL_DECOMP_MMPBSA.dat",1,329)
        internal = np.fabs(internal1 - internal2)
        vanDerWaals = np.fabs(vanDerWaals1- vanDerWaals2)
        electrostatic  = np.fabs(electrostatic1- electrostatic2)
        polarSolvation  = np.fabs(polarSolvation1- polarSolvation2)
        nonPolarSolv  = np.fabs(nonPolarSolv1 - nonPolarSolv2)
        total  = np.fabs(total1 - total2)
        
	return internal, vanDerWaals, electrostatic, polarSolvation, nonPolarSolv, total
	
if __name__ == '__main__':
        fin = open(sys.argv[1],"r")
        internal = {}
        vanDerWaals = {}
        electrostatic = {}
        polarSolvation = {}
        nonPolarSolv = {}
        total = {}
        for line in fin.readlines():
            iv1 = line.strip().split()[1] 
            iv2 = line.strip().split()[2]
            pairNum = int(line.strip().split()[0])

            internal[pairNum], vanDerWaals[pairNum], electrostatic[pairNum], polarSolvation[pairNum], nonPolarSolv[pairNum], total[pairNum] = Extrating(iv1,iv2)	
        fin.close()

        outfile = open(sys.argv[2]+"Internal",'w')
        for key in internal:
            outfile.write(str(key)+"\t"+" ".join(map(str,internal[key]))+"\n")
        outfile.close()

        outfile = open(sys.argv[2]+"VanDerWaals",'w')
        for key in vanDerWaals:
            outfile.write(str(key)+"\t"+" ".join(map(str,vanDerWaals[key]))+"\n")
        outfile.close()

        outfile = open(sys.argv[2]+"electrostatic",'w')
        for key in electrostatic:
            outfile.write(str(key)+"\t"+" ".join(map(str,electrostatic[key]))+"\n")
        outfile.close()

        outfile = open(sys.argv[2]+"PolarSolvation",'w')
        for key in polarSolvation:
            outfile.write(str(key)+"\t"+" ".join(map(str,polarSolvation[key]))+"\n")
        outfile.close()

        outfile = open(sys.argv[2]+"NonPolarSolv",'w')
        for key in nonPolarSolv:
            outfile.write(str(key)+"\t"+" ".join(map(str,nonPolarSolv[key]))+"\n")
        outfile.close()

        outfile = open(sys.argv[2]+"total",'w')
        for key in total:
            outfile.write(str(key)+"\t"+" ".join(map(str,total[key]))+"\n")
        outfile.close()



