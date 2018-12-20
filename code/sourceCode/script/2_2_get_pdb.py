#!/usr/bin/env python
import os,sys,shutil
"""
./2_2_get_pdb.py seqName ./Modeler ~/antigenicTransition/myProgram/database/Modeler
"""
targets_f = open(sys.argv[1],'r')

fastalist = targets_f.read().strip().split('\n')
for line in fastalist:
	target = line.strip() 
	pdbFile = sys.argv[2]+'/'+target+'A-'+target+'B/'+target+'.pdb'
        newpdbFile = sys.argv[3]+'/'+target+".pdb"
	if os.path.isfile(pdbFile):
		shutil.copy(pdbFile,newpdbFile)
	else:
		print target 
targets_f.close()
    
    
