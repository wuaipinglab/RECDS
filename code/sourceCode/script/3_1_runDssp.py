#!/usr/bin/python
import sys,os,shutil


"""
./3_1_runDssp.py ../seqNameMap
"""

dsspdir = "/gluster/home/liquan/bin/dssp-2.0.4-linux-amd64"

def Rundssp(protein, outfile):
	"""
	rung depth
	"""
        print dsspdir+" "+protein+" >"+outfile
	os.system(dsspdir+" "+protein+" >"+outfile)

if __name__ == '__main__':
	
	pdbNameFile = open(sys.argv[1],'r')
        outDir = '/gluster/home/liquan/antigenicTransition/myProgram/database/features/dssp_modeler'
        if not os.path.exists(outDir):
            os.mkdir(outDir)  
	fastalist = pdbNameFile.read().strip().split('\n')
	fastaNum = len(fastalist)
	for i in range(fastaNum):
		fasta = fastalist[i].strip().split()[1]
		Modeler = '/gluster/home/liquan/antigenicTransition/myProgram/database/Modeler/'+fasta+".pdb"
                outdir = '/gluster/home/liquan/antigenicTransition/myProgram/database/features/dssp_modeler/'+fasta+".dssp"
		Rundssp(Modeler,outdir)

