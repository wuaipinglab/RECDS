#!/usr/bin/python

"""
./5_1_runHSE.py ../dData
"""
import sys,os,shutil
sys.path.insert(1,"~/.local/lib/python2.7/site-packages")
#from Bio import
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import HSExposure

sys.path.insert(1,'/gluster/home/liquan/tools/python')
import preparation



if __name__ == '__main__':
	
	pdbNameFile = open(sys.argv[1],'r')

    outDir = "/gluster/home/liquan/antigenicTransition/myProgram/database/features/HSE"

	modelDir = '/gluster/home/liquan/antigenicTransition/myProgram/database/Modeler'

	fastalist = pdbNameFile.read().strip().split('\n')
	fastaNum = len(fastalist)
	for i in range(fastaNum):
	        fasta = fastalist[i].strip().split()[1]
                workDir = preparation.check_slash(outDir)+fasta
                preparation.mkdir(workDir)

                parser = PDBParser(PERMISSIVE=1)
                structure_id = fasta
                filename = preparation.check_slash(modelDir)+fasta+'.pdb'
                structure = parser.get_structure(structure_id, filename)
                model = structure[0]
                
                #Class to calculate HSE based on the approximate CA-CB vectors.
                hse_ca = HSExposure.HSExposureCA(model)
                # Calculate HSEbeta
                #Class to calculate HSE based on the real CA-CB vectors
                hse_cb = HSExposure.HSExposureCB(model)
                # Calculate classical coordination number
                exp_fs = HSExposure.ExposureCN(model)
                # Print HSEalpha for a residue
                chain_id = 'A'
                hse_ca_file = open(preparation.check_slash(workDir)+"hse_ca.txt",'w')
                hse_cb_file = open(preparation.check_slash(workDir)+"hse_cb.txt",'w')
                exp_fs_file = open(preparation.check_slash(workDir)+"exp_fs.txt",'w')
                for i in range(1,len(list(model.get_residues()))+1):
                    res_id = (' ',i,' ')
                    if hse_ca.has_key((chain_id, res_id)):
                        """
                        @param hse_up_key: key used to store HSEup in the entity.xtra attribute
                        @param hse_down_key: key used to store HSEdown in the entity.xtra attribute
                        @param angle_key: key used to store the angle between CA-CB and CA-pCB in the entity.xtra attribute
                        """
                        hse_u, hse_d, angle = hse_ca[(chain_id,res_id)]
                        hse_ca_file.write(str(i)+'\t'+str(hse_u)+'\t'+str(hse_d)+'\t'+str(angle)+'\n')
                    else:
                        hse_ca_file.write(str(i)+'\n')

                    if hse_cb.has_key((chain_id, res_id)):
                        hse_u, hse_d, angle = hse_cb[(chain_id,res_id)]
                        hse_cb_file.write(str(i)+'\t'+str(hse_u)+'\t'+str(hse_d)+'\t'+str(angle)+'\n')
                    else:
                        hse_cb_file.write(str(i)+'\n')

                    if exp_fs.has_key((chain_id, res_id)):
                        fs_u, fs_d, angle = hse_cb[(chain_id,res_id)]
                        exp_fs_file.write(str(i)+'\t'+str(fs_u)+'\t'+str(fs_d)+'\t'+str(angle)+'\n')
                    else:
                        exp_fs_file.write(str(i)+'\n')
        hse_ca_file.close()
        hse_cb_file.close()
        exp_fs_file.close()

