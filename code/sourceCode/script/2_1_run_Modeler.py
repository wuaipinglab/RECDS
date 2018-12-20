#!/gluster/Bioshare/anaconda2/bin/python
"""
./2_1_run_Modeler.py seqName
"""
import sys
sys.path.insert(1,'/gluster/home/liquan/tools/python')
sys.path.append('/gluster/home/liquan/tools/MODELLER/PDB/')
sys.path.append('/gluster/home/liquan/tools/MODELLER/CommandFunctions')
from pdbchains import PDBchains
from commandfunctions import Execute
import preparation

runModeller = '/gluster/home/liquan/tools/MODELLER/modeller-9.13/runModeller.py'
fastaDir = "/gluster/home/liquan/antigenicTransition/myProgram/database/seq"
outDir = '/gluster/home/liquan/antigenicTransition/myProgram/structure/Modeler'

if __name__=='__main__':
    pdbNameFile = open(sys.argv[1],'r')
    fastalist = pdbNameFile.read().strip().split('\n')
    fastaNum = len(fastalist)

    template = PDBchains()
    template.Read("/gluster/home/liquan/antigenicTransition/myProgram/structure/native/HA/2viu.pdb")

    for i in range(fastaNum):
        fastaName = fastalist[i]
        workName = fastaName+"A-"+fastaName+"B"
        workDir = preparation.check_slash(outDir) + workName 

        isExist = preparation.mkdir(workDir)
        if isExist:
            newSeqFile = preparation.check_slash(workDir)+fastaName+'A.seq'
            head, query = preparation.read_write_seq(preparation.check_slash(fastaDir)+fastaName,newSeqFile)
            template.Write(preparation.check_slash(workDir)+'Template.pdb')
            template.WriteFasta(preparation.check_slash(workDir)+fastaName+'B.seq',[template.sequence[1]],[fastaName+'B'])
            print 'running start'
            cmd = [runModeller,'-tDir',workName,'-iDir', preparation.check_slash(outDir),'-oDir',preparation.check_slash(workDir)+'ModelerTest','-templ','Template.pdb']
            print cmd
            Execute(cmd,raiseStdError = True)
            print 'running end'
            modeler = PDBchains()
            modeler.Read(preparation.check_slash(workDir)+'ModelerTest/'+fastaName[0:4]+".B99990001.pdb")
            modeler.Write(preparation.check_slash(workDir)+fastaName+'.pdb',specific_chains =[0],write_details =True)

