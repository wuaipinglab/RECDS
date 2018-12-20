#!/usr/bin/env python


import sys,os,re
sys.path.insert(1,'/gluster/home/liquan/tools/python')
import scSef,preparation

AMBERHOME="/gluster/home/liquan/tools/amber16"
os.environ['AMBERHOME'] = AMBERHOME
path = os.environ['PATH']
path = path+":"+AMBERHOME+"/bin"
os.environ['PATH'] = path;
if 'PYTHONPATH' in os.environ:
    os.environ['PYTHONPATH']=os.environ['PYTHONPATH']+":"+AMBERHOME+"/lib/python2.7/site-packages"
else:
    os.environ['PYTHONPATH']=AMBERHOME+"/lib/python2.7/site-packages"

if 'LD_LIBRARY_PATH' in os.environ:
    os.environ['LD_LIBRARY_PATH']=os.environ['LD_LIBRARY_PATH']+":"+AMBERHOME+"/lib"
else:
    os.environ['LD_LIBRARY_PATH']=AMBERHOME+"/lib"

"""
6_1_run_Amber.py 10000 /gluster/home/liquan/antigenicTransition/myProgram/H1N1/myProgram/features/amber/10000 /gluster/home/liquan/antigenicTransition/myProgram/H1N1/myProgram/database/Modeler
"""
#### USER INPUT
AMBERHOME="/gluster/home/liquan/tools/amber16"
pdb          =sys.argv[1]   #modeler model name
OUTPUT        =sys.argv[2]
pdbLibrary    = sys.argv[3]   #modeler model dir
pdb4amber =AMBERHOME+"/bin/pdb4amber"
tleap =AMBERHOME+"/bin/tleap"
sander =AMBERHOME+"/bin/sander"
MMPBSA =AMBERHOME+"/bin/MMPBSA.py"
## PROGRAMS and files;
preparation.mkdir(OUTPUT)


##make error file
err = open(OUTPUT+"/Amber_scwrl_err",'w')
os.chdir(OUTPUT)


pdbFile = pdb+".pdb"
print "*******************************************************"
print pdbFile
print pdb4amber+" -i "+pdbLibrary+"/"+pdbFile+" -o propdb1.pdb --nohyd --dry"
os.system(pdb4amber+" -i "+pdbLibrary+"/"+pdbFile+" -o propdb1.pdb --nohyd --dry")

pdb1 = open("propdb1.pdb",'r')
pdb2 = open("propdb.pdb",'w')
for line in pdb1.readlines():
    if re.match('          H',line):
        continue
    pdb2.write(line)
pdb1.close()
pdb2.close()
    
        
#-------------------------------------------->
tleap_in ="logFile leap.log\nsource leaprc.protein.ff14SB\nprot=loadpdb propdb.pdb\n\n##-- make SYSTEM one\nsaveamberparm prot pdb.prmtop pdb.inpcrd\nsavepdb prot pdb.pdb\n\nquit\n"

tleap_file = open("tleap.in",'w')
tleap_file.write(tleap_in)
tleap_file.close()
print tleap+" -f tleap.in"
os.system(tleap+" -f tleap.in")


#-------------------------------------------->
min_gbsa_sander = "Stage 1 - minimization of protein for GBSA calculation\n&cntrl\n    imin=1, maxcyc=50, ncyc=50, gbsa=1,\n    cut=16., igb=2, ntb=0, saltcon=0.1\n    ntr=0/\n"
min_file = open("min_gbsa.sander",'w')
min_file.write(min_gbsa_sander)
min_file.close()
print sander+" -O -i min_gbsa.sander -p pdb.prmtop -c pdb.inpcrd -r pdb.rst -o pdb.out -inf pdb.inf"
os.system(sander+" -O -i min_gbsa.sander -p pdb.prmtop -c pdb.inpcrd -r pdb.rst -o pdb.out -inf pdb.inf")



#-------------------------------------------->
mmpbsa_in = "MMPBSA input file for running per-residue decomp\n&general\n    startframe=1, endframe=1, interval=1,\n    keep_files=0, debug_printlevel=2\n/\n&gb\n    igb=5, saltcon=0.1\n/\n&decomp\n    idecomp=1, csv_format=0,\n/"
mmpbsa_file=open("mmpbsa.in",'w')
mmpbsa_file.write(mmpbsa_in)
mmpbsa_file.close()
print MMPBSA+" -O -i mmpbsa.in -cp pdb.prmtop -y pdb.rst" 
os.system(MMPBSA+" -O -i mmpbsa.in -cp pdb.prmtop -y pdb.rst")

#-------------------------------------------->
if not os.path.exists("FINAL_DECOMP_MMPBSA.dat"):
    err.write(pdb+"\n")

