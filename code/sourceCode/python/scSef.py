#!/usr/bin/env python
import sys,re,os
sys.path.insert(1,'/gluster/home/liquan/tools/python')
import preparation,userConfig

sefddG = 'thermo/SingleMutationScanning'

def run_sef(model,mutation,outdir):
	outdir = preparation.delete_slash(outdir)	
	primaryDir = os.getcwd()
	os.chdir(userConfig.sef)
	#print os.getcwd()
	cmd = 'java '+sefddG+' '+model+' '+mutation+' '+outdir
	#print cmd
	ddG = preparation.run_cmd(cmd) 
	
	outFile = preparation.check_slash(outdir)+mutation+'.ddG'
	fout = open(outFile,'w')
	fout.write(ddG)
	os.chdir(primaryDir)
	#print os.getcwd()
	return outFile
