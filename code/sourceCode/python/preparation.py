#!/usr/bin/env python
from subprocess import Popen, PIPE
import sys,re,os,datetime,random,urllib,urllib2
sys.path.insert(1,'/nfs/amino-home/liquan/CancerDriver/myprogram/python')
import userConfig

def get_uniqueNum():
	nowTime = datetime.datetime.now().strftime('%Y%m%d%H%M%S')
	randomNum = random.randint(0,100)
	if randomNum <= 10:
		randomNum = str(0)+str(randomNum)
	uniqueNum = str(nowTime)+str(randomNum)
	print uniqueNum
	return uniqueNum

def mkdir(path):
    path=path.strip()
    path=path.rstrip('/')
    isExists=os.path.exists(path)
 
    if not isExists:
        print path+' ....'
        os.makedirs(path)
        return True
    else:
        print path+' .....'
        return False

def check_slash(route):
	if not re.search('/$',route.strip()):
		route = route.strip()+'/'
	return route.strip()

def delete_slash(route):
	if re.search('/$',route.strip()):
		route = re.sub('/$','',route.strip())
	return route.strip()

def read_write_seq(seqFile, outFile):
	if not os.path.isfile(seqFile):
		print 'ERROR: '+seqFile+' does not exist!'
	nCh = 0
	newSeq = ''
	head = ''

	fin = open(seqFile,'r')
	fline = fin.readline()
	if not re.match('^>',fline):
		nCh = 1
	while fline:
		if re.match('^>',fline):
			nCh += 1
			if (nCh > 1):
				print 'LOG: submitted protein sequence contained multiple chains. It will use only the first protein chain from the provided file'
				break
			head = fline
		else:
			fline = re.sub('\s+','',fline)
			for aa in fline:
				if re.match('^[ARNDCEQGHILKMFPSTWYV]{1}$',aa):
					newSeq += aa
				else:
					print 'LOG: '+aa+' is not a standard amino acid. It was removed automatically.'
		fline = fin.readline()
	fout = open(outFile,'w')
	if head =='':
		head = '>query\n'
	fout.write(head)
	fout.write(newSeq)
	return re.sub('^>','',head).strip(), newSeq

	
def run_cmd(cmd,mod=0):
	handle = Popen(cmd,shell = True, stdout=PIPE, stderr=PIPE)
	stdout,stderr = handle.communicate()
	if stderr != '':
		print '*******'+stderr
		if mod == 0:
			sys.exit()
	return stdout

def rmfile(myfile):
	if os.path.isfile(myfile):
		os.remove(myfile)


def download(url,local):
	urllib.urlretrieve(url,local)

def readUrl(url):
	content = ''
	try:
		cto = urllib2.urlopen(url)
		content = cto.read()
	except Exception, e:
		print url
		print e
	return content
