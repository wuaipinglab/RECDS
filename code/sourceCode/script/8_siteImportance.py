#!/usr/bin/env python
"""
./oneCombTrain3.py ./dFeatures ./sFeatures ./train/allImportance
"""

import sys, numpy, time

sys.path.insert(1,"~/.local/lib/python2.7/site-packages")
from sklearn.model_selection import KFold
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import GradientBoostingClassifier
from itertools import combinations

_FILE_dFEATURE = sys.argv[1]
_FILE_sFEATURE = sys.argv[2]

if __name__=="__main__":

    featureNum = 6
    with open(_FILE_dFEATURE) as f:
        size = len(f.readline().strip().split())
    ddata = numpy.genfromtxt(_FILE_dFEATURE,usecols=xrange(1,size))
    dlable = numpy.array([1]*ddata.shape[0])

    with open(_FILE_sFEATURE) as f:
        size = len(f.readline().strip().split())
    sdata = numpy.genfromtxt(_FILE_sFEATURE,usecols=xrange(1,size))
    slable = numpy.array([0]*sdata.shape[0])

    data = numpy.vstack((ddata,sdata))
    lable = numpy.append(dlable,slable)

    
    iterations = 50
    aveImp = {}
    aveImpAll = []
    for loop in range(iterations):
        pred = GradientBoostingClassifier(learning_rate=0.05, n_estimators=3000, max_depth=6,max_features="sqrt")
        pred.fit(data,lable)
        importances = pred.feature_importances_
        print "**************************************************************"
        print loop, pred.score(data,lable)

        site_imp = {}
        i = 0
        site_num = 1
        while i < importances.size:
            temp_imp = 0.0
            for j in range(featureNum):
                temp_imp+=importances[i]
                i+=1
            site_imp[site_num] = temp_imp
            site_num+=1
    
        if loop ==0:
            aveImp = site_imp
            aveImpAll = importances
        else:
            for key in aveImp:
                aveImp[key]+=site_imp[key]
            aveImpAll += importances
        sortedSite = sorted(site_imp.iteritems(),key=lambda asd:asd[1],reverse=True)

        outf = open(sys.argv[3]+"_V"+str(loop),"w")
        for f in range(1,site_num):
            outf.write(str(f)+"\t"+str(site_imp[f])+"\n")
        outf.close()

    for key in aveImp:
        aveImp[key]=float(aveImp[key])/iterations
    aveSortedSite = sorted(aveImp.iteritems(),key=lambda asd:asd[1],reverse=True)
    print "**************************************************************"
    print "average!!"

    outf = open(sys.argv[3]+"_ave","w")
    for f in range(1,site_num):
        outf.write(str(f)+"\t"+str(aveImp[f])+"\n")
    outf.close()

    aveImpAll = aveImpAll/float(iterations)
    outf = open(sys.argv[3]+"_all_ave","w")
    f=0
    while f <aveImpAll.size:
        outf.write(str(aveImpAll[f])+"\n")
        f += 1
    outf.close()

