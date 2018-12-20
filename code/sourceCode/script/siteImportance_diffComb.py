#!/usr/bin/env python
"""
./bestFeatureImportance.py ./d16Features ./s16Features ./train/siteImportant 5
"""

import sys, numpy, time

sys.path.insert(1,"~/.local/lib/python2.7/site-packages")
from sklearn.model_selection import KFold
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import GradientBoostingClassifier
from itertools import combinations

_FILE_dFEATURE = sys.argv[1]
_FILE_sFEATURE = sys.argv[2]
keySite = [145,155,156,158,159,189,193]

def accuracy(lables,predictions):
    tp = 0
    fp = 0
    tn = 0
    fn = 0
    for i in range(len(lables)):
        if lables[i] == 1:
            if lables[i] == predictions[i]:
                tp += 1
            else:
                fp += 1
        if lables[i] == 0:
            if lables[i] == predictions[i]:
                tn += 1
            else:
                fn += 1
    print tp, fp, tn, fn
    mcc = float(tp*tn-fp*fn)/numpy.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
    acc = float(tp+tn)/(tp+fp+tn+fn)
    return mcc, acc

def topN(sortedList, topNum):
    siteList = []
    rank = {}
    for i in range(topNum):
        siteList.append(sortedList[i][0])
        rank[sortedList[i][0]]=i+1
    matchNum = 0
    matchRank = 0
    for key in keySite:
        if key in siteList:
            matchNum +=1
            matchRank += rank[key]
            print key, rank[key]
    if matchNum == 0:
        return matchNum, matchNum
    return matchNum, float(matchRank)/matchNum

if __name__=="__main__":
    allfeatures = range(1,17)
    featureNum = int(sys.argv[4])
    diffCombs = combinations(allfeatures,featureNum)
    _FOLD = 5

    for diffComb in diffCombs:
        featureList = []
        for i in range(329):
            for j in diffComb:
                featureList.append(16*i+j)
        featureList = tuple(featureList)

        combName = "_".join(map(str,diffComb))
        ddata = numpy.genfromtxt(_FILE_dFEATURE,usecols=featureList)
        dlable = numpy.array([1]*ddata.shape[0])

        sdata = numpy.genfromtxt(_FILE_sFEATURE,usecols=featureList)
        slable = numpy.array([0]*sdata.shape[0])

        data = numpy.vstack((ddata,sdata))
        lable = numpy.append(dlable,slable)

        
        kfold = StratifiedKFold(n_splits=_FOLD, shuffle=True)
        last_predicts = []
        last_lables = []
        last_index = []
        for train, test in kfold.split(data,lable):
            pred = GradientBoostingClassifier(learning_rate=0.05, n_estimators=3000, max_depth=6,max_features="sqrt")
            data_train, data_test, lable_train, lable_test = data[train], data[test], lable[train], lable[test]
            pred.fit(data_train, lable_train)
            pred_test = pred.predict(data_test)
            for j in range(len(pred_test)):
                last_predicts.append(pred_test[j])
                last_lables.append(lable_test[j])
                last_index.append(test[j])
        mcc, acc = accuracy(last_lables, last_predicts)
        print "**************************************************************"
        print combName, mcc, acc 

        pred = GradientBoostingClassifier(learning_rate=0.05, n_estimators=3000, max_depth=6,max_features="sqrt")
        pred.fit(data,lable)
        importances = pred.feature_importances_  
        print combName, pred.score(data,lable)

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

        sortedSite = sorted(site_imp.iteritems(),key=lambda asd:asd[1],reverse=True)
        matchNum, matchRank = topN(sortedSite, 10)
        print "10\t"+str(matchNum)+"\t"+str(matchRank)

        matchNum, matchRank = topN(sortedSite, 20)
        print "20\t"+str(matchNum)+"\t"+str(matchRank)

        outf = open(sys.argv[3]+"_"+combName,"w")
        for f in range(1,site_num):
            outf.write(str(f)+"\t"+str(site_imp[f])+"\n")


    


