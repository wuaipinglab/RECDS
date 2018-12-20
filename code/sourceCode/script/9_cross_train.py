#!/usr/bin/env python
"""
./train.py dFeatures sFeatures crossTrain.out
./train.py /gluster/home/liquan/antigenicTransition/myProgram/script/underFeatures/df12 /gluster/home/liquan/antigenicTransition/myProgram/script/underFeatures/sf12 crossTrainF12.out
"""

import sys, numpy, time

sys.path.insert(1,"~/.local/lib/python2.7/site-packages")
from sklearn.model_selection import KFold
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import GradientBoostingClassifier

_FILE_dFEATURE = sys.argv[1]
_FILE_sFEATURE = sys.argv[2]
_FILE_out = sys.argv[3]

_ITERATION = 50
_FOLD = 5

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




if __name__=="__main__":
    with open(_FILE_dFEATURE) as f:
        size = len(f.readline().strip().split())
    ddata = numpy.genfromtxt(_FILE_dFEATURE,usecols=xrange(1,size))
    dlable = numpy.array([1]*ddata.shape[0])

    with open(_FILE_sFEATURE) as f:
        size = len(f.readline().strip().split())
    sdata = numpy.genfromtxt(_FILE_sFEATURE,usecols=xrange(1,size))
    slable = numpy.array([0]*sdata.shape[0])

    data= numpy.vstack((ddata,sdata))
    lable = numpy.append(dlable,slable)

    mcc_max = 0
    accuracy_max = 0
    best_predicts = []
    best_lables = []
    best_index = []

    for i in range(_ITERATION):
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
        print i, mcc, acc 
        if mcc_max < mcc:
            mcc_max = mcc
            accuracy_max = acc
            best_predicts = last_predicts
            best_lables = last_lables
            best_index = last_index

    print "*********************************"
    print "best", mcc_max, accuracy_max
    outf = open(_FILE_out,"w")
    for i in range(len(best_predicts)):
        outf.write(str(best_index[i])+"\t"+str(best_lables[i])+"\t"+str(best_predicts[i])+"\n")
    outf.close()



    


