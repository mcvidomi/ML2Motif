from sklearn.model_selection import train_test_split
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import LeaveOneOut
from shogun.Classifier import *
from shogun.Evaluation import *
from shogun.Features import *
from shogun.Kernel import *
import pickle
import numpy as np
import copy
import pdb
from poim import *

class SVM:
    '''
    --- use shogun toolbox ---
    feature: weighted degree string kernel
    '''

    def __init__(self,x=None,y=None,kernel_degree=8):
        self.x = x
        self.y = y
        self.kernel_degree = kernel_degree

        if not (self.x is None):
            self.L = len(x[0])
        if not (self.x is None) and type(self.x) == np.ndarray:
            self.x = self.x.tolist()
        if not (self.y is None) and type(self.y) != np.ndarray:
            self.y = np.array(self.y)

    def svm_save(self, file=None):
        if file is None:
            file = "SVM.pkl"
        fobj = open(file, 'wb')
        pickle.dump([self.svm], fobj)
        fobj.close()

    def save(self,file=None):
        if file is None:
            file = "SVM.pkl"
        fobj = open(file, 'wb')
        pickle.dump(self.__dict__, fobj)
        fobj.close()

    def load(self,file=None):
        if file is None:
            file = "SVM.pkl"
        fobj = open(file, 'rb')
        self.__dict__ = pickle.load(fobj)
        fobj.close()

    def svm_load(self, file=None):
        if file is None:
            file = "SVM.pkl"
        fobj = open(file, 'rb')
        self.svm = pickle.load(fobj)
        fobj.close()

    def train(self, C=1):
        if not(self.y is None) and not(type(self.y) is np.ndarray):
            self.y = np.array(self.y)
        '''
        train SVM
        :param C: regularization parameter
        :return: -
        '''
        feats_train = StringCharFeatures(DNA)
        feats_train.set_features(self.x) # ask if x is list
        labels = BinaryLabels(self.y) # ask if y is array
        self.C = C
        self.kernel = WeightedDegreePositionStringKernel(feats_train, feats_train, self.kernel_degree)
        self.svm = LibSVM(self.C, self.kernel, labels)
        self.svm.train()
        self.ma = compute_poims(self.svm, self.kernel, 7, self.L)
        self.save()

    def apply(self,x):
        if x != None and type(x) == np.ndarray:
            x = x.tolist()
        featstest = StringCharFeatures(DNA)
        featstest.set_features(x)
        output = self.svm.apply(featstest)
        return output.get_values(), output.get_labels()

    def cross_validation(self,x,y,C=1,k=5,train_size=0.7,test_size=0.3,stratify=None):
        acc = []

        for i in range(k):
            self.x, X_test, self.y, y_test = train_test_split(x,y,
                                                              test_size=test_size,
                                                              train_size=train_size,
                                                              stratify=stratify)
            self.train(C)
            acc.append(self.accuracy(X_test,y_test))
        return acc

    def reduce_x(self, idx, x=None, out=False):
        if out == False:
            if x == None:
                if self.x != None and type(self.x) == list:
                    self.x = np.array(self.x)
                self.x = self.x[idx]
                self.x = self.x.tolist()
            else:
                self.x=x
                self.reduce_x(idx)
        else:
            if x != None and type(x) == list:
                x = np.array(x)
            return x[idx]

    def reduce_y(self, idx, y=None, out=False):
        if out == False:
            if y == None:
                if self.y != None and type(self.y) == list:
                    self.y = np.array(self.y)
                self.y = self.y[idx]
            else:
                self.y = y
                self.reduce_y(idx)
        else:
            if y != None and type(y) == list:
                y = np.array(y)
            return y[idx]


    def stratified_coss_validation(self,x,y,nsplits=5,C=1,shuffle=False,random_state=None):
        acc = []
        skf = StratifiedKFold(n_splits=nsplits,random_state=random_state,shuffle=shuffle)
        for train_idx, test_idx in skf.split(x,y):
            self.reduce_x(train_idx,x)
            self.reduce_y(train_idx,y)
            self.train(C)
            acc.append(self.accuracy(self.reduce_x(test_idx,x,True),self.reduce_y(test_idx,y,True)))
        return acc


    def leave_one_out(self,x,y,C=1):
        loo = LeaveOneOut()
        acc = []
        for train_idx, test_idx in loo.split(x):
            self.reduce_x(train_idx, x)
            self.reduce_y(train_idx, y)
            self.train(C)
            acc.append(self.accuracy(self.reduce_x(test_idx, x, True), self.reduce_y(test_idx, y, True)))
            print acc
        return acc


    def tune_C(self,x,y,nsplits,values):
        val_opt = values[0]
        acc_opt = 0
        for val in values:
            acc = np.mean(self.stratified_coss_validation(x,y,nsplits,C=val))
            if acc > acc_opt:
                acc_opt=acc
                val_opt=val
        return val_opt, acc_opt


    def accuracy(self,xtest,ytest):
        if ytest != None and type(ytest) != np.ndarray:
            ytest = np.array(ytest)
        eval = AccuracyMeasure()
        _, y_pred = self.apply(xtest)
        return eval.evaluate(BinaryLabels(y_pred), BinaryLabels(ytest))


    def category(self,x,y):

        acc, y_pred = self.apply(x)
        if y != None and type(y) != np.ndarray:
            y = np.array(y)
        category = {}
        category['tp'] = np.where((y == 1) * (y_pred == 1))[0]    # TP
        category['fp'] = np.where((y == -1) * (y_pred == 1))[0]   # FP
        category['fn'] = np.where((y == 1) * (y_pred == -1))[0]   # FN
        category['tn'] = np.where((y == -1) * (y_pred == -1))[0]  # TN
        return category