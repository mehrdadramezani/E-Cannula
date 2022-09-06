import scipy.io as sio
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sn

from sklearn.feature_selection import RFECV
from sklearn.svm import SVC
from sklearn.dummy import DummyClassifier

from sklearn.metrics import accuracy_score, confusion_matrix, precision_score, recall_score
from sklearn.model_selection import cross_val_score, cross_val_predict

from matplotlib import pyplot as plt
from matplotlib import gridspec
import time
from scipy.stats import ranksums
import pandas as pd




# load the data
mat_contents    = sio.loadmat('SWR_dFoF.mat')
cell_assemblies = sio.loadmat('Cell_assemblies.mat')['clust_ID_final'][0]

X = mat_contents['cluster_input_feature'][:,:]
y = np.ravel(mat_contents['cluster_target_label'])
classnum = len(np.unique(y))

n_cell_assemblies = len(cell_assemblies)
n_cell = np.shape(X)[1]

cell_id = np.zeros((n_cell,2),dtype=int)
cell_id[:,0] = np.arange(1,n_cell+1)
for i in range(n_cell_assemblies):
    cell_id[cell_assemblies[i][:,0]-1,1] = i+1
    
assembly_size = np.zeros(n_cell_assemblies,)
for i in range(n_cell_assemblies):
    assembly_size[i] = len(cell_assemblies[i][:,0])


# resampling to balance the data
from sklearn.utils import resample
balanced_input = []
balanced_target = []
min_size = min([sum(y==i) for i in range(1,5)])
for i in range(1,classnum+1):
    if min_size < 20:
        balanced_input.append(resample(X[y==i], replace=True, n_samples=20, random_state=123))
        balanced_target.append(resample(y[y==i], replace=True, n_samples=20, random_state=456))
    else:
        balanced_input.append(resample(X[y==i], replace=False, n_samples=20, random_state=123))
        balanced_target.append(resample(y[y==i], replace=False, n_samples=20, random_state=456))
X_balanced = np.concatenate(balanced_input)
y_balanced = np.concatenate(balanced_target)
use_balanced = True

# construct classifier and use SVM RFE
estimator = SVC(C=10,kernel="linear",class_weight='balanced',random_state=1)
selector = RFECV(estimator, step=3, cv=5, n_jobs = 10, scoring='accuracy')
# start fitting the model
selector = selector.fit(X_balanced, y_balanced)

# compute the cross-validated prediction result
np.random.seed(0)
if use_balanced:
    X_select = selector.transform(X_balanced)
    Y_pred = cross_val_predict(estimator,X_select,y_balanced,cv=10,n_jobs=5)
    # classification measurement of the model
    acc = [float(len(np.intersect1d(np.where(Y_pred==i), np.where(y_balanced==i))))/len(np.where(y_balanced==i)[0]) for i in range(1,classnum+1)] #accuracy_score(y, Y_pred)
    prc = precision_score(y_balanced, Y_pred, average=None)
    rec = recall_score(y_balanced, Y_pred, average=None)
    #confmat = confusion_matrix(y_balanced,Y_pred)