# -*- coding: utf-8 -*-
"""
Created on Thu March 1 23:17:41 2016

@author: nguyen
"""
import numpy as np
import math
from netCDF4 import Dataset

# Constants
Nh = 96
Nt = 37
sspacing = 3
tspacing = 4

HTLS_sknots = np.arange(0,Nh,sspacing)
HTHS_sknots = np.arange(0,Nh,1)
LTHS_tknots = np.arange(0,Nh,tspacing)
Nl = len(HTLS_sknots)
Ns = len(LTHS_tknots)

Dh = Nh*Nh
Dl = Nl*Nl

N = Nt*Ns

#Load all training data
Xh_tr = np.zeros((N, Dh))
Xl_tr = np.zeros((N, Dl))
ncfile1 = Dataset('/data/ISOTROPIC/data/data_downsampled4.nc','r')
for t in range(Nt):
    count = 0
    for i in LTHS_tknots:
        xh = np.array(ncfile1.variables['velocity_x'][t,0:Nh,0:Nh,i])
        xl = xh[0:-1:sspacing,0:-1:sspacing] # xh[np.meshgrid(HTLS_sknots,HTLS_sknots)]
        Xh_tr[t*Ns + count,:] = np.reshape(xh,(1, Dh))
        Xl_tr[t*Ns + count,:] = np.reshape(xl,(1, Dl))
        count = count + 1
ncfile1.close()


# normalized: centered, variance 1
mea_l = np.zeros(Dl)
sig_l = np.zeros(Dl)
for k in range(Dl):
    mea_l[k] = np.mean(Xl_tr[:,k])
    sig_l[k] = np.std(Xl_tr[:,k])
    Xl_tr[:,k] = (Xl_tr[:,k]-mea_l[k])/sig_l[k] 
    
mea_h = np.zeros(Dh)
sig_h = np.zeros(Dh)
for k in range(Dh):
    mea_h[k] = np.mean(Xh_tr[:,k])
    sig_h[k] = np.std(Xh_tr[:,k])
    Xh_tr[:,k] = (Xh_tr[:,k]-mea_h[k])/sig_h[k]     

# Learning curve 
from sklearn.kernel_ridge import KernelRidge
from sklearn.learning_curve import learning_curve
from sklearn import cross_validation

estimator = KernelRidge(kernel='rbf', alpha=math.pow(10, -3.5), gamma=10e-4)
cv = cross_validation.ShuffleSplit(np.shape(Xl_tr)[0], n_iter=50, test_size=0.1, random_state=0)
train_sizes, train_MSE, test_MSE = learning_curve(estimator, Xl_tr, Xh_tr, 
                                                        cv=cv, n_jobs=4,
                                                        train_sizes=np.linspace(.1, 1.0, 20),
                                                        scoring = "mean_squared_error")
                                                        
# save to .mat file
import scipy.io as io
io.savemat('/data/ISOTROPIC/regression/KRR_learningcurve.mat', 
           dict(train_sizes=train_sizes, train_MSE=train_MSE, test_MSE=test_MSE))
