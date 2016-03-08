# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 14:20:44 2016

@author: nguyen
"""
import numpy as np
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
        xl = xh[0:-1:sspacing,0:-1:sspacing]
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


from sklearn.kernel_ridge import KernelRidge

"""
# Use GridSearchCV to get a rough idea of parameters
from sklearn.grid_search import GridSearchCV

kr = GridSearchCV(KernelRidge(kernel='poly', degree=3), cv=5,
                  param_grid={"alpha": np.append(0,np.logspace(-10, -1, 11))})
kr.fit(Xl_tr, Xh_tr)

print('\n Best estimator:', kr.best_estimator_)
"""

"""                                                        
# Validation curve for lambda
from sklearn.learning_curve import validation_curve
lambdas_range = np.logspace(-10, 0, 44)
train_MSE_lambda, test_MSE_lambda = validation_curve(KernelRidge(kernel='poly', degree=2), 
                                                     Xl_tr, Xh_tr, param_name="alpha", param_range=lambdas_range, 
                                                     scoring = "mean_squared_error", cv=10)                                                      
print('\n Finish validation curve for lambda:')
"""
# Validation curve for poly degree
from sklearn.learning_curve import validation_curve
degree_range = np.linspace(1, 10, 10)
train_MSE_degree, test_MSE_degree = validation_curve(KernelRidge(kernel='poly'), 
                                                     Xl_tr, Xh_tr, param_name="degree", param_range=degree_range, 
                                                     scoring = "mean_squared_error", cv=10)                                                      
print('\n Finish validation curve for poly degree')

# save to .mat file
import scipy.io as io
io.savemat('/data/ISOTROPIC/regression/KRR_poly_validationcurve.mat', 
           dict(degree_range=degree_range, train_MSE_degree = train_MSE_degree, test_MSE_degree = test_MSE_degree))             