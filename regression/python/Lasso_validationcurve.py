# -*- coding: utf-8 -*-
"""
Created on Thu March 1 23:17:41 2016

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


############## LassoCV ########################################################
from sklearn.linear_model import MultiTaskLassoCV
n_alphas = 5
alphas = np.logspace(-10, 0, n_alphas)

lasso = MultiTaskLassoCV(alphas = alphas, cv = 5, fit_intercept=False, normalize=False,n_jobs=3)
lasso.fit(Xl_tr, Xh_tr)

Lasso_lambda_opt = lasso.alpha_

print('\n Optimal lambda:', Lasso_lambda_opt)
############ Validation curve #################################################
"""
# validation curve
from sklearn.linear_model import Lasso
from sklearn.learning_curve import validation_curve

lambdas_range= np.append(0, np.logspace(0, 6, 28))
train_MSE, test_MSE = validation_curve(Lasso(),Xl_tr, Xh_tr, param_name="alpha", param_range=lambdas_range, 
                                             scoring = "mean_squared_error", cv=10)

# API always tries to maximize a loss function, so MSE is actually in the flipped sign
train_MSE_mean = np.mean(-train_MSE, axis=1)
train_MSE_std = np.std(-train_MSE, axis=1)
test_MSE_mean = np.mean(-test_MSE, axis=1)
test_MSE_std = np.std(-test_MSE, axis=1)  

# save to .mat file
import scipy.io as io
io.savemat('/data/ISOTROPIC/regression/Lasso_crossvalidation.mat', 
           dict(Lasso_lambda_opt=Lasso_lambda_opt, lambdas_range=lambdas_range, train_MSE=train_MSE, test_MSE = test_MSE))
"""
           