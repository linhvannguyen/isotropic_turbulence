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


## RR
from sklearn import linear_model
n_alphas = 200
lambdas_RR = np.append(0,np.logspace(-3, 6, n_alphas))

clf = linear_model.Ridge( fit_intercept=False, normalize=False)
coefs_RR_lambdas = []
for a in lambdas_RR:
    clf.set_params(alpha=a)
    clf.fit(Xl_tr, Xh_tr[:,0:8])
    coefs_RR_lambdas.append(clf.coef_[0:8,0]) 
print('\n Finish computing coefficients as funcs of lambda')

## Lasso
n_alphas = 200
lambdas_Lasso = np.append(0,np.logspace(-3, 6, n_alphas))

clf = linear_model.Lasso( fit_intercept=False, normalize=False)
coefs_Lasso_lambdas = []
for a in lambdas_Lasso:
    clf.set_params(alpha=a)
    clf.fit(Xl_tr, Xh_tr[:,0:8])
    coefs_Lasso_lambdas.append(clf.coef_[0:8,0]) 
print('\n Finish computing coefficients as funcs of lambda')



# save to .mat file
import scipy.io as io
io.savemat('/data/ISOTROPIC/regression/RR_Lasso_coeffs_lambda.mat', 
           dict(coefs_RR_lambdas=coefs_RR_lambdas, lambdas_RR=lambdas_RR,
                coefs_Lasso_lambdas=coefs_Lasso_lambdas, lambdas_Lasso=lambdas_Lasso))