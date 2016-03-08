# -*- coding: utf-8 -*-
"""
Created on Thu March 1 23:17:41 2016

@author: nguyen
"""
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
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


############## RidgeCV ########################################################
# Load optimal lambda
import scipy.io as sio
mf = sio.loadmat('/data/ISOTROPIC/regression/RR_crossvalidation.mat', squeeze_me=True, struct_as_record=False)
RR_lambda_opt = mf['RR_lambda_opt']
print'Optimal lambda:', RR_lambda_opt

from sklearn.linear_model import Ridge
ridge = Ridge(alpha=RR_lambda_opt, fit_intercept=False, normalize=False)
ridge.fit(Xl_tr, Xh_tr)
print np.shape(ridge.coef_)

coefs_RR=np.reshape(ridge.coef_[:,Dl/2+Nl/2],(96,96))

ridge = Ridge(alpha=0.0, fit_intercept=False, normalize=False)
ridge.fit(Xl_tr, Xh_tr)

coefs_LSE=np.reshape(ridge.coef_[:,Dl/2+Nl/2],(96,96))

##
fig1 = plt.figure(figsize=(8, 7))
plt.imshow(coefs_RR,interpolation='none')    
plt.show()

fig2 = plt.figure(figsize=(8, 7))
plt.imshow(coefs_LSE,interpolation='none')    
plt.show()

# save to .mat file
sio.savemat('/data/ISOTROPIC/regression/RR_sample_coefficients.mat', dict(coefs_RR=coefs_RR, coefs_LSE=coefs_LSE))