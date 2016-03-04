# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 23:17:41 2016

@author: nguyen
"""
import numpy as np
from scipy import interpolate
from sklearn.learning_curve import learning_curve

def plot_learning_curve(estimator, plt, X, y, ylim=None, cv=None, n_jobs=1, train_sizes=np.linspace(.1, 1.0, 5)):
    """
    Generate a simple plot of the test and traning learning curve.

    Parameters
    ----------
    estimator : object type that implements the "fit" and "predict" methods
        An object of that type which is cloned for each validation.

    plt : current matplotlib plot

    X : array-like, shape (n_samples, n_features)
        Training vector, where n_samples is the number of samples and
        n_features is the number of features.

    y : array-like, shape (n_samples) or (n_samples, n_features), optional
        Target relative to X for classification or regression;
        None for unsupervised learning.

    ylim : tuple, shape (ymin, ymax), optional
        Defines minimum and maximum yvalues plotted.

    cv : integer, cross-validation generator, optional
        If an integer is passed, it is the number of folds (defaults to 3).
        Specific cross-validation objects can be passed, see
        sklearn.cross_validation module for the list of possible objects

    n_jobs : integer, optional
        Number of jobs to run in parallel (default 1).
    """
    if ylim is not None:
        plt.ylim(*ylim)
    plt.xlabel("Number of training examples")
    plt.ylabel("Score")
    train_sizes, train_scores, test_scores = learning_curve(estimator, X, y, cv=cv, 
                                                            n_jobs=n_jobs, train_sizes=train_sizes)
    
    train_scores_mean = np.mean(train_scores, axis=1)
    train_scores_std = np.std(train_scores, axis=1)
    test_scores_mean = np.mean(test_scores, axis=1)
    test_scores_std = np.std(test_scores, axis=1)

    plt.fill_between(train_sizes, train_scores_mean - train_scores_std,
                     train_scores_mean + train_scores_std, alpha=0.1, color="r")
    plt.fill_between(train_sizes, test_scores_mean - test_scores_std,
                     test_scores_mean + test_scores_std, alpha=0.1, color="g")
    plt.plot(train_sizes, train_scores_mean, 'o-', color="r", label="Training score")
    plt.plot(train_sizes, test_scores_mean, 'o-', color="g", label="Cross-validation score")
    
    plt.grid()
    plt.legend(loc="best")
    return plt
    
    

def interp2 (x, y, z, xnew, ynew, kind='cubic'):
    f = interpolate.interp2d(x, y, z, kind=kind)
    return f(xnew, ynew)

def NRMSE (xref, xrec):
    err = np.sqrt(np.sum(np.square(xref.ravel()-xrec.ravel())))/np.sqrt(np.sum(np.square(xref.ravel())))
    return err
    
def gen_grids (sspacing, ispacing, Nh, Np):
    # ispacing = sspacing in case extract LR as predictors
    # ispacing = 1 in case extract HR to predict HR (say LF to predict HF)
    HTLS_sknots = np.arange(0, Nh,sspacing)
    Nl = len(HTLS_sknots)
    zz, yy = np.meshgrid(HTLS_sknots, HTLS_sknots)
    
    if ispacing == sspacing:
        offset = np.arange(-Np*sspacing,(Np+1)*sspacing+1,sspacing)
    elif ispacing == 1:
        offset = np.arange(-Np*sspacing,(Np+1)*sspacing+1,1)
    patchsize_o = len(offset) 
    # center HR patches: grid of size [Nl, Nl, sspacing+1, sspacing+ 1] is position of Nl x Nl HR patches
    # each is of size (sspacing+1) x (sspacing+1), each take (i,j), i,j = 0, ..., Nl-1 as the bottom-left point
    zz_o = np.reshape(zz,(Nl, Nl, 1, 1))
    zz_o = np.tile(zz_o, [1, 1, sspacing+1, sspacing+1])
    yy_o = np.reshape(yy,(Nl, Nl, 1, 1))
    yy_o = np.tile(yy_o, [1, 1, sspacing+1, sspacing+1])
    
    off_z, off_y = np.meshgrid(np.arange(sspacing+1),np.arange(sspacing+1))
    off_z = np.reshape(off_z,(1,1,sspacing+1,sspacing+1))
    off_z = np.tile(off_z, [Nl, Nl, 1, 1])
    off_y = np.reshape(off_y,(1,1,sspacing+1,sspacing+1))
    off_y = np.tile(off_y, [Nl, Nl, 1, 1])
    
    zz_o = zz_o + off_z
    zz_o[zz_o < 0] = Nh + zz_o[zz_o < 0]
    zz_o[zz_o > Nh-1] = zz_o[zz_o > Nh-1] - Nh    
    
    yy_o = yy_o + off_y    
    yy_o[yy_o < 0] = Nh + yy_o[yy_o < 0]
    yy_o[yy_o > Nh-1] = yy_o[yy_o > Nh-1] - Nh   

    # LR neighbors: grid of size [Nl, Nl, 4, 4] is position of Nl x Nl neihbor LR
    # of each HR patches
    zz_i = np.reshape(zz,(Nl, Nl, 1, 1))
    zz_i = np.tile(zz_i, [1, 1, patchsize_o, patchsize_o])
    yy_i = np.reshape(yy,(Nl, Nl, 1, 1))
    yy_i = np.tile(yy_i, [1, 1, patchsize_o, patchsize_o])
    
    off_z, off_y = np.meshgrid(offset, offset)
    off_z = np.reshape(off_z,(1,1,patchsize_o,patchsize_o))
    off_z = np.tile(off_z, [Nl, Nl, 1, 1])
    off_y = np.reshape(off_y,(1,1,patchsize_o,patchsize_o))
    off_y = np.tile(off_y, [Nl, Nl, 1, 1])
    
    zz_i = zz_i + off_z
    zz_i[zz_i < 0] = Nh + zz_i[zz_i < 0]
    zz_i[zz_i > Nh-1] = zz_i[zz_i > Nh-1] - Nh    
    
    yy_i = yy_i + off_y    
    yy_i[yy_i < 0] = Nh + yy_i[yy_i < 0]
    yy_i[yy_i > Nh-1] = yy_i[yy_i > Nh-1] - Nh       
    return (zz_o, yy_o, zz_i, yy_i)