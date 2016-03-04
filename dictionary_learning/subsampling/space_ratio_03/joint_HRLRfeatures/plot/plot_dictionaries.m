clear all; close all; clc;
path(path,'/home/nguyen/Documents/Codes/isotropic/DictionaryLearning/funcs');

%% INITIAL PARAMS
Nh=96; % Number of grid-point in HR field
scale=3;
Nl=Nh/scale;

size_l=8; 
patchsize_h = size_l*scale; % size of HR patches
patchsize_l = patchsize_h; % 4x4 LR patches

dim_h=patchsize_h^2;
dim_l=4*patchsize_l^2;
dim_l_pca=81;

nc = netcdf('/data/ISOTROPIC/NS3D_DNS_FORCE_384/TRAININGPLANES_downsampled4_HR.nc','r');
Nx = nc('Nt').itsDimsize;
close(nc)

overlap=7; % at LR

params.lambda=0.2; % params.lambda=1*1e-5, overlap=7
params.K=2*(dim_h+dim_l_pca);

%% LOAD DICTIONARY 
ODL_FILENAME=strcat('/data/ISOTROPIC/DICTIONARYLEARNING/final/downsample4_3/joint_HRLRfeatures/dictionaries/patchsize08/DICTIONARY_patchesHRLR_patchesLRfeatures_K',num2str(params.K,'%.4d'),'_lambda',strrep(num2str(params.lambda,'%.2f'),'.',''),'.mat');
load(ODL_FILENAME);

D_HR = D_HR./repmat(sqrt(sum(D_HR.^2, 1)), dim_h, 1);
D_HR=sqrt(dim_h/dim_l_pca)*D_HR; % scale the HR dictionary (very important!!!)
D_LR = D_LR./repmat(sqrt(sum(D_LR.^2, 1)), dim_l_pca, 1);

%% PLOT DICTIONARY
[~,ids]=sort(sum(CoefMatrix.^2,2),'descend');
fig_HR=figure(); set(fig_HR, 'Color', 'w'); displayPatches(D_HR(:,ids(1:10:end))); axis off;
DL_FILENAME=strcat('./figures/DICTIONARY_patchsize08_patchesHRLR_patchesLRfeatures_K',num2str(params.K,'%.4d'),'_lambda',strrep(num2str(params.lambda,'%.2f'),'.',''));
export_fig(DL_FILENAME,'-eps');
close();
