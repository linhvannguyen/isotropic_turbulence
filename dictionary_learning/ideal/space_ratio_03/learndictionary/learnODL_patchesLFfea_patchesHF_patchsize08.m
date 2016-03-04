clear all; close all; clc;

%% INITIAL PARAMS
space_spacing=3;
time_spacing=4;

patchsize_l = 8; % 8x8 LR patches
patchsize_h = patchsize_l*space_spacing; % size of HR patches
dim_patch=patchsize_h^2;

num_patch=16*16; % number of patches extracted from each plane

filename_ref='/data/ISOTROPIC/data/data_downsampled4.nc';
nc = netcdf(filename_ref,'r');
Nt = nc('Nt').itsDimsize;
Nh = nc('Nx').itsDimsize;
Nl=Nh/space_spacing;

u_HR_all=nc{'velocity_x'}(1:Nt,1:Nh,1:Nh,1:space_spacing:Nh);
close(nc)

LTHS_idt=1:time_spacing:Nh;

%% Load patches
PATCHES_FILENAME=strcat('/data/ISOTROPIC/dictionary_learning/ideal/space_ratio_03/trainingpatches_sspacing',num2str(space_spacing,'%.1d'),...
    '_tspacing',num2str(time_spacing,'%.1d'),'_patchsize',num2str(patchsize_l,'%.2d'),'_numpatch',num2str(Nt*numel(LTHS_idt)*num_patch,'%.6d'),'.mat');
load(PATCHES_FILENAME, 'patches_lf_features_pca','patches_hf_all', 'V_pca');

% SUBSTRACT MEAN AND PRE-NORMALIZE ALL PATCHES
dim_pca = size(V_pca,2);
patches_lf_features_pca = patches_lf_features_pca - repmat(mean(patches_lf_features_pca,1),dim_pca,1); 
patches_lf_features_pca = patches_lf_features_pca./repmat(sqrt(sum(patches_lf_features_pca.^2,1)), dim_pca, 1);
patches_hf_all = patches_hf_all - repmat(mean(patches_hf_all,1),dim_patch,1);
patches_hf_all = patches_hf_all./repmat(sqrt(sum(patches_hf_all.^2,1)), dim_patch, 1);

%% ONLINE DICTIONARY LEARNING
params.mode=2; % min_{D in C} (1/n) sum_{i=1}^n (1/2)||x_i-Dalpha_i||_2^2 + lambda||alpha_i||_1 + lambda_2||alpha_i||_2^2
params.K=2*dim_patch; 
params.lambda=0.2;
params.lambda2=0;
params.numThreads=3; % number of threads
params.iter=2000;  % max number of iterations.

% LEARN DICT
fprintf(['Train dictionary for lambda=',num2str(params.lambda,'%.3f'),', K=',num2str(params.K,'%.4d'),'\n']);

[D_LF_FEA,~] = mexTrainDL(patches_lf_features_pca,params);
CoefMatrix=mexLasso(patches_lf_features_pca,D_LF_FEA,params);

D_HF = (patches_hf_all * CoefMatrix')/(full(CoefMatrix * CoefMatrix'));

ODL_FILENAME=strcat('/data/ISOTROPIC/dictionary_learning/ideal/space_ratio_03/DICTIONARY_coupleLFfeaHF_K',num2str(params.K,'%.4d'),'_lambda',strrep(num2str(params.lambda,'%.2f'),'.',''),'.mat');
save(ODL_FILENAME,'D_LF_FEA','D_HF','CoefMatrix');
