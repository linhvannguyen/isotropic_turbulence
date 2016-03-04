clear all; close all; clc;

%% INITIAL PARAMS
space_spacing=3;
time_spacing=4;

size_l=8; 
patchsize_h = size_l*space_spacing; % size of HR patches
patchsize_l = patchsize_h; % 4x4 LR patches

num_patch=8*8; % number of patches 
dim_h=patchsize_h^2;
dim_l=4*patchsize_l^2;

filename_ref='/data/ISOTROPIC/data/data_downsampled4.nc';
nc = netcdf(filename_ref,'r');
Nt = nc('Nt').itsDimsize;
Nh = nc('Nx').itsDimsize;
Nl=Nh/space_spacing;

close(nc)

LTHS_idt=1:time_spacing:Nh;
[gridz_HR, gridy_HR] = meshgrid(1:Nh,1:Nh);
gridz_LR = gridz_HR(1:space_spacing:Nh,1:space_spacing:Nh);
gridy_LR = gridy_HR(1:space_spacing:Nh,1:space_spacing:Nh);

%% Load patches
PATCHES_FILENAME=strcat('/data/ISOTROPIC/dictionary_learning/downsampling/space_ratio_03/trainingpatches_couplefeatures_spaceratio',num2str(space_spacing,'%.1d'),'_timeratio',num2str(time_spacing,'%.1d'),'_patchsize',num2str(size_l,'%.2d'),'_numpatch',num2str(Nt*numel(LTHS_idt)*num_patch,'%.6d'),'.mat');

load(PATCHES_FILENAME, 'patches_HRLR_all','patches_LR_features_all');

% SUBSTRACT MEAN AND PRE-NORMALIZE ALL PATCHES
start=tic();
patches_LR_features_all=patches_LR_features_all-repmat(mean(patches_LR_features_all),[dim_l 1]);
patches_LR_features_all=patches_LR_features_all ./ repmat(sqrt(sum(patches_LR_features_all.^2)),[dim_l 1]); 

patches_HRLR_all=patches_HRLR_all-repmat(mean(patches_HRLR_all),[dim_h 1]);
patches_HRLR_all=patches_HRLR_all ./ repmat(sqrt(sum(patches_HRLR_all.^2)),[dim_h 1]); 

fprintf(['Processing data in ',num2str(toc(start),'%.3f'),' seconds \n']);

%% Dimension reduction
start=tic();
C = patches_LR_features_all * patches_LR_features_all';
[V, D] = eig(C);
D = diag(D);
D = cumsum(D) / sum(D);
% k = find(D >= 1e-3, 1); % ignore 0.1% energy

dim_l_pca=81;
V_pca = V(:, end-dim_l_pca+1:end); % ignore 0.23% energy
patches_LR_features_pca = V_pca' * patches_LR_features_all;

patches_LR_features_pca=patches_LR_features_pca-repmat(mean(patches_LR_features_pca),[dim_l_pca 1]);
patches_LR_features_pca=patches_LR_features_pca ./ repmat(sqrt(sum(patches_LR_features_pca.^2)),[dim_l_pca 1]); 

fprintf(['PCA in ',num2str(toc(start),'%.3f'),' seconds \n']);

clearvars patches_LR_features_all;

%% ONLINE DICTIONARY LEARNING
params.mode=2; % min_{D in C} (1/n) sum_{i=1}^n (1/2)||x_i-Dalpha_i||_2^2 + lambda||alpha_i||_1 + lambda_2||alpha_i||_2^2
params.lambda=0.2;
params.lambda2=0;
params.K=2*(dim_h+dim_l_pca);
params.numThreads=4; % number of threads
params.iter=1000;  % max number of iterations.

% LEARN DICT
fprintf(['Train dictionary for lambda=',num2str(params.lambda,'%.3f'),', K=',num2str(params.K,'%.4d'),'\n']);

[D_LR,~] = mexTrainDL(patches_LR_features_pca,params);
CoefMatrix=mexLasso(patches_LR_features_pca,D_LR,params);

D_HR = (patches_HRLR_all * CoefMatrix')/(full(CoefMatrix * CoefMatrix'));

ODL_FILENAME=strcat('/data/ISOTROPIC/dictionary_learning/downsampling/space_ratio_03/DICTIONARY_couplefeatures_patchesHR_patchesLR_joint_K',num2str(params.K,'%.4d'),'_lambda',strrep(num2str(params.lambda,'%.2f'),'.',''),'.mat');
save(ODL_FILENAME,'D_HR','D_LR','V_pca','CoefMatrix');