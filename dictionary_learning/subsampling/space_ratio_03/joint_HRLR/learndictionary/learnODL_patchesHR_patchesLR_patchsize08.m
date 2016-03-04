clear all; close all; clc;

%% INITIAL PARAMS
space_spacing=3;
time_spacing=4;

patchsize_l = 8; % 8x8 LR patches
patchsize_h = patchsize_l*space_spacing; % size of HR patches

num_patch=8*8; % number of patches extracted from each plane
dim_h=patchsize_h^2;
dim_l=patchsize_l^2;

filename_ref='/data/ISOTROPIC/data/data_downsampled4.nc';
nc = netcdf(filename_ref,'r');
Nt = nc('Nt').itsDimsize;
Nh = nc('Nx').itsDimsize;
Nl=Nh/space_spacing;

u_HR_all=nc{'velocity_x'}(1:Nt,1:Nh,1:Nh,1:space_spacing:Nh);
close(nc)

LTHS_idt=1:time_spacing:Nh;
[gridz_HR, gridy_HR] = meshgrid(1:Nh,1:Nh);
gridz_LR = gridz_HR(1:space_spacing:Nh,1:space_spacing:Nh);
gridy_LR = gridy_HR(1:space_spacing:Nh,1:space_spacing:Nh);

%% Load patches
PATCHES_FILENAME=strcat('/data/ISOTROPIC/dictionary_learning/space_ratio_03/trainingpatches_coupleHRLR_spaceratio',num2str(space_spacing,'%.1d'),'_timeratio'...
    ,num2str(time_spacing,'%.1d'),'_patchsize',num2str(patchsize_l,'%.2d'),'_numpatch',num2str(Nt*numel(LTHS_idt)*num_patch,'%.6d'),'.mat');
load(PATCHES_FILENAME, 'patches_HR_all','patches_LR_all');

% SUBSTRACT MEAN AND PRE-NORMALIZE ALL PATCHES
patches_HR_all = patches_HR_all - repmat(mean(patches_HR_all,1),dim_h,1);
patches_HR_all = patches_HR_all./repmat(sqrt(sum(patches_HR_all.^2,1)), dim_h, 1);
patches_LR_all = patches_LR_all - repmat(mean(patches_LR_all,1),dim_l,1);
patches_LR_all = patches_LR_all./repmat(sqrt(sum(patches_LR_all.^2,1)), dim_l, 1);

% JOINT PATCHES
patches_all = [(1/sqrt(dim_h))*patches_HR_all; (1/sqrt(dim_l))*patches_LR_all];
patch_norm = sqrt(sum(patches_all.^2, 1));
patches_all = patches_all(:, patch_norm > 1e-5); 
patches_all = patches_all - repmat(mean(patches_all,1),dim_h+dim_l,1);
patches_all = patches_all./repmat(patch_norm, dim_h+dim_l, 1);

clearvars patches_HR_all patches_LR_all;


%% ONLINE DICTIONARY LEARNING
params.mode=2; % min_{D in C} (1/n) sum_{i=1}^n (1/2)||x_i-Dalpha_i||_2^2 + lambda||alpha_i||_1 + lambda_2||alpha_i||_2^2
params.K=2*(dim_h+dim_l); 
params.lambda=0.05;
params.lambda2=0;
params.numThreads=4; % number of threads
params.iter=1000;  % max number of iterations.

fprintf(['Train dictionary for lambda=',num2str(params.lambda,'%.3f'),', K=',num2str(params.K,'%.4d'),'\n']);

[D,~] = mexTrainDL(patches_all,params);
CoefMatrix=mexLasso(patches_all,D,params);

D_HR = D(1:dim_h,:);
D_LR = D(dim_h+1:end,:);

ODL_FILENAME=strcat('/data/ISOTROPIC/dictionary_learning/space_ratio_03/DICTIONARY_coupleHRLR_patchesHR_patchesLR_joint_K',num2str(params.K,'%.4d'),'_lambda',strrep(num2str(params.lambda,'%.2f'),'.',''),'.mat');
save(ODL_FILENAME,'D_HR','D_LR','CoefMatrix');
