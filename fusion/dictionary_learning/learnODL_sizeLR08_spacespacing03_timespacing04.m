clear all; close all; clc;

%% INITIAL PARAMS
space_spacing=3;
time_spacing=4;

size_l=8; 
patchsize_h = size_l*space_spacing; % size of HR patches
dim_h=patchsize_h^2;
num_patch=8*8; % number of patches 
params.K=2*dim_h;

% size_l=4; 
% patchsize_h = size_l*space_spacing; % size of HR patches
% dim_h=patchsize_h^2;
% num_patch=16*16; % number of patches 
% params.K=2*dim_h;


% size_l=2; 
% patchsize_h = size_l*space_spacing; % size of HR patches
% dim_h=patchsize_h^2;
% num_patch=16*16; % number of patches 
% params.K=4*dim_h;

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
filename_patches=strcat('/data/ISOTROPIC/fusion/patches_smallscales_sizeLR',num2str(size_l,'%.2d'),'_timespacing',num2str(time_spacing,'%.2d'),'_spacespacing',num2str(space_spacing,'%.2d'),'.mat');
load(filename_patches, 'patches_smallscales_all');

% SUBSTRACT MEAN AND PRE-NORMALIZE ALL PATCHES
patches_smallscales_all=patches_smallscales_all-repmat(mean(patches_smallscales_all),[dim_h 1]);
patches_smallscales_all=patches_smallscales_all ./ repmat(sqrt(sum(patches_smallscales_all.^2)),[dim_h 1]); 

%% ONLINE DICTIONARY LEARNING
params.mode=2; % min_{D in C} (1/n) sum_{i=1}^n (1/2)||x_i-Dalpha_i||_2^2 + lambda||alpha_i||_1 + lambda_2||alpha_i||_2^2
params.lambda=0.05;
params.lambda2=0;
params.numThreads=2; % number of threads
params.iter=2000;  % max number of iterations.

% LEARN DICT
fprintf(['Train dictionary for lambda=',num2str(params.lambda,'%.3f'),', K=',num2str(params.K,'%.4d'),'\n']);

[D,~] = mexTrainDL(patches_smallscales_all,params);
CoefMatrix=mexLasso(patches_smallscales_all,D,params);

filename_ODL=strcat('/data/ISOTROPIC/fusion/dictionary_smallscales_sizeLR',num2str(size_l,'%.2d'),'_timespacing',num2str(time_spacing,'%.2d'),'_spacespacing',num2str(space_spacing,'%.2d'),'.mat');
save(filename_ODL,'D','CoefMatrix');