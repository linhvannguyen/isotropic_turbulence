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

%% AF-HF
% load /data/ISOTROPIC/dictionary_learning/ideal/space_ratio_03/DICTIONARY_coupleAFHF_K1152_lambda020.mat
% 
% D_AF = D_AF./repmat(sqrt(sum(D_AF.^2, 1)), dim_patch, 1);
% D_HF = D_HF./repmat(sqrt(sum(D_HF.^2, 1)), dim_patch, 1);
% 
% [~,ids]=sort(sum(CoefMatrix.^2,2),'descend');
% fig_AF=figure(); set(fig_AF, 'Color', 'w'); ImD=displayPatches(D_AF(:,ids(1:10:end))); axis off;
% fig_HF=figure(); set(fig_HF, 'Color', 'w'); ImD2=displayPatches(D_HF(:,ids(1:10:end))); axis off;

%% LF-HF
% load /data/ISOTROPIC/dictionary_learning/ideal/space_ratio_03/DICTIONARY_coupleLFHF_K1152_lambda020.mat
% 
% D_LF = D_LF./repmat(sqrt(sum(D_LF.^2, 1)), dim_patch, 1);
% D_HF = D_HF./repmat(sqrt(sum(D_HF.^2, 1)), dim_patch, 1);
% 
% [~,ids]=sort(sum(CoefMatrix.^2,2),'descend');
% fig_LF=figure(); set(fig_LF, 'Color', 'w'); ImD=displayPatches(D_LF(:,ids(1:10:end))); axis off;
% fig_HF=figure(); set(fig_HF, 'Color', 'w'); ImD2=displayPatches(D_HF(:,ids(1:10:end))); axis off;


%% AF-LF
load /data/ISOTROPIC/dictionary_learning/ideal/space_ratio_03/DICTIONARY_coupleAFLF_K1152_lambda005.mat

D_AF = D_AF./repmat(sqrt(sum(D_AF.^2, 1)), dim_patch, 1);
D_LF = D_LF./repmat(sqrt(sum(D_LF.^2, 1)), dim_patch, 1);

[~,ids]=sort(sum(CoefMatrix.^2,2),'descend');
fig_AF=figure(); set(fig_AF, 'Color', 'w'); ImD=displayPatches(D_AF(:,ids(1:10:end))); axis off;
fig_LF=figure(); set(fig_LF, 'Color', 'w'); ImD2=displayPatches(D_LF(:,ids(1:10:end))); axis off;
