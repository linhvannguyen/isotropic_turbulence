% clear all; close all; clc;

%% INITIAL PARAMS
space_spacing=3;
time_spacing=4;

% size_l=8; 
% patchsize_h = size_l*space_spacing; % size of HR patches
% dim_h=patchsize_h^2;
% num_patch=8*8; % number of patches 

size_l=4; 
patchsize_h = size_l*space_spacing; % size of HR patches
dim_h=patchsize_h^2;
num_patch=16*16; % number of patches 


% size_l=2; 
% patchsize_h = size_l*space_spacing; % size of HR patches
% dim_h=patchsize_h^2;
% num_patch=16*16; % number of patches 

filename_ref='/data/ISOTROPIC/data/data_downsampled4.nc';
nc = netcdf(filename_ref,'r');
Nt = nc('Nt').itsDimsize;
Nh = nc('Nx').itsDimsize;
Nl=Nh/space_spacing;
close(nc)


overlap=7; % at LR

LTHS_idt=1:time_spacing:Nh;
[gridz_HR, gridy_HR] = meshgrid(1:Nh,1:Nh);
gridz_LR = gridz_HR(1:space_spacing:Nh,1:space_spacing:Nh);
gridy_LR = gridy_HR(1:space_spacing:Nh,1:space_spacing:Nh);

%% LOAD DICTIONARY 
filename_ODL=strcat('/data/ISOTROPIC/fusion/dictionary_smallscales_sizeLR',num2str(size_l,'%.2d'),'_timespacing',num2str(time_spacing,'%.2d'),'_spacespacing',num2str(space_spacing,'%.2d'),'.mat');
load(filename_ODL);

D = D./repmat(sqrt(sum(D.^2, 1)), dim_h, 1);

%% PLOT DICTIONARY
[~,ids]=sort(sum(CoefMatrix.^2,2),'descend');
fig_HR=figure(); set(fig_HR, 'Color', 'w'); displayPatches(D(:,ids(1:1:end))); axis off;
