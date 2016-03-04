clear all; close all; clc;

%% INITIAL PARAMS
space_spacing=3;
time_spacing=4;

dim_l_pca=81;
size_l=8; 
patchsize_h = size_l*space_spacing; % size of HR patches
patchsize_l = patchsize_h; % 4x4 LR patches

num_patch=8*8; % number of patches 
dim_h=patchsize_h^2;
dim_l=4*patchsize_l^2;

params_train.lambda=0.2; % params.lambda=1*1e-5, overlap=7
params_train.K=2*(dim_h+dim_l_pca);
filename_SR=strcat('/data/ISOTROPIC/dictionary_learning/downsampling/space_ratio_03/SR_couplefeatures_K',num2str(params_train.K,'%.4d'),'_lambda',strrep(num2str(params_train.lambda,'%.2f'),'.',''),'.nc');
filename_ref='/data/ISOTROPIC/data/data_downsampled4.nc';
nc = netcdf(filename_ref,'r');
Nt = nc('Nt').itsDimsize;
Nx = nc('Nx').itsDimsize;
Ny = nc('Ny').itsDimsize;
Nz = nc('Nz').itsDimsize;
close(nc)



Nh=Nx;
Nl=Nh/space_spacing;
midplane_ids=3:time_spacing:Nh;


NRMSE= @(org,rec) sqrt(sum((org(:)-rec(:)).^2))/sqrt(sum((org(:)).^2));

%% Load data
nc1=netcdf(filename_ref,'r');
x_ref_all = nc1{'velocity_x'}(:,:,:,midplane_ids);
close(nc1);

x_interp_all=zeros(size(x_ref_all));
for t=1:Nt
    for i=1:numel(midplane_ids)
        x_ref=squeeze(x_ref_all(t,:,:,i));
        x_interp_all(t,:,:,i)=resize_nguyen(resize_nguyen(x_ref, 1/space_spacing,'bicubic'), space_spacing,'bicubic');       
    end
end
nc2=netcdf(filename_SR,'r');
x_SR_all = nc2{'Zhat_all'}(:,:,:,:);
close(nc2);

NRMSE_interp = NRMSE(x_ref_all(:),x_interp_all(:))
NRMSE_SR = NRMSE(x_ref_all(:),x_SR_all(:))