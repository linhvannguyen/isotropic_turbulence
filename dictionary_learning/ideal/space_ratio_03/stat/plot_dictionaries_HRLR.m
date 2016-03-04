clear all; close all; clc;
path(path,'/home/nguyen/Documents/Codes/isotropic/dictionary_learning/funcs/')

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

params.lambda=0.2; % params.lambda=0.001, overlap=9
params.K=2*(dim_h+dim_l); 


%% LOAD DICTIONARY 
ODL_FILENAME=strcat('/data/ISOTROPIC/dictionary_learning/space_ratio_03/DICTIONARY_patchesHR_patchesLR_joint_K',num2str(params.K,'%.4d'),'_lambda',strrep(num2str(params.lambda,'%.2f'),'.',''),'.mat');
load(ODL_FILENAME);

D_HR = D_HR./repmat(sqrt(sum(D_HR.^2, 1)), dim_h, 1);
D_HR=sqrt(dim_h)/sqrt(dim_l)*D_HR; % scale the HR dictionary
D_LR = D_LR./repmat(sqrt(sum(D_LR.^2, 1)), dim_l, 1);

%% PLOT DICTIONARY
[~,ids]=sort(sum(CoefMatrix.^2,2),'descend');
fig_HR=figure(); set(fig_HR, 'Color', 'w'); ImD=displayPatches(D_HR(:,ids(1:10:end))); axis off;
DL_FILENAME=strcat('../figures/DICTIONARY_patchsize08_HR_patchesHR_patchesLR_K',num2str(params.K,'%.4d'),'_lambda',strrep(num2str(params.lambda,'%.2f'),'.',''));
export_fig(DL_FILENAME,'-depsc');
close();

fig_LR=figure(); set(fig_LR, 'Color', 'w'); ImD2=displayPatches(D_LR(:,ids(1:10:end))); axis off;
DL_FILENAME=strcat('../figures/DICTIONARY_patchsize08_LR_patchesHR_patchesLR_K',num2str(params.K,'%.4d'),'_lambda',strrep(num2str(params.lambda,'%.2f'),'.',''));
export_fig(DL_FILENAME,'-depsc');
close();
