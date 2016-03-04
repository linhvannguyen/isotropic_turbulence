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

%% Load patches
PATCHES_FILENAME=strcat('/data/ISOTROPIC/dictionary_learning/space_ratio_03/trainingpatches_coupleHRLR_spaceratio',num2str(space_spacing,'%.1d'),'_timeratio'...
    ,num2str(time_spacing,'%.1d'),'_patchsize',num2str(patchsize_l,'%.2d'),'_numpatch',num2str(Nt*numel(LTHS_idt)*num_patch,'%.6d'),'.mat');
load(PATCHES_FILENAME, 'patches_HR_all','patches_LR_all');

% SUBSTRACT MEAN AND PRE-NORMALIZE ALL PATCHES
patches_HR_all = patches_HR_all - repmat(mean(patches_HR_all,2),1,size(patches_HR_all,2));
patches_HR_all = patches_HR_all./repmat(sqrt(sum(patches_HR_all.^2)), dim_h, 1);
patches_LR_all = patches_LR_all - repmat(mean(patches_LR_all,2),1,size(patches_LR_all,2));
patches_LR_all = patches_LR_all./repmat(sqrt(sum(patches_LR_all.^2)), dim_l, 1);

%%
ids=randperm(size(patches_HR_all,2));

for i=1:10
    patch_LR=patches_LR_all(:,ids(i));
    patch_LR=(patch_LR-min(patch_LR))/(max(patch_LR)-min(patch_LR));
    patch_LR=reshape(patch_LR,8,8);
    
    patch_HR=patches_HR_all(:,ids(i));
    patch_HR=(patch_HR-min(patch_HR))/(max(patch_HR)-min(patch_HR));
    patch_HR=reshape(patch_HR,24,24);

    fig=figure();
    set(fig, 'Position', [200 200 900 1000]);
    set(fig,'color','w')
    imagesc(patch_HR); caxis([0 1]); 
    set(gca,'YDir','normal'); axis off; axis equal;

    PATCHES_HR=strcat('./figures/patchHR_',num2str(i,'%.2d'));
    export_fig(PATCHES_HR,'-depsc');
    close();
    
    fig=figure();
    set(fig, 'Position', [200 200 900 1000]);
    set(fig,'color','w')
    imagesc(patch_LR); caxis([0 1]); 
    set(gca,'YDir','normal'); axis off; axis equal;

    PATCHES_LR=strcat('./figures/patchLR_',num2str(i,'%.2d'));
    export_fig(PATCHES_LR,'-depsc');
%     export_fig test.png
    close();    
end