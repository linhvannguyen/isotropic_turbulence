clear all; close all; clc;

%% INITIAL PARAMS
space_spacing=3;
time_spacing=4;

dim_l_pca=81;
size_l=8; 
patchsize_h = size_l*space_spacing; % size of HR patches
patchsize_l = patchsize_h;

num_patch=8*8; % number of patches 
dim_h=patchsize_h^2;
dim_l=4*patchsize_l^2;

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
ODL_FILENAME=strcat('/data/ISOTROPIC/dictionary_learning/subsampling/space_ratio_03/DICTIONARY_alltranslations_K1400_lambda005.mat');
load(ODL_FILENAME);
dim_l_pca = size(V_pca,2);

params.lambda=0.2; % params.lambda=1*1e-5, overlap=7
params.K=2*(dim_h+dim_l_pca);

D_HR = D(1:dim_h,:);
D_LR = D(1+dim_h:end,:);

D_HR = D_HR./repmat(sqrt(sum(D_HR.^2, 1)), dim_h, 1);
D_HR=sqrt(dim_h/dim_l_pca)*D_HR; % scale the HR dictionary (very important!!!)
D_LR = D_LR./repmat(sqrt(sum(D_LR.^2, 1)), dim_l_pca, 1);

%% PLOT DICTIONARY
% [~,ids]=sort(sum(CoefMatrix.^2,2),'descend');
fig_HR=figure(); set(fig_HR, 'Color', 'w'); displayPatches(D_HR(:,1:5:end)); axis off;

%% FILTERS
% O = zeros(1, space_spacing-1);
O=0;
hf1 = [-1,O,1]; vf1 = hf1';
hf2 = [1,O,-2,O,1]; vf2 = hf2';

%% GRID OF FIRST NODES IN LR AND HR PATCHES (enlarged version)
gridz_l = 1: size_l-overlap : Nl+1;
gridy_l = 1: size_l-overlap : Nl+1;
gridz_h = 1: space_spacing*(size_l-overlap) : Nh+1;
gridy_h = 1: space_spacing*(size_l-overlap) : Nh+1;  
grids={gridz_h,gridy_h};
sizes={Nl+size_l,Nl+size_l,Nh+patchsize_h,Nh+patchsize_h};

%% SR
params.mode=2; % min_{D in C} (1/n) sum_{i=1}^n (1/2)||x_i-Dalpha_i||_2^2 + lambda||alpha_i||_1 + lambda_2||alpha_i||_2^2
params.lambda=1*1e-5; 
params.lambda2=0;
scaling_factor=1; % is this sqrt(patchsize_h*patchsize_h/dim_l_pca) (currently dim_l_pca=100; best: 2.15)
fprintf('Increasing resolution...\n');

% Load data
t=10; i=3;
nc = netcdf(filename_ref,'r');
X_HR_org=nc{'velocity_x'}(t,1:Nh,1:Nh,i);
close(nc)

left=patchsize_h/2; right = patchsize_h/2;
bottom=patchsize_h/2; top = patchsize_h/2;

[gridz_HS_enlarged, gridy_HS_enlarged] = meshgrid(1:Nh+left+right,1:Nh+bottom+top);
gridz_LS_enlarged=gridz_HS_enlarged(1:space_spacing:end,1:space_spacing:end);
gridy_LS_enlarged=gridy_HS_enlarged(1:space_spacing:end,1:space_spacing:end);

X_HR_org_enlarged=enlarge_2D(X_HR_org,left, right, bottom, top);


% LR field
X_LR_enlarged=resize_nguyen(X_HR_org_enlarged, 1/space_spacing,'bicubic');

X_LR_interp_enlarged=resize_nguyen(X_LR_enlarged, space_spacing,'bicubic');
% X_LR_interp_enlarged=interp2(gridz_LS_enlarged,gridy_LS_enlarged, X_LR_enlarged, gridz_HS_enlarged, gridy_HS_enlarged,'spline');

% SR
start=tic();
 
X_HR_rec_smallscales_enlarged=SR_nguyen(X_LR_interp_enlarged, D_HR, D_LR, params,'FLUC_FEATURES_INTERP',grids, sizes, {hf1, vf1, hf2, vf2}, V_pca,scaling_factor);
X_HR_rec_enlarged = X_LR_interp_enlarged+X_HR_rec_smallscales_enlarged;

% REMOVE BOLDER (translation is not complete near bolder)
X_HR_rec=X_HR_rec_enlarged(bottom+1:bottom+Nh,left+1:left+Nh);
X_LR_interp=X_LR_interp_enlarged(bottom+1:bottom+Nh,left+1:left+Nh);

NRMSE1=sqrt(sum((X_LR_interp(:)-X_HR_org(:)).^2))/sqrt(sum(X_HR_org(:).^2));
NRMSE2=sqrt(sum((X_HR_rec(:)-X_HR_org(:)).^2))/sqrt(sum(X_HR_org(:).^2));
fprintf(['SR snapshot t=',num2str(t,'%.3d'),' in ',num2str(toc(start)/60,'%.3f'),' minute(s), improves NRMSE by ', num2str((1-NRMSE2/NRMSE1)*100,'%.2f'),' percents\n']);


%% PLOT
fig=figure();
set(fig, 'Position', [200 200 700 650]);
set(fig,'color','w')

ax(1)=subplot(2,2,1,'position',[0.05 0.525 0.425 0.425]); 
imagesc(X_HR_org); caxis([-3 3]); 
set(gca,'YDir','normal'); axis off; axis equal;
title('Reference')

ax(2)=axes('position',[0.05 0.05 0.425 0.425]); 
imagesc(X_LR_interp); caxis([-3 3]); set(gca,'YDir','normal'); axis off; axis equal;
title('Interpolation')

ax(3)=axes('position',[0.5 0.525 0.425 0.425]); 
imagesc(X_HR_rec); caxis([-3 3]); set(gca,'YDir','normal'); axis off; axis equal;
title('SR')

ax(4)=axes('position',[0.5 0.05 0.425 0.425]); 
imagesc(X_HR_rec-X_LR_interp); caxis([-0.2 0.2]); set(gca,'YDir','normal'); axis off; axis equal;
title('SR')

%%
E = estimate_spect_2D(X_HR_rec-X_LR_interp); 
