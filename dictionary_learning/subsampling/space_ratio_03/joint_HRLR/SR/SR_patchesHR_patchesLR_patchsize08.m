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

%% LOAD DICTIONARY 
overlap=7; % at LR

params.lambda=0.05;
params.K=2*(dim_h+dim_l); 

ODL_FILENAME=strcat('/data/ISOTROPIC/dictionary_learning/subsampling/space_ratio_03/DICTIONARY_coupleHRLR_patchesHR_patchesLR_joint_K',num2str(params.K,'%.4d'),'_lambda',strrep(num2str(params.lambda,'%.2f'),'.',''),'.mat');
load(ODL_FILENAME);

D_HR = D_HR./repmat(sqrt(sum(D_HR.^2, 1)), dim_h, 1);
D_HR=sqrt(dim_h)/sqrt(dim_l)*D_HR; % scale the HR dictionary
D_LR = D_LR./repmat(sqrt(sum(D_LR.^2, 1)), dim_l, 1);

[~,ids]=sort(sum(CoefMatrix.^2,2),'descend');
fig_HR=figure(); set(fig_HR, 'Color', 'w'); ImD=displayPatches(D_HR(:,ids(1:10:end))); axis off;

fig_LR=figure(); set(fig_HR, 'Color', 'w'); ImD2=displayPatches(D_LR(:,ids(1:10:end))); axis off;


%% GRID OF FIRST NODES IN LR AND HR PATCHES (enlarged version)
gridz_l = 1: patchsize_l-overlap : Nl+1;
gridy_l = 1: patchsize_l-overlap : Nl+1;
gridz_h = 1: space_spacing*(patchsize_l-overlap) : Nh+1;
gridy_h = 1: space_spacing*(patchsize_l-overlap) : Nh+1;  
grids={gridz_l,gridy_l,gridz_h,gridy_h};
sizes={Nl+patchsize_l,Nl+patchsize_l,Nh+patchsize_h,Nh+patchsize_h};

%% SR
params.mode=2; % min_{D in C} (1/n) sum_{i=1}^n (1/2)||x_i-Dalpha_i||_2^2 + lambda||alpha_i||_1 + lambda_2||alpha_i||_2^2
params.lambda = 1e-6; 
params.lambda2 = 0;

t=10; planeid=2;
nc = netcdf(filename_ref,'r');
X_HR_ref=nc{'velocity_x'}(t,:,:,planeid);
X_LR=X_HR_ref(1:space_spacing:end,1:space_spacing:end);
close(nc);

fprintf('Increasing resolution...\n');
        
% LR field
left=patchsize_h/2; right = patchsize_h/2-(Nh-gridz_LR(1,end));
bottom=patchsize_h/2; top = patchsize_h/2-(Nh-gridy_LR(end,1));

[gridz_HS_enlarged, gridy_HS_enlarged] = meshgrid(1:Nh+left+right,1:Nh+bottom+top);
gridz_LS_enlarged=gridz_HS_enlarged(1:space_spacing:end,1:space_spacing:end);
gridy_LS_enlarged=gridy_HS_enlarged(1:space_spacing:end,1:space_spacing:end);

X_HR_ref_enlarged=enlarge_2D(X_HR_ref,left, right, bottom, top);
X_LR_enlarged = X_HR_ref_enlarged(1:space_spacing:end,1:space_spacing:end);
X_HR_interp_enlarged=interp2(gridz_LS_enlarged, gridy_LS_enlarged, X_LR_enlarged, gridz_HS_enlarged, gridy_HS_enlarged,'spline');
X_HR_interp=X_HR_interp_enlarged(bottom+1:bottom+Nh,left+1:left+Nh);

% SR
start=tic();
X_HR_rec_enlarged = SR_nguyen(X_LR_enlarged, D_HR, D_LR, params,'HRLR', grids, sizes);
X_HR_rec = X_HR_rec_enlarged(patchsize_h/2+1:end-patchsize_h/2,patchsize_h/2+1:end-patchsize_h/2);

NRMSE1=sqrt(sum((X_HR_interp(:)-X_HR_ref(:)).^2))/sqrt(sum(X_HR_ref(:).^2));
NRMSE2=sqrt(sum((X_HR_rec(:)-X_HR_ref(:)).^2))/sqrt(sum(X_HR_ref(:).^2));
fprintf(['SR snapshot t=',num2str(t,'%.3d'),' in ',num2str(toc(start)/60,'%.3f'),' minute(s), improves NRMSE by ', num2str((1-NRMSE2/NRMSE1)*100,'%.2f'),' percents\n']);


%% PLOT
fsize=16;
fname='CMU Serif';
set(0,'DefaultAxesFontName', fname)
set(0,'DefaultAxesFontSize', fsize)
set(0,'DefaultTextFontname', fname)
set(0,'DefaultTextFontSize', fsize)

fig=figure();
set(fig, 'Position', [200 200 700 650]);
set(fig,'color','w')

ax(1)=subplot(2,2,1,'position',[0.05 0.525 0.425 0.425]); 
imagesc(X_LR); caxis([-4 4]); 
set(gca,'YDir','normal'); axis off; axis equal;
title('LR')

ax(2)=axes('position',[0.05 0.05 0.425 0.425]); 
imagesc(X_HR_interp); caxis([-4 4]); set(gca,'YDir','normal'); axis off; axis equal;
title('Interpolation')

ax(3)=axes('position',[0.5 0.525 0.425 0.425]); 
imagesc(X_HR_rec); caxis([-4 4]); set(gca,'YDir','normal'); axis off; axis equal;
title('SR')

ax(4)=axes('position',[0.5 0.05 0.425 0.425]); 
imagesc(X_HR_ref); caxis([-4 4]); set(gca,'YDir','normal'); axis off; axis equal;
title('Reference')

% export_fig('reconstruction_patchsize08_patchesHR_patchesLRprefilbicubic','-eps');
% close();
