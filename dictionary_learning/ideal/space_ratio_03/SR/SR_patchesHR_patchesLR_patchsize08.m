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

u_all=nc{'velocity_x'}(1:Nt,1:Nh,1:Nh,1:space_spacing:Nh);
close(nc)

LTHS_idt=1:time_spacing:Nh;

%% LOAD DICTIONARY 
overlap=7; % at LR

params.lambda=0.05; % params.lambda=0.001, overlap=9
params.K=2*dim_patch; 
load /data/ISOTROPIC/dictionary_learning/ideal/space_ratio_03/DICTIONARY_coupleAFLF_K1152_lambda005.mat
D_AF = D_AF./repmat(sqrt(sum(D_AF.^2, 1)), dim_patch, 1);
D_LF = D_LF./repmat(sqrt(sum(D_LF.^2, 1)), dim_patch, 1);

[~,ids]=sort(sum(CoefMatrix.^2,2),'descend');
fig_AF=figure(); set(fig_AF, 'Color', 'w'); ImD=displayPatches(D_AF(:,ids(1:10:end))); axis off;
fig_LF=figure(); set(fig_LF, 'Color', 'w'); ImD2=displayPatches(D_LF(:,ids(1:10:end))); axis off;

%% GRID OF FIRST NODES IN LR AND HR PATCHES (enlarged version)
gridz_l = 1: patchsize_l-overlap : Nl+1;
gridy_l = 1: patchsize_l-overlap : Nl+1;
gridz_h = 1: space_spacing*(patchsize_l-overlap) : Nh+1;
gridy_h = 1: space_spacing*(patchsize_l-overlap) : Nh+1;  
grids={gridz_h,gridy_h};
sizes={Nh+patchsize_h,Nh+patchsize_h};

%% SR
params.mode=2; % min_{D in C} (1/n) sum_{i=1}^n (1/2)||x_i-Dalpha_i||_2^2 + lambda||alpha_i||_1 + lambda_2||alpha_i||_2^2
params.lambda = 1e-4; 
params.lambda2 = 0;

t=10; planeid=3;
nc = netcdf(filename_ref,'r');
xref=nc{'velocity_x'}(t,:,:,planeid);
close(nc);

fprintf('Increasing resolution...\n');
        
% LR field
left=patchsize_h/2; right = patchsize_h/2;
bottom=patchsize_h/2; top = patchsize_h/2;

[gridz_HS_enlarged, gridy_HS_enlarged] = meshgrid(1:Nh+left+right,1:Nh+bottom+top);
gridz_LS_enlarged=gridz_HS_enlarged(1:space_spacing:end,1:space_spacing:end);
gridy_LS_enlarged=gridy_HS_enlarged(1:space_spacing:end,1:space_spacing:end);

xref_enlarged=enlarge_2D(xref,left, right, bottom, top);
xref_LF_enlarged = filter_2D(xref_enlarged, space_spacing, space_spacing);
xref_HF_enlarged=xref_enlarged-xref_LF_enlarged;
xref_HF=xref_HF_enlarged(bottom+1:bottom+Nh,left+1:left+Nh);
xref_LF=xref_LF_enlarged(bottom+1:bottom+Nh,left+1:left+Nh);

% SR
start=tic();
xrec_enlarged = SR_nguyen(xref_LF_enlarged, D_AF, D_LF, params,'LFHF', grids, sizes);
xrec = xrec_enlarged(patchsize_h/2+1:end-patchsize_h/2,patchsize_h/2+1:end-patchsize_h/2);

NRMSE1=sqrt(sum((xref_LF(:)-xref(:)).^2))/sqrt(sum(xref(:).^2));
NRMSE2=sqrt(sum((xrec(:)-xref(:)).^2))/sqrt(sum(xref(:).^2));
fprintf(['SR snapshot t=',num2str(t,'%.3d'),' in ',num2str(toc(start)/60,'%.3f'),' minute(s), improves NRMSE by ', num2str((1-NRMSE2/NRMSE1)*100,'%.2f'),' percents\n']);

xrec_HF = xref-xrec;
NRMSE3=sqrt(sum((xrec_HF(:)-xref_HF(:)).^2))/sqrt(sum(xref_HF(:).^2));

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
imagesc(xref_LF); caxis([-4 4]); 
set(gca,'YDir','normal'); axis off; axis equal;
title('LR')

ax(2)=axes('position',[0.05 0.05 0.425 0.425]); 
imagesc(xref_HF); caxis([-1 1]); set(gca,'YDir','normal'); axis off; axis equal;
title('Interpolation')

ax(3)=axes('position',[0.5 0.525 0.425 0.425]); 
imagesc(xrec); caxis([-4 4]); set(gca,'YDir','normal'); axis off; axis equal;
title('SR')

ax(4)=axes('position',[0.5 0.05 0.425 0.425]); 
imagesc(xref-xrec); caxis([-1 1]); set(gca,'YDir','normal'); axis off; axis equal;
title('Reference')

% export_fig('reconstruction_patchsize08_patchesHR_patchesLRprefilbicubic','-eps');
% close();
