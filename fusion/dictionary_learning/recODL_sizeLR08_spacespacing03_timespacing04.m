% clear all; close all; clc;

%% INITIAL PARAMS
space_spacing=3;
time_spacing=4;

% size_l=8; 
% patchsize_h = size_l*space_spacing; % size of HR patches
% dim_h=patchsize_h^2;
% params_train.K=2*dim_h;

size_l=4; 
patchsize_h = size_l*space_spacing; % size of HR patches
dim_h=patchsize_h^2;
params_train.K=2*dim_h;


% size_l=2; 
% patchsize_h = size_l*space_spacing; % size of HR patches
% dim_h=patchsize_h^2;
% params_train.K=4*dim_h;


filename_ref='/data/ISOTROPIC/data/data_downsampled4.nc';
filename_fusion=strcat('/data/ISOTROPIC/fusion/FusedData_diagCov_timespacing_',num2str(time_spacing,'%.2d'),'_spacespacing_',num2str(space_spacing,'%.2d'),'.nc');

nc1 = netcdf(filename_ref,'r');
Nt = nc1('Nt').itsDimsize;
Nh = nc1('Nx').itsDimsize;
Nl=Nh/space_spacing;
close(nc1)


overlap=size_l-1; % at LR
% params_train.lambda=0.2;
params_train.lambda=0.05;

LTHS_idt=1:time_spacing:Nh;
[gridz_HR, gridy_HR] = meshgrid(1:Nh,1:Nh);
gridz_LR = gridz_HR(1:space_spacing:Nh,1:space_spacing:Nh);
gridy_LR = gridy_HR(1:space_spacing:Nh,1:space_spacing:Nh);

%% LOAD DICTIONARY 
filename_ODL=strcat('/data/ISOTROPIC/fusion/dictionary_smallscales_sizeLR',num2str(size_l,'%.2d'),'_timespacing',num2str(time_spacing,'%.2d'),'_spacespacing',num2str(space_spacing,'%.2d'),'.mat');
load(filename_ODL);

%% PLOT DICTIONARY
[~,ids]=sort(sum(CoefMatrix.^2,2),'descend');
fig_HR=figure(); set(fig_HR, 'Color', 'w'); displayPatches(D(:,ids(1:2:end))); axis off;

%%
t=1; i=3;

nc1 = netcdf(filename_ref,'r');
x_ref=nc1{'velocity_x'}(t,1:Nh,1:Nh,i);
close(nc1)
x_ref_largescales = filter_2D(x_ref, space_spacing, space_spacing);
x_ref_smallscales = x_ref - x_ref_largescales;

nc2 = netcdf(filename_fusion,'r');
x_fusion=nc2{'Zhat_all'}(t,1:Nh,1:Nh,i);
close(nc2)
x_fusion_largescales = filter_2D(x_fusion, space_spacing, space_spacing);
x_fusion_smallscales = x_fusion - x_fusion_largescales;

% fig1=figure(); pcolor(x_ref_smallscales); shading flat; caxis([-1,1]);
% fig2=figure(); pcolor(x_fusion_smallscales); shading flat; caxis([-1,1]);

%%

[gridz_l, gridy_l]= meshgrid(-size_l+2: size_l-overlap : Nl+1,-size_l+2: size_l-overlap : Nl+1);
gridz_h = space_spacing.*gridz_l + 1;  
gridy_h = space_spacing.*gridy_l + 1;  
num_trans_y=size(gridz_l,1);
num_trans_z=size(gridz_l,2);

dZ = 0:patchsize_h-1;
dY = 0:patchsize_h-1;
dZ = reshape(dZ, [patchsize_h 1]);
dY = reshape(dY, [patchsize_h 1]);

gridz_h = repmat(gridz_h, [1 1 patchsize_h]) + permute(repmat(dZ, [1 num_trans_y num_trans_z]),[2 3 1]);
gridy_h = repmat(gridy_h, [1 1 patchsize_h]) + permute(repmat(dY, [1 num_trans_y num_trans_z]),[2 3 1]);

gridz_h(gridz_h<1) = Nh + gridz_h(gridz_h<1); 
gridy_h(gridy_h<1) = Nh + gridy_h(gridy_h<1); 
gridz_h(gridz_h>Nh) = gridz_h(gridz_h>Nh) - Nh; 
gridy_h(gridy_h>Nh) = gridy_h(gridy_h>Nh) - Nh;

% [gridz_l, gridy_l]= meshgrid(-size_l+2: size_l-overlap : Nl+1,-size_l+2: size_l-overlap : Nl+1);
% gridz_h = space_spacing.*gridz_l + 1;  
% gridy_h = space_spacing.*gridy_l + 1;  
% num_trans=size(gridz_l,1);

% [dZ,dY] = meshgrid(0:patchsize_h-1,0:patchsize_h-1);
% dZ = reshape(dZ, [1 1 patchsize_h patchsize_h]);
% dY = reshape(dY, [1 1 patchsize_h patchsize_h]);
% 
% gridz_h = repmat(gridz_h, [1 1 patchsize_h patchsize_h]) + repmat(dZ, [num_trans num_trans 1 1]);
% gridy_h = repmat(gridy_h, [1 1 patchsize_h patchsize_h]) + repmat(dY, [num_trans num_trans 1 1]);
% 
% 
% % boundary condition by periodicity
% gridz_h(gridz_h<1) = Nh + gridz_h(gridz_h<1); 
% gridy_h(gridy_h<1) = Nh + gridy_h(gridy_h<1); 
% gridz_h(gridz_h>Nh) = gridz_h(gridz_h>Nh) - Nh; 
% gridy_h(gridy_h>Nh) = gridy_h(gridy_h>Nh) - Nh;

% extract_patch = @(f)f(gridz_h + (gridy_h-1)*Nh);
% patches = extract_patch(x_fusion_smallscales);

%%
params_rec.mode=2; % min_{D in C} (1/n) sum_{i=1}^n (1/2)||x_i-Dalpha_i||_2^2 + lambda||alpha_i||_1 + lambda_2||alpha_i||_2^2
params_rec.lambda=0.05; 
params_rec.lambda2=0;
params_rec.numThreads=4;

% prepare to reconstruct
x_rec = zeros(Nh, Nh); 
cntMat = zeros(Nh, Nh); 

% collect patches
patches=zeros(patchsize_h*patchsize_h,num_trans_y*num_trans_z); 
for ii = 1:num_trans_y
    for jj = 1:num_trans_z
        zz = squeeze(gridz_h(ii,jj,:));
        yy = squeeze(gridy_h(ii,jj,:));
        patch = x_fusion_smallscales(yy, zz);
        patches(:,(ii-1)*num_trans_y+jj) = patch(:);
    end
end

m_patches=mean(patches);
norm_patches=sqrt(sum(patches.^2));

patches=patches - repmat(m_patches,[dim_h 1]);
patches=patches ./ repmat(norm_patches,[dim_h 1]); 

CoefMatrix=mexLasso(patches,D,params_rec);

rec_patches=D*CoefMatrix; 
fprintf(['NRMSE of feature reconstruction:',num2str(sqrt(sum((patches(:)-rec_patches(:)).^2))/sqrt(sum(patches(:).^2)),'%.3f'),'\n']);


% put HR patches back
for ii = 1:num_trans_y
    for jj = 1:num_trans_z
        zz = squeeze(gridz_h(ii,jj,:));
        yy = squeeze(gridy_h(ii,jj,:));
        patch = reshape(rec_patches(:,(ii-1)*num_trans_y+jj),[patchsize_h, patchsize_h]);
        patch = patch*norm_patches((ii-1)*num_trans_y+jj) + m_patches((ii-1)*num_trans_y+jj);
        x_rec(yy, zz) = x_rec(yy, zz) + patch;
        cntMat(yy, zz) = cntMat(yy, zz) + 1;
    end
end
% fill in the empty with bicubic interpolation
idx = (cntMat < 1); 
cntMat(idx) = 1; x_rec(idx) = 0;
x_rec = x_rec./cntMat;        

%%
NRMSE1=sqrt(sum((x_fusion_smallscales(:)-x_ref_smallscales(:)).^2))/sqrt(sum(x_ref_smallscales(:).^2));
NRMSE2=sqrt(sum((x_rec(:)-x_ref_smallscales(:)).^2))/sqrt(sum(x_ref_smallscales(:).^2));
fprintf(['SR improves NRMSE by ', num2str((1-NRMSE2/NRMSE1)*100,'%.2f'),' percents\n']);

fig1=figure(); pcolor(x_ref_smallscales); shading flat; caxis([-1,1]);
fig2=figure(); pcolor(x_fusion_smallscales); shading flat; caxis([-1,1]);
fig3=figure(); pcolor(x_rec); shading flat; caxis([-1,1]);

%%
x_rec_allscales = x_rec + x_fusion_largescales;
NRMSE1=sqrt(sum((x_fusion(:)-x_ref(:)).^2))/sqrt(sum(x_ref(:).^2));
NRMSE2=sqrt(sum((x_rec_allscales(:)-x_ref(:)).^2))/sqrt(sum(x_ref(:).^2));
fprintf(['SR improves NRMSE by ', num2str((1-NRMSE2/NRMSE1)*100,'%.2f'),' percents\n']);

