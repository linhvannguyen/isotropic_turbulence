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

filename_ref='/data/ISOTROPIC/data/data_downsampled4.nc';
nc = netcdf(filename_ref,'r');
Nt = nc('Nt').itsDimsize;
Nh = nc('Nx').itsDimsize;
Nl=Nh/space_spacing;
close(nc)

midplane_ids=3:time_spacing:Nh;

overlap=7; % at LR

params_train.lambda=0.2; % params.lambda=1*1e-5, overlap=7
params_train.K=2*(dim_h+dim_l_pca);

LTHS_idt=1:time_spacing:Nh;
[gridz_HR, gridy_HR] = meshgrid(1:Nh,1:Nh);
gridz_LR = gridz_HR(1:space_spacing:Nh,1:space_spacing:Nh);
gridy_LR = gridy_HR(1:space_spacing:Nh,1:space_spacing:Nh);

%% LOAD DICTIONARY 
ODL_FILENAME=strcat('/data/ISOTROPIC/dictionary_learning/downsampling/space_ratio_03/DICTIONARY_couplefeatures_patchesHR_patchesLR_joint_K',num2str(params_train.K,'%.4d'),'_lambda',strrep(num2str(params_train.lambda,'%.2f'),'.',''),'.mat');
load(ODL_FILENAME);

D_HR = D_HR./repmat(sqrt(sum(D_HR.^2, 1)), dim_h, 1);
D_HR=sqrt(dim_h/dim_l_pca)*D_HR; % scale the HR dictionary (very important!!!)
D_LR = D_LR./repmat(sqrt(sum(D_LR.^2, 1)), dim_l_pca, 1);

%% FILTERS
% O = zeros(1, scale-1);
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
params_rec.mode=2; % min_{D in C} (1/n) sum_{i=1}^n (1/2)||x_i-Dalpha_i||_2^2 + lambda||alpha_i||_1 + lambda_2||alpha_i||_2^2
params_rec.lambda=1*1e-5; 
params_rec.lambda2=0;
scaling_factor=1; % is this sqrt(patchsize_h*patchsize_h/dim_l_pca) (currently dim_l_pca=100; best: 2.15)

left=patchsize_h/2; right = patchsize_h/2;
bottom=patchsize_h/2; top = patchsize_h/2;
[gridz_HS_enlarged, gridy_HS_enlarged] = meshgrid(1:Nh+left+right,1:Nh+bottom+top);
gridz_LS_enlarged=gridz_HS_enlarged(1:space_spacing:end,1:space_spacing:end);
gridy_LS_enlarged=gridy_HS_enlarged(1:space_spacing:end,1:space_spacing:end);

filenam_SR=strcat('/data/ISOTROPIC/dictionary_learning/downsampling/space_ratio_03/SR_couplefeatures_K',num2str(params_train.K,'%.4d'),'_lambda',strrep(num2str(params_train.lambda,'%.2f'),'.',''),'.nc');

fprintf('\n ********** START SUPERRESOLUTION ************ \n');
% nc = netcdf(filenam_SR,'clobber');
% nc('Nt')=0;
% nc('Nx')=numel(midplane_ids);
% nc('Ny')=Nh;
% nc('Nz')=Nh;
% nc{'Zhat_all'}=ncfloat('Nt','Nz','Ny','Nx');

nc1 = netcdf(filename_ref,'r');

for t=1:Nt
    fprintf('Superresolution: %.2d-th block.\n', t);
    Zhat_oneblock=zeros(Nh,Nh,numel(midplane_ids));
    count=0;
    for i=1:numel(midplane_ids)
        start=tic();
        plane_id = midplane_ids(i);

        X_HR_org=nc1{'velocity_x'}(t,1:Nh,1:Nh,plane_id);
        X_HR_org_enlarged=enlarge_2D(X_HR_org,left, right, bottom, top);

        % LR field
        X_LR_enlarged=resize_nguyen(X_HR_org_enlarged, 1/space_spacing,'bicubic');

        X_LR_interp_enlarged=resize_nguyen(X_LR_enlarged, space_spacing,'bicubic');

        % SR
        Zhat_smallscales_enlarged=SR_nguyen(X_LR_interp_enlarged, D_HR, D_LR, params_rec,'FLUC_FEATURES_INTERP',grids, sizes, {hf1, vf1, hf2, vf2}, V_pca,scaling_factor);
        Zhat_enlarged = X_LR_interp_enlarged+Zhat_smallscales_enlarged;

        % REMOVE BOLDER (translation is not complete near bolder)
        Zhat=Zhat_enlarged(bottom+1:bottom+Nh,left+1:left+Nh);
        X_LR_interp=X_LR_interp_enlarged(bottom+1:bottom+Nh,left+1:left+Nh);

        NRMSE1=sqrt(sum((X_LR_interp(:)-X_HR_org(:)).^2))/sqrt(sum(X_HR_org(:).^2));
        NRMSE2=sqrt(sum((Zhat(:)-X_HR_org(:)).^2))/sqrt(sum(X_HR_org(:).^2));
        
        Zhat_oneblock(:,:,i)=Zhat;
        fprintf(['SR snapshot t=',num2str(plane_id,'%.3d'),' in ',num2str(toc(start)/60,'%.3f'),' minute(s), improves NRMSE by ', num2str((1-NRMSE2/NRMSE1)*100,'%.2f'),' percents\n']);
    end
%     nc{'Zhat_all'}(t,:,:,:) = Zhat_oneblock;
end
% close(nc); close(nc1);