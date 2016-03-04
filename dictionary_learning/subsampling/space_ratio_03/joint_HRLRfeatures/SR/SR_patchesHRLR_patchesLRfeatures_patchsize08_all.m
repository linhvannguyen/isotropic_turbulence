clear all; close all; clc;
path(path,'/home/nguyen/Documents/Codes/isotropic/DictionaryLearning/final/funcs')

%% INITIAL PARAMS
Nh=96; % Number of grid-point in HR field
scale=3;
Nl=Nh/scale;

size_l=8; 
patchsize_h = size_l*scale; % size of HR patches
patchsize_l = patchsize_h; % 4x4 LR patches

dim_h=patchsize_h^2;
dim_l=4*patchsize_l^2;
dim_l_pca=81;

nc = netcdf('/data/ISOTROPIC/NS3D_DNS_FORCE_384/TRAININGPLANES_downsampled4_HR.nc','r');
Nx = nc('Nt').itsDimsize;
close(nc)

overlap=7; % at LR

params.lambda=0.2; % params.lambda=1*1e-5, overlap=7
params.K=2*(dim_h+dim_l_pca);

%% LOAD DICTIONARY 
ODL_FILENAME=strcat('/data/ISOTROPIC/DICTIONARYLEARNING/final/downsample4_3/joint_HRLRfeatures/dictionaries/patchsize08/DICTIONARY_patchesHRLR_patchesLRfeatures_K',num2str(params.K,'%.4d'),'_lambda',strrep(num2str(params.lambda,'%.2f'),'.',''),'.mat');
load(ODL_FILENAME);

D_HR = D_HR./repmat(sqrt(sum(D_HR.^2, 1)), dim_h, 1);
D_HR=sqrt(dim_h/dim_l_pca)*D_HR; % scale the HR dictionary (very important!!!)
D_LR = D_LR./repmat(sqrt(sum(D_LR.^2, 1)), dim_l_pca, 1);

%% PLOT DICTIONARY
% [~,ids]=sort(sum(CoefMatrix.^2,2),'descend');
% fig_HR=figure(); set(fig_HR, 'Color', 'w'); displayPatches(D_HR(:,ids(1:2:end))); axis off;
% DL_FILENAME=strcat('DICTIONARY_patchsize08_patchesHRLR_patchesLRfeatures_K',num2str(params.K,'%.4d'),'_lambda',strrep(num2str(params.lambda,'%.2f'),'.',''));
% export_fig(DL_FILENAME,'-eps');
% close();

%% FILTERS
% O = zeros(1, scale-1);
O=0;
hf1 = [-1,O,1]; vf1 = hf1';
hf2 = [1,O,-2,O,1]; vf2 = hf2';

%% GRID OF FIRST NODES IN LR AND HR PATCHES (enlarged version)
gridz_l = 1: size_l-overlap : Nl+1;
gridy_l = 1: size_l-overlap : Nl+1;
gridz_h = 1: scale*(size_l-overlap) : Nh+1;
gridy_h = 1: scale*(size_l-overlap) : Nh+1;  
grids={gridz_h,gridy_h};
sizes={Nl+size_l,Nl+size_l,Nh+patchsize_h,Nh+patchsize_h};

%% SR
params.mode=2; % min_{D in C} (1/n) sum_{i=1}^n (1/2)||x_i-Dalpha_i||_2^2 + lambda||alpha_i||_1 + lambda_2||alpha_i||_2^2
params.lambda=1*1e-5; 
params.lambda2=0;
scaling_factor=1; % is this sqrt(patchsize_h*patchsize_h/dim_l_pca) (currently dim_l_pca=100; best: 2.15)


fprintf('Increasing resolution...\n');

ttds=1:200;
% ttds=201:592;

nc1 = netcdf('/data/ISOTROPIC/DICTIONARYLEARNING/final/downsample4_3/joint_HRLRfeatures/SRrec/SR_patchesLRfeatures_patchsize08_all_1.nc','clobber');
% nc1 = netcdf('/data/ISOTROPIC/DICTIONARYLEARNING/final/downsample4_3/joint_HRLRfeatures/SRrec/SR_patchesLRfeatures_patchsize08_all_2.nc','clobber');

nc1('Nt')=0; 
nc1('Ny')=Nh;
nc1('Nz')=Nh;
nc1{'X_HR_rec'}=ncfloat('Nt','Ny','Nz');

nc2 = netcdf('/data/ISOTROPIC/NS3D_DNS_FORCE_384/RECONSTRUCTINGPLANES_downsampled4_HR.nc','r');

for t=1:numel(ttds) 
    % Load data
    X_HR_org=nc2{'velocity_x_HR'}(ttds(t),:,:);
    
    % LR field
    X_LR=resize_nguyen(X_HR_org, 1/scale,'bicubic');
    X_LR_enlarged = [X_LR(Nl-size_l/2+1:Nl,Nl-size_l/2+1:Nl),X_LR(Nl-size_l/2+1:Nl,:),X_LR(Nl-size_l/2+1:Nl,1:size_l/2);...
    X_LR(:,Nl-size_l/2+1:Nl),X_LR,X_LR(:,1:size_l/2);...
    X_LR(1:size_l/2,Nl-size_l/2+1:Nl),X_LR(1:size_l/2,:),X_LR(1:size_l/2,1:size_l/2)];

    X_LR_interp=resize_nguyen(X_LR_enlarged, scale,'bicubic');

    % SR
    start=tic();

    X_HR_rec_smallscales=SR_nguyen(X_LR_interp, D_HR, D_LR, params,'FLUC_FEATURES_INTERP',grids, sizes, {hf1, vf1, hf2, vf2}, V_pca,scaling_factor);
    X_HR_rec = X_LR_interp+X_HR_rec_smallscales;

    % REMOVE BOLDER (translation is not complete near bolder)
    X_HR_rec=X_HR_rec(patchsize_h/2+1:end-patchsize_h/2,patchsize_h/2+1:end-patchsize_h/2);
    X_LR_interp=X_LR_interp(patchsize_h/2+1:end-patchsize_h/2,patchsize_h/2+1:end-patchsize_h/2);
    
    % SAVE FILE
    nc1{'X_HR_rec'}(t,:,:)=X_HR_rec;
       
    NRMSE1=sqrt(sum((X_LR_interp(:)-X_HR_org(:)).^2))/sqrt(sum(X_HR_org(:).^2));
    NRMSE2=sqrt(sum((X_HR_rec(:)-X_HR_org(:)).^2))/sqrt(sum(X_HR_org(:).^2));
    fprintf(['SR snapshot t=',num2str(t,'%.3d'),' in ',num2str(toc(start)/60,'%.3f'),' minute(s), improves NRMSE by ', num2str((1-NRMSE2/NRMSE1)*100,'%.2f'),' percents\n']);
end
close(nc2); close(nc1);
