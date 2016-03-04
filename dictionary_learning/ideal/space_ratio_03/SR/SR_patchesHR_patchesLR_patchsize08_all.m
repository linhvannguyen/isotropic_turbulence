clear all; close all; clc;
path(path,'/home/nguyen/Documents/Codes/isotropic/DictionaryLearning/final/funcs')

%% INITIAL PARAMS
Nh=96; % Number of grid-point in HR field
scale=3;
Nl=Nh/scale;
patchsize_l = 8; % 8x8 LR patches
patchsize_h = patchsize_l*scale; % size of HR patches

dim_h=patchsize_h^2;
dim_l=patchsize_l^2;

nc2 = netcdf('/data/ISOTROPIC/NS3D_DNS_FORCE_384/TRAININGPLANES_downsampled4_HR.nc','r');
Nx = nc2('Nt').itsDimsize;
close(nc2)

overlap=7; % at LR
num_patch=113664; % total number of patches

params.lambda=0.2; % params.lambda=0.001, overlap=9
params.K=2*(dim_h+dim_l);  

%% LOAD DICTIONARY 
ODL_FILENAME=strcat('/data/ISOTROPIC/DICTIONARYLEARNING/final/downsample4_3/joint_HRLR/dictionaries/patchsize08/DICTIONARY_patchesHR_patchesLR_joint_K',num2str(params.K,'%.4d'),'_lambda',strrep(num2str(params.lambda,'%.2f'),'.',''),'.mat');
load(ODL_FILENAME);

D_HR = D_HR./repmat(sqrt(sum(D_HR.^2, 1)), dim_h, 1);
D_HR=sqrt(dim_h)/sqrt(dim_l)*D_HR; % scale the HR dictionary
D_LR = D_LR./repmat(sqrt(sum(D_LR.^2, 1)), dim_l, 1);

%% GRID OF FIRST NODES IN LR AND HR PATCHES (enlarged version)
gridz_l = 1: patchsize_l-overlap : Nl+1;
gridy_l = 1: patchsize_l-overlap : Nl+1;
gridz_h = 1: scale*(patchsize_l-overlap) : Nh+1;
gridy_h = 1: scale*(patchsize_l-overlap) : Nh+1;  
grids={gridz_l,gridy_l,gridz_h,gridy_h};
sizes={Nl+patchsize_l,Nl+patchsize_l,Nh+patchsize_h,Nh+patchsize_h};

%% SR
params.mode=2; % min_{D in C} (1/n) sum_{i=1}^n (1/2)||x_i-Dalpha_i||_2^2 + lambda||alpha_i||_1 + lambda_2||alpha_i||_2^2
params.lambda=1*1e-6; 
params.lambda2=0;

fprintf('Increasing resolution...\n');

nc1 = netcdf('/data/ISOTROPIC/DICTIONARYLEARNING/final/downsample4_3/joint_HRLR/SRrec/SR_patchesHRLR_patchsize08_all.nc','clobber');
nc1('Nt')=0; 
nc1('Ny')=Nh;
nc1('Nz')=Nh;
nc1{'X_HR_rec'}=ncfloat('Nt','Ny','Nz');


ttds=1:592;
nc2 = netcdf('/data/ISOTROPIC/NS3D_DNS_FORCE_384/RECONSTRUCTINGPLANES_downsampled4_HR.nc','r');
for t=1:numel(ttds) 
    % Load data
    X_HR_org=nc2{'velocity_x_HR'}(ttds(t),:,:);
    
    % LR field
    X_LR=resize_nguyen(X_HR_org, 1/scale,'bicubic');
    X_LR_enlarged = [X_LR(Nl-patchsize_l/2+1:Nl,Nl-patchsize_l/2+1:Nl),X_LR(Nl-patchsize_l/2+1:Nl,:),X_LR(Nl-patchsize_l/2+1:Nl,1:patchsize_l/2);...
        X_LR(:,Nl-patchsize_l/2+1:Nl),X_LR,X_LR(:,1:patchsize_l/2);...
        X_LR(1:patchsize_l/2,Nl-patchsize_l/2+1:Nl),X_LR(1:patchsize_l/2,:),X_LR(1:patchsize_l/2,1:patchsize_l/2)];

    X_LR_interp=resize_nguyen(X_LR_enlarged, scale,'bicubic');
    X_LR_interp=X_LR_interp(patchsize_h/2+1:end-patchsize_h/2,patchsize_h/2+1:end-patchsize_h/2);
    
    % SR
    start=tic();
    X_HR_rec = SR_nguyen(X_LR_enlarged, D_HR, D_LR, params,'HRLR', grids, sizes);
    X_HR_rec = X_HR_rec(patchsize_h/2+1:end-patchsize_h/2,patchsize_h/2+1:end-patchsize_h/2);
    
    nc1{'X_HR_rec'}(t,:,:)=X_HR_rec;
    
    NRMSE1=sqrt(sum((X_LR_interp(:)-X_HR_org(:)).^2))/sqrt(sum(X_HR_org(:).^2));
    NRMSE2=sqrt(sum((X_HR_rec(:)-X_HR_org(:)).^2))/sqrt(sum(X_HR_org(:).^2));
    fprintf(['SR snapshot t=',num2str(t,'%.3d'),' in ',num2str(toc(start)/60,'%.3f'),' minute(s), improves NRMSE by ', num2str((1-NRMSE2/NRMSE1)*100,'%.2f'),' percents\n']);
end
close(nc2); close(nc1)

