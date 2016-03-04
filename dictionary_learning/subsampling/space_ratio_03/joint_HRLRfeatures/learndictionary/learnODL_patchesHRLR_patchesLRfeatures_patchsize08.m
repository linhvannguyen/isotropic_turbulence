% clear all; close all; clc;

%% INITIAL PARAMS
space_spacing=3;
time_spacing=4;

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

size_l=8; 
patchsize_l = size_l; 
patchsize_h = patchsize_l*space_spacing;
dim_h=patchsize_h^2;
dim_l=patchsize_l^2;

left=patchsize_h/2; right = patchsize_h/2-(Nh-gridz_LR(1,end));
bottom=patchsize_h/2; top = patchsize_h/2-(Nh-gridy_LR(end,1));

[gridz_HS_enlarged, gridy_HS_enlarged] = meshgrid(1:Nh+left+right,1:Nh+bottom+top);
gridz_LS_enlarged=gridz_HS_enlarged(1:space_spacing:end,1:space_spacing:end);
gridy_LS_enlarged=gridy_HS_enlarged(1:space_spacing:end,1:space_spacing:end);

%% Load patches
% size_l=8; 
% patchsize_h = size_l*space_spacing; % size of HR patches
% patchsize_l = patchsize_h; % 4x4 LR patches
% 
% num_patch=8*8; % number of patches 
% dim_h=patchsize_h^2;
% dim_l=4*patchsize_l^2;

% PATCHES_FILENAME=strcat('/data/ISOTROPIC/dictionary_learning/subsampling/space_ratio_03/trainingpatches_couplefeatures_spaceratio',num2str(space_spacing,'%.1d'),'_timeratio',num2str(time_spacing,'%.1d'),'_patchsize',num2str(size_l,'%.2d'),'_numpatch',num2str(Nt*numel(LTHS_idt)*num_patch,'%.6d'),'.mat');
% load(PATCHES_FILENAME, 'patches_HRLR_all','patches_LR_features_all');
% 
% % SUBSTRACT MEAN AND PRE-NORMALIZE ALL PATCHES
% start=tic();
% m_LR=mean(patches_LR_features_all,1);
% for i=1:size(patches_LR_features_all,2)
%     patches_LR_features_all(:,i) = patches_LR_features_all(:,i) - m_LR(i);
% end
% norm_LR=sqrt(sum(patches_LR_features_all.^2,1));
% for i=1:dim_l
%     patches_LR_features_all(i,:) = patches_LR_features_all(i,:)./norm_LR;
% end
% 
% m_HR=mean(patches_HRLR_all,1);
% for i=1:size(patches_HRLR_all,2)
%     patches_HRLR_all(:,i) = patches_HRLR_all(:,i) - m_HR(i);
% end
% norm_HR=sqrt(sum(patches_HRLR_all.^2,1));
% for i=1:dim_h
%     patches_HRLR_all(i,:) = patches_HRLR_all(i,:)./norm_HR;
% end
% fprintf(['Processing data in ',num2str(toc(start),'%.3f'),' seconds \n']);


%% Filter
O = zeros(1, space_spacing-1);
% O=0;
hf1 = [-1,O,1];
vf1 = hf1';
hf2 = [1,O,-2,O,1];
vf2 = hf2';

filters={hf1,vf1,hf2,vf2};
npatches={Nl-size_l+1,Nl-size_l+1}; % max: Nl-size_l+1
npatches_all=npatches{1}*npatches{2};

%% PCA 
% start_all=tic();
% fprintf('\n=========== STEP 1: PCA =========== \n');
% 
% start=tic();
% C = zeros(4*dim_h,4*dim_h);
% 
% for t=1:Nt
%     fprintf(['\n Proceeding block: ', num2str(t,'%.2d'),' \n']);
%     index=1;
%     patches_features_all=zeros(4*dim_h, numel(LTHS_idt)*npatches_all);
%     for i=1:numel(LTHS_idt)
%         X_HR_ref=squeeze(u_HR_all(t,:,:,i)); 
% 
%         X_HR_ref_enlarged=enlarge_2D(X_HR_ref,left, right, bottom, top);
%         X_LR = X_HR_ref_enlarged(1:space_spacing:end,1:space_spacing:end);
%         X_HR_interp=interp2(gridz_LS_enlarged, gridy_LS_enlarged, X_LR, gridz_HS_enlarged, gridy_HS_enlarged,'spline');
%         X_feas=extract_features(X_HR_interp,filters);
%         
%         X_feas=X_feas(bottom+1:bottom+Nh,left+1:left+Nh,:);
%         X_HR_interp=X_HR_interp(bottom+1:bottom+Nh,left+1:left+Nh);
% 
%         patches_features = extract_random_patches(X_feas, size_l, space_spacing, npatches);
%         patches_features_all(:,index:index+npatches_all-1) = patches_features; 
%         index=index+npatches_all;
%     end
%     
%     C = C + patches_features_all * patches_features_all'; 
%     clearvars patches_features_all;
% end 
% 
% % PCA
% [V, D] = eig(C);
% D = diag(D);
% D = cumsum(D) / sum(D);
% 
% k = find(D >= 1e-4, 1); % ignore 0.01% energy
% dim_l_pca=numel(D)-k;
% 
% 
% V_pca = V(:, end-dim_l_pca+1:end); % ignore 0.23% energy
% 
% save pca_features.mat V_pca dim_l_pca
% 
% fprintf(['\nDimension of truncated data: ',num2str(dim_l_pca,'%.3d'),' seconds \n']);
% 
% fprintf(['\nFinish step 1 in ',num2str(toc(start),'%.3f'),' seconds \n']);

%% Extract patches and online training
load pca_features.mat;
fprintf('\n=========== STEP 2: ONLINE TRAINING =========== \n');

params.mode=2; % min_{D in C} (1/n) sum_{i=1}^n (1/2)||x_i-Dalpha_i||_2^2 + lambda||alpha_i||_1 + lambda_2||alpha_i||_2^2
params.lambda=0.2;
params.lambda2=0;
params.K=2*(dim_h+dim_l_pca);
params.numThreads=4; % number of threads
params.iter=10;  % max number of iterations.

for t=1:Nt

    fprintf(['\nTraining on block: ', num2str(t,'%.2d'),' \n']);
    
    index=1;
    patches_HRLR_all=zeros(dim_h,numel(LTHS_idt)*npatches_all);
    patches_features_all=zeros(4*dim_h,numel(LTHS_idt)*npatches_all);

    for i=1:numel(LTHS_idt)
        X_HR_ref=squeeze(u_HR_all(t,:,:,i)); 

        X_HR_ref_enlarged=enlarge_2D(X_HR_ref,left, right, bottom, top);
        X_LR = X_HR_ref_enlarged(1:space_spacing:end,1:space_spacing:end);
        X_HR_interp=interp2(gridz_LS_enlarged, gridy_LS_enlarged, X_LR, gridz_HS_enlarged, gridy_HS_enlarged,'spline');
        X_feas=extract_features(X_HR_interp,filters);
        
        X_feas=X_feas(bottom+1:bottom+Nh,left+1:left+Nh,:);
        X_HR_interp=X_HR_interp(bottom+1:bottom+Nh,left+1:left+Nh);

        patches_HRLR = extract_random_patches(X_HR_ref-X_HR_interp, size_l, space_spacing, npatches);
        patches_features = extract_random_patches(X_feas, size_l, space_spacing, npatches);
        patches_features_all(:,index:index+npatches_all-1) = patches_features; 
        patches_HRLR_all(:,index:index+npatches_all-1)= patches_HRLR;
        
        index=index+npatches_all;
    end
    
    % Dimensionality reduction
    patches_LR_features_pca = V_pca' * patches_features_all;

    % Normalize
    patches_LR_features_pca=patches_LR_features_pca-repmat(mean(patches_LR_features_pca),[dim_l_pca 1]);
    patches_LR_features_pca=patches_LR_features_pca ./ repmat(sqrt(sum(patches_LR_features_pca.^2)),[dim_l_pca 1]); 

    patches_HRLR_all=patches_HRLR_all-repmat(mean(patches_HRLR_all),[dim_h 1]);
    patches_HRLR_all=patches_HRLR_all ./ repmat(sqrt(sum(patches_HRLR_all.^2)),[dim_h 1]);   

    % join patches
    patches_all = [(1/sqrt(dim_h))*patches_HRLR_all; (1/sqrt(dim_l_pca))*patches_LR_features_pca];
    patches_all=patches_all-repmat(mean(patches_all),[dim_h+dim_l_pca 1]);
    patches_all=patches_all ./ repmat(sqrt(sum(patches_all.^2)),[dim_h+dim_l_pca 1]);   
    
    start = tic();
    % Online dictionary learning
    if (t==1)
        [D,model] = mexTrainDL(patches_all,params);
    else
        [D,model] = mexTrainDL(patches_all,params, model); % retrain the model
    end
    
    CoefMatrix=mexLasso(patches_all,D,params);
    R=mean(0.5*sum((patches_all-D*CoefMatrix).^2)+params.lambda*sum(abs(CoefMatrix)));
 
    fprintf(['\nFinish training on block: ',num2str(t,'%.2d'), ' in ', num2str(toc(start),'%.3f'),' seconds, objective function: ', num2str(R,'%.3f'), ' \n']);   
end

fprintf(['\n=========== COMPLETED (', num2str(toc(start_all)/3600,'%.3f'),' hours) =========== \n']);

%% ONLINE DICTIONARY LEARNING
ODL_FILENAME=strcat('/data/ISOTROPIC/dictionary_learning/subsampling/space_ratio_03/DICTIONARY_couplefeatures_patchesHR_patchesLR_joint_K',num2str(params.K,'%.4d'),'_lambda',strrep(num2str(params.lambda,'%.2f'),'.',''),'.mat');
save(ODL_FILENAME,'D','V_pca');