clear all; close all; clc;

%% INITIAL PARAMS
space_spacing=3;
time_spacing=4;

patchsize_l = 4; % 8x8 LR patches
patchsize_h = patchsize_l*space_spacing; % size of HR patches
dim_patch=patchsize_h^2;

num_patch=16*16; % number of patches extracted from each plane

filename_ref='/data/ISOTROPIC/data/data_downsampled4.nc';
nc = netcdf(filename_ref,'r');
Nt = nc('Nt').itsDimsize;
Nh = nc('Nx').itsDimsize;
Nl=Nh/space_spacing;

u_HR_all=nc{'velocity_x'}(1:Nt,1:Nh,1:Nh,1:space_spacing:Nh);
close(nc)

LTHS_idt=1:time_spacing:Nh;

%% Random indices at HR or LR
xh = randperm(Nh-2*patchsize_h-1) + patchsize_h; % starting point of HR patches, avoid border by convolution
yh = randperm(Nh-2*patchsize_h-1) + patchsize_h;

xh=xh(1:sqrt(num_patch));
yh=yh(1:sqrt(num_patch));

[Xh,Yh] = meshgrid(xh,yh);
xrow_h = Xh(:); ycol_h = Yh(:); 

%% FILTERS
O = zeros(1, space_spacing-1);
% O=0;
hf1 = [-1,O,1];
vf1 = hf1';
hf2 = [1,O,-2,O,1];
vf2 = hf2';

%% EXTRACT PATCHES
patches_af_all=zeros(dim_patch,num_patch*Nt*numel(LTHS_idt));
patches_lf_all = zeros(dim_patch,num_patch*Nt*numel(LTHS_idt));
patches_hf_all = zeros(dim_patch,num_patch*Nt*numel(LTHS_idt));
patches_lf_features_all = zeros(4*dim_patch,num_patch*Nt*numel(LTHS_idt));
index=1;
for t=1:Nt
    t
    for i=1:numel(LTHS_idt)
        x_ref = squeeze(u_HR_all(t,:,:,i));       
        x_lf = filter_2D(x_ref, space_spacing, space_spacing);
        x_hf = x_ref - x_lf;

        x_lf_Fea(:, :, 1) = conv2(x_lf, hf1, 'same');
        x_lf_Fea(:, :, 2) = conv2(x_lf, vf1, 'same');
        x_lf_Fea(:, :, 3) = conv2(x_lf,hf2,'same');
        x_lf_Fea(:, :, 4) = conv2(x_lf,vf2,'same');
     
        % Extract random patches: HR and corresponding subsampled LR
        patches_af = zeros(dim_patch, length(xrow_h));
        patches_lf = zeros(dim_patch, length(xrow_h));
        patches_lf_features = zeros(4*dim_patch, length(xrow_h));
        patches_hf = zeros(dim_patch, length(xrow_h));

        for ii = 1:num_patch
            row_h = xrow_h(ii); col_h = ycol_h(ii);

            % all frequency
            patch_af = x_ref(row_h:row_h+patchsize_h-1,col_h:col_h+patchsize_h-1);
            patches_af(:,ii)=patch_af(:);

            % low frequency
            patch_lf = x_lf(row_h:row_h+patchsize_h-1,col_h:col_h+patchsize_h-1);
            patches_lf(:,ii) = patch_lf(:);  
            
            % low frequency
            patch_lf_feature1 = x_lf_Fea(row_h:row_h+patchsize_h-1,col_h:col_h+patchsize_h-1,1);
            patch_lf_feature2 = x_lf_Fea(row_h:row_h+patchsize_h-1,col_h:col_h+patchsize_h-1,2);
            patch_lf_feature3 = x_lf_Fea(row_h:row_h+patchsize_h-1,col_h:col_h+patchsize_h-1,3);
            patch_lf_feature4 = x_lf_Fea(row_h:row_h+patchsize_h-1,col_h:col_h+patchsize_h-1,4);
            patches_lf_features(:,ii) = [patch_lf_feature1(:); patch_lf_feature2(:); patch_lf_feature3(:); patch_lf_feature4(:)];

            % high frequency
            patch_hf = x_hf(row_h:row_h+patchsize_h-1,col_h:col_h+patchsize_h-1);
            patches_hf(:,ii) = patch_hf(:);                
        end
        
        patches_af_all(:,index:index+num_patch-1)= patches_af;
        patches_lf_all(:,index:index+num_patch-1) = patches_lf;
        patches_hf_all(:,index:index+num_patch-1) = patches_hf;
        patches_lf_features_all(:,index:index+num_patch-1) = patches_lf_features;
        index=index+num_patch;
    end
end
PATCHES_FILENAME=strcat('/data/ISOTROPIC/dictionary_learning/ideal/space_ratio_03/trainingpatches_sspacing',num2str(space_spacing,'%.1d'),'_tspacing',num2str(time_spacing,'%.1d'),'_patchsize',num2str(patchsize_l,'%.2d'),'_numpatch',num2str(Nt*numel(LTHS_idt)*num_patch,'%.6d'),'.mat');
save(PATCHES_FILENAME,'patches_af_all','patches_lf_all','patches_hf_all','patches_lf_features_all');
