clear all; close all; clc;

%% INITIAL PARAMS
space_spacing=3;
time_spacing=4;
filename_ref='/data/ISOTROPIC/data/data_downsampled4.nc';
nc = netcdf(filename_ref,'r');
Nt = nc('Nt').itsDimsize;
Nh = nc('Nx').itsDimsize;

LTHS_idt=1:time_spacing:Nh;
u_HR_all=nc{'velocity_x'}(1:Nt,1:Nh,1:Nh,LTHS_idt);
close(nc)

Nl=Nh/space_spacing;

size_l=8; 
patchsize_l = size_l; 
patchsize_h = patchsize_l*space_spacing;

num_patch=8*8; % number of patches extracted from each plane
dim_h=patchsize_h^2;
dim_l=patchsize_l^2;

[gridz_HR, gridy_HR] = meshgrid(1:Nh,1:Nh);
gridz_LR = gridz_HR(1:space_spacing:Nh,1:space_spacing:Nh);
gridy_LR = gridy_HR(1:space_spacing:Nh,1:space_spacing:Nh);

%% Random indices at HR or LR

xl = randperm(Nl-2*size_l-1)+ size_l; % taking care of bolder (convolution)
yl = randperm(Nl-2*size_l-1)+ size_l; % taking care of bolder (convolution)
xl=xl(1:sqrt(num_patch));
yl=yl(1:sqrt(num_patch));

[Xl,Yl] = meshgrid(xl,yl);
 
xh=zeros(size(xl)); yh=xh; % starting point of HR patches
for i=1:numel(xl)
    % (a,b) in LR corresponding to ((a-1)*spacing+1:a*spacing,(b-1)*spacing+1:b*spacing) in HR
    xh(1,i) = (xl(i)-1)*space_spacing+1;
    yh(1,i) = (yl(i)-1)*space_spacing+1;
end
[Xh,Yh] = meshgrid(xh,yh);

xrow_l = Xl(:); ycol_l = Yl(:); 
xrow_h = Xh(:); ycol_h = Yh(:); 

%% FILTERS 
% O = zeros(1, scale-1);
O=0;
hf1 = [-1,O,1];
vf1 = hf1';
hf2 = [1,O,-2,O,1];
vf2 = hf2';

%% EXTRACT PATCHES
num_planes=Nt*numel(LTHS_idt);

patches_LR_features_all = zeros(4*dim_h,num_patch*num_planes);
patches_HRLR_all = zeros(dim_h,num_patch*num_planes);

[gridz_h,gridy_h]=meshgrid(1:Nh,1:Nh);
gridz_l=gridz_h(2:space_spacing:end,2:space_spacing:end);
gridy_l=gridy_h(2:space_spacing:end,2:space_spacing:end);

plane_id=0;
for t=1:Nt
    t
    for i=1:numel(LTHS_idt)
        plane_id=plane_id+1;
        X_HR_ref=squeeze(u_HR_all(t,:,:,i));    
        X_LR=resize_nguyen(X_HR_ref, 1/space_spacing,'bicubic');

        X_LR_interp = resize_nguyen(X_LR, space_spacing,'bicubic');

        X_LR_Fea(:, :, 1) = conv2(X_LR_interp, hf1, 'same');
        X_LR_Fea(:, :, 2) = conv2(X_LR_interp, vf1, 'same');
        X_LR_Fea(:, :, 3) = conv2(X_LR_interp,hf2,'same');
        X_LR_Fea(:, :, 4) = conv2(X_LR_interp,vf2,'same');


        % Extract random patches: HR and corresponding subsampled LR
        patches_LR_features = zeros(4*dim_h, length(xrow_l));
        patches_HRLR = zeros(dim_h, length(xrow_l));

        for ii = 1:num_patch
            row_h = xrow_h(ii); col_h = ycol_h(ii);
            row_l = xrow_l(ii); col_l = ycol_l(ii);

            % LR (features)
            Lpatch_Fea1 = X_LR_Fea(row_h:row_h+patchsize_h-1,col_h:col_h+patchsize_h-1,1);
            Lpatch_Fea2 = X_LR_Fea(row_h:row_h+patchsize_h-1,col_h:col_h+patchsize_h-1,2);
            Lpatch_Fea3 = X_LR_Fea(row_h:row_h+patchsize_h-1,col_h:col_h+patchsize_h-1,3);
            Lpatch_Fea4 = X_LR_Fea(row_h:row_h+patchsize_h-1,col_h:col_h+patchsize_h-1,4);
            patches_LR_features(:,ii) = [Lpatch_Fea1(:); Lpatch_Fea2(:); Lpatch_Fea3(:); Lpatch_Fea4(:)];

            % HR-LR (bicubic interp)
            Lpatch_HRLR = X_HR_ref(row_h:row_h+patchsize_h-1,col_h:col_h+patchsize_h-1)- X_LR_interp(row_h:row_h+patchsize_h-1,col_h:col_h+patchsize_h-1);
            patches_HRLR(:,ii) = Lpatch_HRLR(:);              
        end
        patches_LR_features_all(:,num_patch*(plane_id-1)+1:num_patch*plane_id) = patches_LR_features; 
        patches_HRLR_all(:,num_patch*(plane_id-1)+1:num_patch*plane_id)= patches_HRLR;
    end
end

PATCHES_FILENAME=strcat('/data/ISOTROPIC/dictionary_learning/downsampling/space_ratio_03/trainingpatches_couplefeatures_spaceratio',num2str(space_spacing,'%.1d'),'_timeratio',num2str(time_spacing,'%.1d'),'_patchsize',num2str(size_l,'%.2d'),'_numpatch',num2str(Nt*numel(LTHS_idt)*num_patch,'%.6d'),'.mat');
save(PATCHES_FILENAME,'patches_LR_features_all','patches_HRLR_all');
