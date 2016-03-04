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

%% Random indices at HR or LR
xl = randperm(Nl-patchsize_l+1); % starting point of LR patches
yl = randperm(Nl-patchsize_l+1); % total possible number of patches: xl*yl
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

%% EXTRACT PATCHES

patches_HR_all=zeros(dim_h,num_patch*Nt*numel(LTHS_idt));
patches_LR_all = zeros(dim_l,num_patch*Nt*numel(LTHS_idt));

index=1;
for t=1:Nt
    t
    for i=1:numel(LTHS_idt)
        X_HR=squeeze(u_HR_all(t,:,:,i));       
        X_LR = X_HR(1:space_spacing:Nh,1:space_spacing:Nh);
        fig1=figure(); imagesc(X_HR); caxis([-3,3]);
        fig2=figure(); imagesc(X_LR); caxis([-3,3]);
        % Extract random patches: HR and corresponding subsampled LR
        patches_HR = zeros(dim_h, length(xrow_h));
        patches_LR = zeros(dim_l, length(xrow_l));
        for ii = 1:num_patch
            row_h = xrow_h(ii); col_h = ycol_h(ii);
            row_l = xrow_l(ii); col_l = ycol_l(ii);

            % HR
            Hpatch = X_HR(row_h:row_h+patchsize_h-1,col_h:col_h+patchsize_h-1);
            patches_HR(:,ii)=Hpatch(:);

            % LR
            Lpatch = X_LR(row_l:row_l+patchsize_l-1,col_l:col_l+patchsize_l-1);
            patches_LR(:,ii) = Lpatch(:);           
        end
        
        patches_HR_all(:,index:index+num_patch-1)= patches_HR;
        patches_LR_all(:,index:index+num_patch-1) = patches_LR;

        fig3=figure(); imagesc(reshape(patches_HR_all(:,index),patchsize_h,patchsize_h)); caxis([-3,3]);
        fig4=figure(); imagesc(reshape(patches_LR_all(:,index),patchsize_l,patchsize_l)); caxis([-3,3]);
        
        index=index+num_patch;
    end
end
PATCHES_FILENAME=strcat('/data/ISOTROPIC/dictionary_learning/subsampling/space_ratio_03/trainingpatches_coupleHRLR_spaceratio',num2str(space_spacing,'%.1d'),'_timeratio',num2str(time_spacing,'%.1d'),'_patchsize',num2str(patchsize_l,'%.2d'),'_numpatch',num2str(Nt*numel(LTHS_idt)*num_patch,'%.6d'),'.mat');
save(PATCHES_FILENAME,'patches_HR_all','patches_LR_all');
