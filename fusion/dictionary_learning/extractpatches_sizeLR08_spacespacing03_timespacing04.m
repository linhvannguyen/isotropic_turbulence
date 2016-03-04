clear all; close all; clc;

%% INITIAL PARAMS
space_spacing=03; % subsampling ration in space
time_spacing=04; % subsampling ration in time (from 40Hz to 4Hz)

filename_ref='/data/ISOTROPIC/data/data_downsampled4.nc';

nc = netcdf(filename_ref,'r');
Nt = nc('Nt').itsDimsize;
Nh = nc('Nx').itsDimsize;
u_HR_all=nc{'velocity_x'}(1:Nt,1:Nh,1:Nh,1:time_spacing:Nh);
close(nc)

%% Define parameters
HTHS_idy=1:Nh; % row indices of HTHS in space 
HTHS_idz=1:Nh; % column indices of HTHS in space
HTLS_idt=1:1:Nh; % indices of all HTHS snapshots in time

HTLS_idy= 1:space_spacing:Nh; % row indices of HTLS in space
HTLS_idz= 1:space_spacing:Nh; % column indices of HTLS in space

LTHS_idt=1:time_spacing:Nh;


Nl=Nh/space_spacing;

% size_l=8; 
% patchsize_h = size_l*space_spacing; % size of HR patches
% num_patch=8*8; % number of patches 

size_l=4; 
patchsize_h = size_l*space_spacing; % size of HR patches
num_patch=16*16; % number of patches 

% size_l=2; 
% patchsize_h = size_l*space_spacing; % size of HR patches
% num_patch=16*16; % number of patches 

dim_h=patchsize_h^2;

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


%% EXTRACT PATCHES
num_planes=Nt*numel(LTHS_idt);
patches_smallscales_all = zeros(dim_h,num_patch*num_planes);

index=1;
for t=1:Nt
    t
    for i=1:numel(LTHS_idt)
        plane_id=LTHS_idt(i);
        x_HR_ref=squeeze(u_HR_all(t,:,:,i));    
        x_HR_smallscales = x_HR_ref - filter_2D(x_HR_ref, space_spacing, space_spacing);
        
        % Extract random patches: HR and corresponding subsampled LR
        patches_smallscales = zeros(dim_h, length(xrow_l));

        for ii = 1:num_patch
            row_h = xrow_h(ii); col_h = ycol_h(ii);

            patch_smallscales = x_HR_smallscales(row_h:row_h+patchsize_h-1,col_h:col_h+patchsize_h-1);
            patches_smallscales(:,ii) = patch_smallscales(:);              
        end
        patches_smallscales_all(:,index:index+num_patch-1)= patches_smallscales;
        index=index+num_patch;
    end
end

filename_patches=strcat('/data/ISOTROPIC/fusion/patches_smallscales_sizeLR',num2str(size_l,'%.2d'),'_timespacing',num2str(time_spacing,'%.2d'),'_spacespacing',num2str(space_spacing,'%.2d'),'.mat');

save(filename_patches,'patches_smallscales_all');
