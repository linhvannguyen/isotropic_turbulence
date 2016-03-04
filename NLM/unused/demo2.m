% clear all; close all; clc;

%% INITIAL PARAMS
space_spacing=3;
time_spacing=4;

simpatch_haftsize = 6; % the half width of the patche
simpatch_fullsize = 2*simpatch_haftsize+1; % the full width.
accpatch_haftsize = 0; 
accpatch_fullsize = 2*accpatch_haftsize + 1; 
neighbor_haftsize = 5; % the half width of the patch
neighbor_fullsize = 2*neighbor_haftsize+1; % the full width.
dim_simpatch=simpatch_fullsize*simpatch_fullsize;
dim_accpatch=accpatch_fullsize*accpatch_fullsize;
dim_neighbor=neighbor_fullsize*neighbor_fullsize;

filename_ref='/data/ISOTROPIC/data/data_downsampled4.nc';
nc = netcdf(filename_ref,'r');
Nt = nc('Nt').itsDimsize;
Nz = nc('Nz').itsDimsize;
Ny = nc('Ny').itsDimsize;
Nx = nc('Nx').itsDimsize;
u_HR_all=nc{'velocity_x'}(1:Nt,1:Nz,1:Ny,1:space_spacing:Nx);
close(nc)
dim_im=Ny*Nz;

%% Load data
u_org=squeeze(u_HR_all(1,:,:,1));
sigma = .4;
u_noisy= u_org + randn(Ny,Nz)*sigma;
tau = .4;
% extract_patch = @(f, rows, cols)f(rows + (cols-1)*size(rows,1));

%% Extract patches
% grids of reference pixels
[gridz,gridy] = meshgrid(1:Nz,1:Ny);

% neighboring offsets
[dz,dy] = meshgrid(-neighbor_haftsize:neighbor_haftsize,-neighbor_haftsize:neighbor_haftsize);
dz = reshape(dz, [1 1 neighbor_fullsize neighbor_fullsize]);
dy = reshape(dy, [1 1 neighbor_fullsize neighbor_fullsize]);
gridz = repmat(gridz, [1 1 neighbor_fullsize neighbor_fullsize]) + repmat(dz, [Ny Nz 1 1]);
gridy = repmat(gridy, [1 1 neighbor_fullsize neighbor_fullsize]) + repmat(dy, [Ny Nz 1 1]);

% similarity patch offsets
[dz,dy] = meshgrid(-simpatch_haftsize:simpatch_haftsize,-simpatch_haftsize:simpatch_haftsize);
dz = reshape(dz, [1 1 1 1 simpatch_fullsize simpatch_fullsize]);
dy = reshape(dy, [1 1 1 1 simpatch_fullsize simpatch_fullsize]);
gridz = repmat(gridz, [1 1 1 1 simpatch_fullsize simpatch_fullsize]) + repmat(dz, [Ny Nz neighbor_fullsize neighbor_fullsize 1 1]);
gridy = repmat(gridy, [1 1 1 1 simpatch_fullsize simpatch_fullsize]) + repmat(dy, [Ny Nz neighbor_fullsize neighbor_fullsize 1 1]);

% handle boundary condition by periodicity
gridy = border_per(gridy,1,Ny); 
gridz = border_per(gridz,1,Nz); 

start=tic();
% Reference patches
% Define the patch matrix P of size (n,n,w1,w1). 
% Each P(i,j,:,:) represent an (w1,w1) patch extracted around pixel (i,j) in the image. 
% u (ids) where ids is 4D, u (MxN) is 2D, will take points in 4D matrix, each has id from 1 to MxN (row first, column after)  
ref_patches = extract_patch(u_noisy, squeeze(gridy(:,:,neighbor_haftsize+1,neighbor_haftsize+1,:,:)), ...
    squeeze(gridz(:,:,neighbor_haftsize+1,neighbor_haftsize+1,:,:)));
ref_patches = reshape(ref_patches,[Ny Nz 1 1 simpatch_fullsize simpatch_fullsize]);
ref_patches = repmat(ref_patches, [1 1 neighbor_fullsize neighbor_fullsize  1 1]);

% Moving patches
% Define the patch matrix P of size (n,n,m,m,w1,w1). 
% Each P(i,j,k,l:,:) represent an (w1,w1) patch extracted around pixel (k,l), which is
% in the searching region around pixel (i,j) in the image. 
moving_patches = extract_patch(u_noisy, gridy, gridz); 
fprintf(['\n Load data in ',num2str(toc(start),'%.3f'),' seconds \n'])

%% Display some example of patches
% y=zeros(16); z=y;
% for i=1:16
%     y(i) = floor(rand*(Ny-1)+1 );
%     z(i) = floor(rand*(Nz-1)+1 );    
% end
% 
% fig1=figure();
% for i=1:16 
%     subplot(4,4,i);
%     imagesc( squeeze(ref_patches(y(i),z(i),1,1,:,:)));
%     axis off;
% end
% 
% fig2=figure();
% for i=1:16
%     subplot(4,4,i);
%     imagesc( squeeze(moving_patches(y(i),z(i),neighbor_haftsize+1,neighbor_haftsize+1,:,:)));
%     axis off;
% end


%% Compute distance
D = (ref_patches - moving_patches).^2;
D = mean(D,6); 
D = mean(D,5);

K=exp(-D/(2.0*tau*tau));
clearvars ref_patches moving_patches;
%%
% fig3=figure();
% for i=1:16
%     subplot(4,4,i);
%     imagesc( squeeze(K(y(i),z(i),:,:)));
%     axis off; axis equal;
% end

%%
x_2D=zeros(Ny,Nz); 
W_2D=zeros(Ny,Nz);

acc_ids = simpatch_haftsize+1-accpatch_haftsize:simpatch_haftsize+1+accpatch_haftsize;
for i=1:Ny
    for j=1:Nz
        patches=extract_patch(u_noisy, squeeze(gridy(i,j,:,:,acc_ids,acc_ids)), squeeze(gridz(i,j,:,:,acc_ids,acc_ids)));
        w = squeeze(K(i,j,:,:));
        rows_ref = squeeze(gridy(i,j,neighbor_haftsize+1,neighbor_haftsize+1,acc_ids,acc_ids));
        cols_ref = squeeze(gridz(i,j,neighbor_haftsize+1,neighbor_haftsize+1,acc_ids,acc_ids));
        [x_2D, W_2D]=put_patch(x_2D, W_2D, patches, rows_ref, cols_ref,w);
    end
end 
x_2D=x_2D./W_2D;
toc()
%%
fig4=figure();
subplot(1,3,1);imagesc(u_noisy); axis equal;
subplot(1,3,2);imagesc(x_2D); axis equal;
subplot(1,3,3);imagesc(u_org); axis equal;