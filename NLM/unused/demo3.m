% clear all; close all; clc;

%% INITIAL PARAMS
space_spacing=3;
time_spacing=4;

simpatch_haftsize = 6; % the half width of the patche
simpatch_fullsize = 2*simpatch_haftsize+1; % the full width.
accpatch_haftsize = 6; 
accpatch_fullsize = 2*accpatch_haftsize + 1; 
neighbor_haftsize = 5; % the half width of the patch
neighbor_fullsize = 2*neighbor_haftsize+1; % the full width.
dim_simpatch=simpatch_fullsize*simpatch_fullsize;
dim_accpatch=accpatch_fullsize*accpatch_fullsize;
dim_neighbor=neighbor_fullsize*neighbor_fullsize;

acc_ids_1D = simpatch_haftsize+1-accpatch_haftsize:simpatch_haftsize+1+accpatch_haftsize;
acc_ids = reshape(1:dim_simpatch,[simpatch_fullsize,simpatch_fullsize]);
acc_ids = acc_ids(acc_ids_1D, acc_ids_1D);
acc_ids = acc_ids(:);

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
tau = .35;
%% Extract patches
% grids of reference pixels
[gridz,gridy] = meshgrid(1:Nz,1:Ny);

% neighboring offsets
[dz,dy] = meshgrid(-neighbor_haftsize:neighbor_haftsize,-neighbor_haftsize:neighbor_haftsize);
dz = reshape(dz(:), [1 1 dim_neighbor]);
dy = reshape(dy(:), [1 1 dim_neighbor]);
gridz = repmat(gridz, [1 1 dim_neighbor]) + repmat(dz, [Ny Nz 1]);
gridy = repmat(gridy, [1 1 dim_neighbor]) + repmat(dy, [Ny Nz 1]);

% similarity patch offsets
[dz,dy] = meshgrid(-simpatch_haftsize:simpatch_haftsize,-simpatch_haftsize:simpatch_haftsize);
dz = reshape(dz(:), [1 1 1 dim_simpatch]);
dy = reshape(dy(:), [1 1 1 dim_simpatch]);
gridz = repmat(gridz, [1 1 1 dim_simpatch]) + repmat(dz, [Ny Nz dim_neighbor 1]);
gridy = repmat(gridy, [1 1 1 dim_simpatch]) + repmat(dy, [Ny Nz dim_neighbor 1]);

% handle boundary condition by periodicity
gridy = border_per(gridy,1,Ny); 
gridz = border_per(gridz,1,Nz); 

start=tic();
% Reference patches of size (n,n,w1,w1). 
ref_patches = extract_patch(u_noisy, squeeze(gridy(:,:,ceil(dim_neighbor/2),:,:)), squeeze(gridz(:,:,ceil(dim_neighbor/2),:,:)));
ref_patches = reshape(ref_patches,[Ny Nz 1 dim_simpatch]);
ref_patches = repmat(ref_patches, [1 1 dim_neighbor  1]);

% Moving patches of size (n,n,m,m,w1,w1). 
moving_patches = extract_patch(u_noisy, gridy, gridz); 
fprintf(['\n Load data in ',num2str(toc(start),'%.3f'),' seconds \n'])

%% PCA
% start=tic();
% resh = @(P)reshape(P, [dim_im*dim_neighbor dim_simpatch])'; 
% 
% remove_mean = @(Q)Q - repmat(mean(Q), [dim_simpatch 1]);
% moving_patches = remove_mean(resh(moving_patches));
% C = moving_patches*moving_patches';
% [V,ener] = eig(C); ener = diag(ener);
% % [ener,~] = sort(ener, 'descend'); 
% ener = diag(ener);
% ener = cumsum(ener) / sum(ener);
% 
% k = find(ener >= 1e-2, 1); % ignore 1% energy
% dim_pca=numel(ener)-k;
% V_pca = V(:, end-dim_pca+1:end); 
% ref_patches = V(:,1:dim_pca)' * remove_mean(resh(ref_patches));
% moving_patches= V(:,1:dim_pca)' * remove_mean(resh(moving_patches));
% 
% resh_i = @(P)reshape(P, [Ny Nz dim_neighbor dim_pca]); 
% ref_patches=resh_i(ref_patches');
% moving_patches=resh_i(moving_patches');
% fprintf(['\n PCA in ',num2str(toc(start),'%.3f'),' seconds \n'])

%% Compute distance
D = (ref_patches - moving_patches).^2;
D = mean(D,4); 

K=exp(-D/(2.0*tau*tau));

%%
start=tic();

x_2D=zeros(Ny,Nz); 
W_2D=zeros(Ny,Nz);

for i=1:Ny
    for j=1:Nz
        patches=extract_patch(u_noisy, squeeze(gridy(i,j,:,acc_ids)), squeeze(gridz(i,j,:,acc_ids)));
        w = squeeze(K(i,j,:));
        rows_ref = squeeze(gridy(i,j,ceil(dim_neighbor/2),acc_ids));
        cols_ref = squeeze(gridz(i,j,ceil(dim_neighbor/2),acc_ids));
        [x_2D, W_2D]=put_patch(x_2D, W_2D, patches, rows_ref, cols_ref,w);
    end
end 
x_2D=x_2D./W_2D;
fprintf(['\n Reconstruction in ',num2str(toc(start),'%.3f'),' seconds \n'])

%% Display some example of patches
y=zeros(16,1);
z=zeros(16,1);
for i=1:16
    y(i) = floor(rand*(Ny-1)+1);
    z(i) = floor(rand*(Nz-1)+1);    
end

fig1=figure();
for i=1:16 
    subplot(4,4,i);
    imagesc(reshape(ref_patches(y(i),z(i),1,:),[simpatch_fullsize simpatch_fullsize]));
    axis off;
end

fig2=figure();
for i=1:16
    subplot(4,4,i);
    imagesc(reshape(moving_patches(y(i),z(i),ceil(dim_neighbor/2),:),[simpatch_fullsize simpatch_fullsize]));
    axis off;
end

fig3=figure();
for i=1:16
    subplot(4,4,i);
    imagesc(reshape(K(y(i),z(i),:),[neighbor_fullsize neighbor_fullsize]));
    axis off; axis equal;
end


fig4=figure();
subplot(1,3,1);imagesc(u_noisy); caxis([-3,3]); axis equal;
subplot(1,3,2);imagesc(x_2D); caxis([-3,3]); axis equal;
subplot(1,3,3);imagesc(u_org); caxis([-3,3]); axis equal;