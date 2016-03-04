clear all; close all; clc;

%% INITIAL PARAMS
space_spacing=3;
time_spacing=4;

simpatch_haftsize = 3; % the half width of the patche
simpatch_fullsize = 2*simpatch_haftsize+1; % the full width.
accpatch_haftsize = 1; 
accpatch_fullsize = 2*accpatch_haftsize + 1; 
neighbor_haftsize = 4; % the half width of the patch
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
close(nc)

dim_im=Ny*Nz;

sigma = .4;
tau = .35;

LTHS_idt=1:time_spacing:Nx;
HTHS_idy=1:Ny; % row indices of HTHS in space 
HTHS_idz=1:Nz; % column indices of HTHS in space
HTLS_idt=1:1:Nx; % indices of all HTHS snapshots in time
HTLS_idy= 1:space_spacing:Ny; % row indices of HTLS in space
HTLS_idz= 1:space_spacing:Nz; % column indices of HTLS in space

%% Load data
nc = netcdf(filename_ref,'r');
x_HR_all=nc{'velocity_x'}(1:Nt,1:Nz,1:Ny,1:space_spacing:Nx);
close(nc)

x_org=squeeze(x_HR_all(1,:,:,1));
x_noisy= x_org + randn(Ny,Nz)*sigma;

%%
[gridy, gridz] = grid_gen(Ny,Nz,simpatch_haftsize,neighbor_haftsize);

start=tic(); 
[x_rec,weights] = NLmean(x_noisy, x_noisy, x_noisy, gridy, gridz, simpatch_haftsize,accpatch_haftsize, neighbor_haftsize, tau);
x_rec=x_rec./weights;
fprintf(['\n Reconstruction in ',num2str(toc(start),'%.3f'),' seconds \n'])

%% fig4=figure();
subplot(2,2,1);imagesc(x_noisy); caxis([-3,3]); axis off; axis equal; axis tight;
subplot(2,2,2);imagesc(x_rec); caxis([-3,3]); axis off; axis equal; axis tight;
subplot(2,2,3);imagesc(x_org); caxis([-3,3]); axis off; axis equal; axis tight;
subplot(2,2,4);imagesc(x_org-x_rec); caxis([-3,3]); axis off; axis equal; axis tight;