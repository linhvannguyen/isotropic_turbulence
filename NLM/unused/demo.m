clear all; close all; clc;

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


%% Load data
u_org=squeeze(u_HR_all(1,:,:,1));
sigma = .4;
u_noisy= u_org + randn(Nh,Nh)*sigma;
tau = .3;
sample_point_ids = [40 40]; % sample point to compute distance and kernel

%% We set up large (n,n,w1,w1) matrices to index the the X and Y position of the pixel to extract
patch_haftwidth = 3; % the half width of the patch
patch_width = 2*patch_haftwidth+1; % the full width.

% location of pixels
[Y,X] = meshgrid(1:Nh,1:Nh);

% offsets
[dY,dX] = meshgrid(-patch_haftwidth:patch_haftwidth,-patch_haftwidth:patch_haftwidth);

% location of pixels to extract
dX = reshape(dX, [1 1 patch_width patch_width]);
dY = reshape(dY, [1 1 patch_width patch_width]);
X = repmat(X, [1 1 patch_width patch_width]) + repmat(dX, [Nh Nh 1 1]);
Y = repmat(Y, [1 1 patch_width patch_width]) + repmat(dY, [Nh Nh 1 1]);
 
% We handle boundary condition by reflexion
X(X<1) = 2-X(X<1); 
Y(Y<1) = 2-Y(Y<1); 
X(X>Nh) = 2*Nh-X(X>Nh); 
Y(Y>Nh) = 2*Nh-Y(Y>Nh);

%% Extract patches
% Patch extractor operator
patch = @(f)f(X + (Y-1)*Nh);

% Define the patch matrix P of size (n,n,w1,w1). 
% Each P(i,j,:,:) represent an (w1,w1) patch extracted around pixel (i,j) in the image. 
% u (ids) where ids is 4D, u (MxN) is 2D, will take points in 4D matrix, each has id from 1 to MxN (row first, column after)  
P = patch(u_noisy);

% Display some example of patches
clf;
for i=1:16
    subplot(4,4,i);
    x = floor( rand*(Nh-1)+1 );
    y = floor( rand*(Nh-1)+1 );
    imagesc( squeeze(P(x,y,:,:)));
end
 
%% Dimensionality Reduction with PCA

% Target dimensionality d
dim_pca = 30; 

% Turn the patch matrix into an (w1*w1,n*n) array, so that each P(:,i) is a w1*w1 vector representing a patch
resh = @(P)reshape(P, [Nh*Nh patch_width*patch_width])'; 

% Remove mean
remove_mean = @(Q)Q - repmat(mean(Q), [patch_width*patch_width 1]);

% Compute the mean and the covariance of the points cloud representing the patches
P1 = remove_mean(resh(P));
C = P1*P1';

% Extract the eigenvectors, sorted by decreasing amplitude
[V,D] = eig(C); D = diag(D);
[D,I] = sort(D, 'descend'); V = V(:,I);

% Display the decaying amplitude of the eigenvalues.
clf;
plot(D); axis('tight');

% Display the leading eigenvectors - they look like Fourier modes
clf;
for i=1:16
    subplot(4,4,i);
    imagesc( reshape(V(:,i),[patch_width patch_width]));
end

% Patch dimensionality reduction operator
iresh = @(Q)reshape(Q', [Nh Nh dim_pca]);
descriptor = @(f)iresh( V(:,1:dim_pca)' * remove_mean(resh(P)) );

% Each H(i,j,:) is a d-dimensional descriptor of a patch
H = descriptor(u_noisy);
 
%% Localizing the Non-local Means

% We set a "locality constant" q that set the maximum distance between patches to compare. 
% This allows to speed up computation, and makes NL-means type methods semi-global (to avoid searching in all the image)
q = 3;

% selection use clamp func to handle boundary conditions,
% clamp(1-3:1+3,1,Ny)=[1,1,1,1,2,3,4]; clamp(Ny-3:Ny+3,1,Ny)=[Ny-3,Ny-2,Ny-2,Ny,Ny,Ny,Ny]
selection = @(i){clamp(i(1)-q:i(1)+q, 1,Nh), clamp(i(2)-q:i(2)+q,1,Nh)};

% Compute distance and kernel only within the window
distance = @(i,sel)sum( (H(sel{1},sel{2},:) - repmat(H(i(1),i(2),:), ...
        [length(sel{1}) length(sel{2}) 1])).^2, 3 )/(patch_width*patch_width);
distance = @(i)distance(i,selection(i));
normalize = @(K)K/sum(K(:));
kernel = @(i,tau)normalize( exp( -distance(i)/(2*tau^2) ) );

% Compute a typical example of kernel for some pixel position (x,y)
D = distance(sample_point_ids);
K = kernel(sample_point_ids,tau);

% Display the squared distance and the kernel
clf;
subplot(1,2,1); imagesc(D);
subplot(1,2,2); imageplot(K);

% The NL-filtered value at pixel (x,y) is obtained by averaging the values of f with the weight K
NLval = @(K,sel)sum(sum(K.*u_noisy(sel{1},sel{2})));
NLval = @(i,tau)NLval( kernel(i,tau), selection(i) );

% We apply the filter to each pixel location to perform the NL-means algorithm
[Y,X] = meshgrid(1:Nh,1:Nh);
NLmeans = @(tau)arrayfun(@(i1,i2)NLval([i1 i2],tau), X,Y);

u_denoised=NLmeans(tau);

%%
h1=figure();
pcolor (gridz_HR,gridy_HR,u_noisy); caxis([-3,3]);
shading flat;

h2=figure();
pcolor (gridz_HR,gridy_HR,u_denoised); caxis([-3,3]);
shading flat;

h3=figure(); 
pcolor (gridz_HR,gridy_HR,u_org); caxis([-3,3]);
shading flat;