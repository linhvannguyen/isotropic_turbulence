clear all; close all; clc;

%% INITIAL PARAMS
space_spacing=3;
time_spacing=4;

dim_l_pca=81;
size_l=8; 
patchsize_h = size_l*space_spacing; % size of HR patches
patchsize_l = patchsize_h; % 4x4 LR patches

num_patch=8*8; % number of patches 
dim_h=patchsize_h^2;
dim_l=4*patchsize_l^2;

params_train.lambda=0.2; % params.lambda=1*1e-5, overlap=7
params_train.K=2*(dim_h+dim_l_pca);
filename_SR=strcat('/data/ISOTROPIC/dictionary_learning/downsampling/space_ratio_03/SR_couplefeatures_K',num2str(params_train.K,'%.4d'),'_lambda',strrep(num2str(params_train.lambda,'%.2f'),'.',''),'.nc');
filename_ref='/data/ISOTROPIC/data/data_downsampled4.nc';
nc = netcdf(filename_ref,'r');
Nt = nc('Nt').itsDimsize;
Nx = nc('Nx').itsDimsize;
Ny = nc('Ny').itsDimsize;
Nz = nc('Nz').itsDimsize;
close(nc)


Nh=Nx;
Nl=Nh/space_spacing;
midplane_ids=3:time_spacing:Nh;


LTHS_idt=1:time_spacing:Nh;
[gridz_HR, gridy_HR] = meshgrid(1:Nh,1:Nh);
gridz_LR = gridz_HR(1:space_spacing:Nh,1:space_spacing:Nh);
gridy_LR = gridy_HR(1:space_spacing:Nh,1:space_spacing:Nh);


k_max=Nh/2;
k_1D_HR=[0:k_max -k_max+1:1:-1];
[KZS,KYS] = meshgrid(k_1D_HR,k_1D_HR); 
k_cutoff = k_max/space_spacing;

%% FFT
% fftn to frequency space
nc1=netcdf(filename_ref,'r');
x_ref_all = nc1{'velocity_x'}(:,:,:,midplane_ids);
close(nc1);

nc2=netcdf(filename_SR,'r');
x_SR_all = nc2{'Zhat_all'}(:,:,:,:);
close(nc2);

left=patchsize_h/2; right = patchsize_h/2;
bottom=patchsize_h/2; top = patchsize_h/2;

E_ref=zeros(Nh,1);
E_interp=zeros(Nh,1);
E_SR=zeros(Nh,1);

E_err_ref=zeros(Nh,1);
E_err_interp=zeros(Nh,1);
E_err_SR=zeros(Nh,1);

norm_factor = 1/(Nt*numel(midplane_ids));
for t=1:Nt
    for i=1:numel(midplane_ids)
        x_ref=squeeze(x_ref_all(t,:,:,i));
        x_SR=squeeze(x_SR_all(t,:,:,i));
        x_interp=resize_nguyen(resize_nguyen(x_ref, 1/space_spacing,'bicubic'), space_spacing,'bicubic');
        
        % Spect
        E = estimate_spect_2D(x_ref); 
        E_ref = E_ref + norm_factor*E;

        E = estimate_spect_2D(x_interp);
        E_interp = E_interp + norm_factor*E;
              
        E = estimate_spect_2D(x_SR);
        E_SR = E_SR + norm_factor*E;
              
        % Error
        E = estimate_spect_2D(x_ref - x_interp);
        E_err_interp = E_err_interp + norm_factor*E;
               
        E = estimate_spect_2D(x_ref - x_SR);
        E_err_SR = E_err_SR + norm_factor*E;  
    end
end


%%
fsize=20;
fname='CMU Serif';

%%
% -5/3 line
xline=[4,40];
yline=10^1*[exp((-5/3)*log(4)),exp((-5/3)*log(50))];

figure();

% Change default axes fonts.
set(0,'DefaultAxesFontName', fname)
set(0,'DefaultAxesFontSize', fsize)

% Change default text fonts.
set(0,'DefaultTextFontname', fname)
set(0,'DefaultTextFontSize', fsize)

set(gcf, 'Position', [400 100 800 600]);
set(gcf, 'Color', 'w');

loglog(2:k_max,E_ref(2:k_max),'k-*','LineWidth',1.5); 
hold on;
loglog(2:k_max,E_interp(2:k_max),'g-*','LineWidth',1.5); 
loglog(2:k_max,E_SR(2:k_max),'r-*','LineWidth',1.5); 
loglog(xline,yline,'k-','LineWidth',2)
plot([k_cutoff,k_cutoff] ,[1e-8,1e-2],'k--','LineWidth',2)
hold off
xlim([1,70]);
ylim([1e-8,5]);

text(20,0.5,'-5/3','HorizontalAlignment','right','FontSize',fsize+2)

xlabel('k'); ylabel('E(k)');

filename_fig_spect=strcat('./figures/spectra2d_spacespacing_',num2str(space_spacing,'%.2d'));
export_fig(filename_fig_spect,'-eps','-q101','-a4','-nocrop');
close()

%%
figure();

% Change default axes fonts.
set(0,'DefaultAxesFontName', fname)
set(0,'DefaultAxesFontSize', fsize)

% Change default text fonts.
set(0,'DefaultTextFontname', fname)
set(0,'DefaultTextFontSize', fsize)

set(gcf, 'Position', [400 100 800 600]);
set(gcf, 'Color', 'w');

h1=loglog(2:k_max,E_err_interp(2:k_max)./E_ref(2:k_max),'g-*','LineWidth',1.5); 
hold on;
h2=loglog(2:k_max,E_err_SR(2:k_max)./E_ref(2:k_max),'b-*','LineWidth',1.5); 
plot([k_cutoff,k_cutoff] ,[10^-8,10^0],'k--','LineWidth',2)
plot([1,100] ,[1,1],'k--','LineWidth',2)
hold off
xlim([0,60]);
ylim([1e-6,10]);

xlabel('k'); ylabel('E(k)');

leg=legend([h1 h2],{'$Interp (space)$','$SR$'},'location','northwest');
set(leg,'FontSize',fsize-4);
set(leg,'Interpreter','Latex')

legend boxoff

filename_fig_spect=strcat('./figures/errspectra2d_normalized_spacespacing_',num2str(space_spacing,'%.2d'));
export_fig(filename_fig_spect,'-eps','-q101','-a4','-nocrop');
close()


%%
figure();

% Change default axes fonts.
set(0,'DefaultAxesFontName', fname)
set(0,'DefaultAxesFontSize', fsize)

% Change default text fonts.
set(0,'DefaultTextFontname', fname)
set(0,'DefaultTextFontSize', fsize)

set(gcf, 'Position', [400 100 800 600]);
set(gcf, 'Color', 'w');


h1=loglog(2:k_max,E_err_interp(2:k_max),'g-*','LineWidth',1.5); 
hold on;
h2=loglog(2:k_max,E_err_SR(2:k_max),'b-*','LineWidth',1.5); 
plot([k_cutoff,k_cutoff] ,[10^-8,1.3*10^-3],'r--','LineWidth',2)
hold off
xlim([0,60]);
ylim([5*1e-7,1e-2]);

xlabel('k'); ylabel('E(k)');

leg=legend([h1 h2],{'$Interp (space)$','$SR$'},'location','northwest');
set(leg,'FontSize',fsize-4);
set(leg,'Interpreter','Latex')

legend boxoff

filename_fig_spect=strcat('./figures/errspectra2d_timespacing_',num2str(time_spacing,'%.2d'),'_spacespacing_',num2str(space_spacing,'%.2d'));
export_fig(filename_fig_spect,'-eps','-q101','-a4','-nocrop');
close()
