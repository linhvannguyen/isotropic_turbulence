clear all; close all; clc;
addpath('./funcs/');

%% INITIAL PARAMETERS
space_spacing=03;
time_spacing=04;

filename_ref='/data/ISOTROPIC/data/data_downsampled4.nc';
filename_interp_space=strcat('/data/ISOTROPIC/fusion/Uinterp_spatialspacing_',num2str(space_spacing,'%.2d'),'.nc');
filename_interp_time=strcat('/data/ISOTROPIC/fusion/Uinterp_timespacing_',num2str(time_spacing,'%.2d'),'.nc');
filename_fusion=strcat('/data/ISOTROPIC/fusion/FusedData_diagCov_timespacing_',num2str(time_spacing,'%.2d'),'_spacespacing_',num2str(space_spacing,'%.2d'),'.nc');

nc = netcdf(filename_ref,'r');
Nt = nc('Nt').itsDimsize;
Nx = nc('Nx').itsDimsize;
Ny = nc('Ny').itsDimsize;
Nz = nc('Nz').itsDimsize;
close(nc);

N_HR=Ny;
k_max=N_HR/2;
k_1D_HR=[0:k_max -k_max+1:1:-1];
[KZS,KYS] = meshgrid(k_1D_HR,k_1D_HR); 

k_cutoff = k_max/space_spacing;

%% 2D, ratio:4

% fftn to frequency space
nc1=netcdf(filename_ref,'r');
nc2=netcdf(filename_interp_space,'r');
nc3=netcdf(filename_interp_time,'r');
nc4=netcdf(filename_fusion,'r');

E_ref=zeros(N_HR,1);
E_interp_space=zeros(N_HR,1);
E_interp_time=zeros(N_HR,1);
E_fusion=zeros(N_HR,1);

E_err_ref=zeros(N_HR,1);
E_err_interp_space=zeros(N_HR,1);
E_err_interp_time=zeros(N_HR,1);
E_err_fusion=zeros(N_HR,1);

pos_t=3;
for t=1:Nt
    for i=pos_t:time_spacing:Nx
        x_ref=nc1{'velocity_x'}(t,:,:,i);
        x_interp_space=nc2{'Uinterp'}(t,:,:,i);
        x_interp_time=nc3{'Uinterp'}(t,:,:,i);
        x_fusion=nc4{'Zhat_all'}(t,:,:,i);
        
        % Spect
        E = estimate_spect_2D(x_ref);
        E_ref = E_ref + 1/(Nt*numel(pos_t:time_spacing:Nx))*E;

        E = estimate_spect_2D(x_interp_space);
        E_interp_space = E_interp_space + 1/(Nt*numel(pos_t:time_spacing:Nx))*E;
        
        E = estimate_spect_2D(x_interp_time);
        E_interp_time = E_interp_time + 1/(Nt*numel(pos_t:time_spacing:Nx))*E;
        
        E = estimate_spect_2D(x_fusion);
        E_fusion = E_fusion + 1/(Nt*numel(pos_t:time_spacing:Nx))*E;
        
        
        % Error
        E = estimate_spect_2D(x_ref - x_interp_space);
        E_err_interp_space = E_err_interp_space + 1/(Nt*numel(pos_t:time_spacing:Nx))*E;
        
        E = estimate_spect_2D(x_ref - x_interp_time);
        E_err_interp_time = E_err_interp_time + 1/(Nt*numel(pos_t:time_spacing:Nx))*E;
        
        E = estimate_spect_2D(x_ref - x_fusion);
        E_err_fusion = E_err_fusion + 1/(Nt*numel(pos_t:time_spacing:Nx))*E;
        
    end
end

close(nc1); close(nc2); close(nc3); close(nc4);



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

h1=loglog(2:k_max,E_ref(2:k_max),'k-*','LineWidth',1.5); 
hold on;
h2=loglog(2:k_max,E_interp_space(2:k_max),'g-*','LineWidth',1.5); 
h3=loglog(2:k_max,E_interp_time(2:k_max),'m-*','LineWidth',1.5); 
h4=loglog(2:k_max,E_fusion(2:k_max),'r-*','LineWidth',1.5); 
loglog(xline,yline,'k-','LineWidth',2)
hold off
xlim([1,70]);
ylim([1e-8,5]);

text(20,0.5,'-5/3','HorizontalAlignment','right','FontSize',fsize+2)

xlabel('k'); ylabel('E(k)');

leg=legend([h1 h2 h3 h4],{'Reference','$\mathbf{I}_t \mathbf{x}$','$\mathbf{I}_s \mathbf{y}$','$Fusion$'},'interpreter', 'latex','location','southwest');
set(leg,'FontSize',fsize-4);
legend boxoff


filename_fig_spect=strcat('./figures/spectra2d_timespacing_',num2str(time_spacing,'%.2d'),'_spacespacing_',num2str(space_spacing,'%.2d'));
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

h1=loglog(2:k_max,E_err_interp_space(2:k_max)./E_ref(2:k_max),'g-*','LineWidth',1.5); 
hold on;
h2=loglog(2:k_max,E_err_interp_time(2:k_max)./E_ref(2:k_max),'m-*','LineWidth',1.5); 
h3=loglog(2:k_max,E_err_fusion(2:k_max)./E_ref(2:k_max),'r-*','LineWidth',1.5); 
plot([k_cutoff,k_cutoff] ,[10^-6,10^0],'r--','LineWidth',2)
plot([1,100] ,[1,1],'k--','LineWidth',2)
hold off
xlim([0,60]);
ylim([1e-6,10]);

xlabel('k'); ylabel('E(k)');

leg=legend([h1 h2 h3],{'$Interp (space)$','$Interp (time)$','$Fusion$'},'location','northwest');
set(leg,'FontSize',fsize-4);
set(leg,'Interpreter','Latex')

legend boxoff


filename_fig_spect=strcat('./figures/errspectra2d_normalized_timespacing_',num2str(time_spacing,'%.2d'),'_spacespacing_',num2str(space_spacing,'%.2d'));
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


h1=loglog(2:k_max,E_err_interp_space(2:k_max),'g-*','LineWidth',1.5); 
hold on;
h2=loglog(2:k_max,E_err_interp_time(2:k_max),'m-*','LineWidth',1.5); 
h3=loglog(2:k_max,E_err_fusion(2:k_max),'r-*','LineWidth',1.5); 
plot([k_cutoff,k_cutoff] ,[10^-6,1.5*10^-3],'r--','LineWidth',2)
hold off
xlim([0,60]);
ylim([5*1e-7,1e-2]);

xlabel('k'); ylabel('E(k)');

leg=legend([h1 h2 h3],{'$Interp (space)$','$Interp (time)$','$Fusion$'},'location','northwest');
set(leg,'FontSize',fsize-4);
set(leg,'Interpreter','Latex')

legend boxoff

filename_fig_spect=strcat('./figures/errspectra2d_timespacing_',num2str(time_spacing,'%.2d'),'_spacespacing_',num2str(space_spacing,'%.2d'));
export_fig(filename_fig_spect,'-eps','-q101','-a4','-nocrop');
close()
