clear all; close all; clc;

%% INITIAL PARAMETERS
space_spacing=4;
time_spacing=4;

filename_ref='/data/ISOTROPIC/data/data_downsampled4.nc';
filename_interp_space=strcat('/data/ISOTROPIC/fusion/Uinterp_spatialspacing_',num2str(space_spacing,'%.2d'),'.nc');
filename_interp_time=strcat('/data/ISOTROPIC/fusion/Uinterp_timespacing_',num2str(time_spacing,'%.2d'),'.nc');
filename_fusion=strcat('/data/ISOTROPIC/fusion/FusedData_diagCov_timespacing_',num2str(time_spacing,'%.2d'),'_spacespacing_',num2str(space_spacing,'%.2d'),'.nc');
filename_NLM=strcat('/data/ISOTROPIC/NLM/interpdiff/NLmean_propag2dirs_sspacing',num2str(space_spacing,'%.1d'),'_tspacing',num2str(time_spacing,'%.1d'),'_sim12_acc12_neighbor2_tau0100.nc');

nc = netcdf(filename_ref,'r');
Nt = nc('Nt').itsDimsize;
Nx = nc('Nx').itsDimsize;
Ny = nc('Ny').itsDimsize;
Nz = nc('Nz').itsDimsize;
close(nc);

N_HR=Ny;
LTHS_idt=1:time_spacing:Nx;

NRMSE= @(org,rec) sqrt(sum((org(:)-rec(:)).^2))/sqrt(sum((org(:)).^2));

%% 2D, ratio:4
nc1=netcdf(filename_ref,'r');
nc2=netcdf(filename_interp_space,'r');
nc3=netcdf(filename_interp_time,'r');
nc4=netcdf(filename_fusion,'r');
nc5=netcdf(filename_NLM,'r');

NRMSE_interp_space=zeros(time_spacing+1,1);
NRMSE_interp_time=zeros(time_spacing+1,1);
NRMSE_fusion=zeros(time_spacing+1,1);
NRMSE_NLM=zeros(time_spacing+1,1);

for pos_t=1:time_spacing
    x_ref=nc1{'velocity_x'}(:,:,:,pos_t:time_spacing:LTHS_idt(end));
    
    x_interp_space=nc2{'Uinterp'}(:,:,:,pos_t:time_spacing:LTHS_idt(end));
    NRMSE_interp_space(pos_t) = NRMSE(x_ref(:),x_interp_space(:));
    
    x_interp_time=nc3{'Uinterp'}(:,:,:,pos_t:time_spacing:LTHS_idt(end));
    NRMSE_interp_time(pos_t) = NRMSE(x_ref(:),x_interp_time(:)); clearvars x_interp_time;

    x_fusion=nc4{'Zhat_all'}(:,:,:,pos_t:time_spacing:LTHS_idt(end));
    NRMSE_fusion(pos_t) = NRMSE(x_ref(:),x_fusion(:)); clearvars x_fusion;    

    x_NLM = x_interp_space+ permute(nc5{'x_rec_all'}(:,pos_t:time_spacing:LTHS_idt(end),:,:), [1 3 4 2]);
    NRMSE_NLM(pos_t) = NRMSE(x_ref(:),x_NLM(:)); clearvars x_NLM; clearvars x_interp_space;
     
    clearvars x_ref;
end
NRMSE_interp_space(:) = sum(NRMSE_interp_space)/time_spacing;
close(nc1); close(nc2); close(nc3); close(nc5);
close(nc4);

%%
fsize=20;
fname='CMU Serif';

h=figure;
set(h, 'Position', [200 200 800 500]);
set(h,'color','w')

% Change default axes fonts.
set(0,'DefaultAxesFontName', fname)
set(0,'DefaultAxesFontSize', fsize)

% Change default text fonts.
set(0,'DefaultTextFontname', fname)
set(0,'DefaultTextFontSize', fsize)

% define how the figure inside the plot appear on the paper
set(gcf,'Units','normal');
set(gca,'Position',[0.15 0.2 0.725 0.75]); % [x_leftlowcorner y_leftlowcorner width height]
set(gcf,'Units','pixels');


hold on;
plot(1:time_spacing+1,NRMSE_interp_space,'g-s','LineWidth',1,'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',6)
plot(1:time_spacing+1,NRMSE_interp_time,'m-o','LineWidth',1,'MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',6)
plot(1:time_spacing+1,NRMSE_fusion,'r-d','LineWidth',1,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',6)
plot(1:time_spacing+1,NRMSE_NLM,'b-s','LineWidth',1,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',6)
hold off;
box on

xlabel('$\tau/\delta t$','Interpreter','latex')
ylabel('$\epsilon(\tau)$', 'interpreter', 'latex')

xlim([1 time_spacing+1]);
% ylim([0 0.2]);
% ylim([0 0.5]);

set(gca, 'YTickMode','manual');
set(gca, 'YTick', 0:0.1:0.2);
set (gca, 'YTickLabel', {'0.0','0.1','0.2'},'FontSize',fsize)
% set(gca, 'YTick', 0:0.1:0.5);
% set (gca, 'YTickLabel', {'0.0','0.1','0.2','0.3','0.4','0.5'},'FontSize',fsize)

set(gca, 'XTickMode','manual');
set(gca, 'XTick', 1:2:time_spacing+1);
set (gca, 'XTickLabel', {'0','P/2Q','P/Q'},'FontSize',fsize)

leg =legend ({'$\mathbf{I}_s \mathbf{y}$','$ \mathbf{I}_t \mathbf{x}$','$Fusion$',},'interpreter', 'latex','location','south'); 

set(leg,'FontSize',fsize-2);
legend boxoff
box on

% filename_fig=strcat('./figures/NRMSE_timespacing_',num2str(time_spacing,'%.2d'),'_spacespacing_',num2str(space_spacing,'%.2d'));
% export_fig(filename_fig,'-eps','-q101','-a4','-nocrop');
% close()