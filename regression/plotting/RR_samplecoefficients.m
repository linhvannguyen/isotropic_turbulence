load /data/ISOTROPIC/regression/RR_sample_coefficients.mat;

fsize=28;
fname='CMU Serif';

% h=figure;
% set(h, 'Position', [200 200 800 800]);
% set(h,'color','w')
% 
% % Change default axes fonts.
% set(0,'DefaultAxesFontName', fname)
% set(0,'DefaultAxesFontSize', fsize)
% 
% % Change default text fonts.
% set(0,'DefaultTextFontname', fname)
% set(0,'DefaultTextFontSize', fsize)
% 
% % define how the figure inside the plot appear on the paper
% set(gcf,'Units','normal');
% set(gca,'Position',[0.1 0.1 0.7 0.7]); % [x_leftlowcorner y_leftlowcorner width height]
% set(gcf,'Units','pixels');
% 
% imagesc(coefs);
% 
% caxis([-0.05,0.4]);
% cb=colorbar;
% set(cb,'position',[0.825 0.1 0.03 0.7]) 
% set(cb, 'YTick', 0:0.1:0.4);
% set (cb, 'YTickLabel', {'0.0','0.1','0.2','0.3','0.4'})
% 
% 
% axis off;
% box on;
% 
% export_fig('./figures/RR_samplecoefficients', '-a1','-q101','-eps','-painters'); % leave '-painters' can cause quality problems
% close()


h=figure;
set(h, 'Position', [200 200 1600 800]);
set(h,'color','w')

% Change default axes fonts.
set(0,'DefaultAxesFontName', fname)
set(0,'DefaultAxesFontSize', fsize)

% Change default text fonts.
set(0,'DefaultTextFontname', fname)
set(0,'DefaultTextFontSize', fsize)

% define how the figure inside the plot appear on the paper
set(gcf,'Units','normal');
set(gca,'Position',[0.1 0.1 0.7 0.7]); % [x_leftlowcorner y_leftlowcorner width height]
set(gcf,'Units','pixels');

ax(1)=subplot(2,2,1,'position',[0.1 0.1 0.35 0.7]); % top left
hold on
imagesc(coefs_LSE);
axis off; box on;
caxis([-0.2,0.4]);


ax(3)=axes('position',[0.45 0.1 0.35 0.7]); % top right
hold on
imagesc(coefs_RR);

caxis([-0.2,0.4]);
cb=colorbar;
set(cb,'position',[0.8 0.105 0.015 0.67]) 
set(cb, 'YTick', -0.2:0.2:0.4);
set (cb, 'YTickLabel', {'-0.2','0.0','0.2','0.4'})
axis off; box on;

export_fig('./figures/RR_LSE_coefficients', '-a1','-q101','-eps','-painters'); % leave '-painters' can cause quality problems
close()
