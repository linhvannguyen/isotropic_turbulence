load /data/ISOTROPIC/regression/RR_crossvalidation.mat;
train_MSE_mean = mean(-train_MSE, 2);
train_MSE_std = std(-train_MSE, 0, 2);
test_MSE_mean = mean(-test_MSE, 2);
test_MSE_std = std(-test_MSE, 0, 2);
N=size(lambdas_range,1);
%%
fsize=32;
fname='CMU Serif';

h=figure;
set(h, 'Position', [200 200 800 800]);
set(h,'color','w')

% Change default axes fonts.
set(0,'DefaultAxesFontName', fname)
set(0,'DefaultAxesFontSize', fsize)

% Change default text fonts.
set(0,'DefaultTextFontname', fname)
set(0,'DefaultTextFontSize', fsize)

% define how the figure inside the plot appear on the paper
set(gcf,'Units','normal');
set(gca,'Position',[0.2 0.2 0.75 0.75]); % [x_leftlowcorner y_leftlowcorner width height]
set(gcf,'Units','pixels');

hold on
shadedErrorBar(log10(lambdas_range),train_MSE_mean,train_MSE_std,{'b-','LineWidth',3},1);
shadedErrorBar(log10(lambdas_range),test_MSE_mean,test_MSE_std,{'r-','LineWidth',3},1);
hold off
box on

xlabel('$log_{10}(\lambda)$','Interpreter','latex')
ylabel('$\bar{\epsilon}$', 'interpreter', 'latex')

xlim([-2 4]);
ylim([0 0.6]);

set(gca, 'YTickMode','manual');
set(gca, 'YTick', 0:0.2:0.6);
set (gca, 'YTickLabel', {'0.0','0.2','0.4','0.6'},'FontSize',fsize)

set(gca, 'XTickMode','manual');
set(gca, 'XTick', -2:2:4);
set (gca, 'XTickLabel', {'-2','0','2','4'},'FontSize',fsize)


export_fig('./figures/RR_validationcurve', '-a1','-q101','-eps','-painters'); % leave '-painters' can cause quality problems
close()
