load /data/ISOTROPIC/regression/RR_validationcurve.mat;
train_MSE_mean = mean(-train_MSE, 2);
train_MSE_std = std(-train_MSE, 0, 2);
test_MSE_mean = mean(-test_MSE, 2);
test_MSE_std = std(-test_MSE, 0, 2);
N=size(lambdas_range,1);
%%
fsize=22;
fname='CMU Serif';

h=figure;
set(h, 'Position', [200 200 800 600]);
set(h,'color','w')

% Change default axes fonts.
set(0,'DefaultAxesFontName', fname)
set(0,'DefaultAxesFontSize', fsize)

% Change default text fonts.
set(0,'DefaultTextFontname', fname)
set(0,'DefaultTextFontSize', fsize)

% define how the figure inside the plot appear on the paper
set(gcf,'Units','normal');
set(gca,'Position',[0.15 0.2 0.8 0.75]); % [x_leftlowcorner y_leftlowcorner width height]
set(gcf,'Units','pixels');

hold on
shadedErrorBar(log10(lambdas_range(2:end)),train_MSE_mean(2:end),train_MSE_std(2:end),{'r-','LineWidth',2},1);
shadedErrorBar(log10(lambdas_range(2:end)),test_MSE_mean(2:end),test_MSE_std(2:end),{'b-','LineWidth',2},1);
hold off
box on

xlabel('$log_{10}(\lambda)$','Interpreter','latex')
ylabel('$\bar{\epsilon}$', 'interpreter', 'latex')

xlim([0 6]);
ylim([0 1.2]);

set(gca, 'YTickMode','manual');
set(gca, 'YTick', 0:0.4:1.2);
set (gca, 'YTickLabel', {'0.0','0.4','0.8','1.2'},'FontSize',fsize)

set(gca, 'XTickMode','manual');
set(gca, 'XTick', 0:2:6);
set (gca, 'XTickLabel', {'0','2','4','6'},'FontSize',fsize)


export_fig('./figures/RR_validationcurve', '-a1','-q101','-eps','-painters'); % leave '-painters' can cause quality problems
close()
