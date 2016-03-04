load /data/ISOTROPIC/regression/KRR_validationcurve.mat;
train_MSE_gamma_mean = mean(-train_MSE_gamma, 2);
train_MSE_gamma_std = std(-train_MSE_gamma, 0, 2);
test_MSE_gamma_mean = mean(-test_MSE_gamma, 2);
test_MSE_gamma_std = std(-test_MSE_gamma, 0, 2);

train_MSE_lambda_mean = mean(-train_MSE_lambda, 2);
train_MSE_lambda_std = std(-train_MSE_lambda, 0, 2);
test_MSE_lambda_mean = mean(-test_MSE_lambda, 2);
test_MSE_lambda_std = std(-test_MSE_lambda, 0, 2);

num_lambdas=size(lambdas_range,1);
num_gammas=size(gammas_range,1);

%% Gamma
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
shadedErrorBar(log10(gammas_range),train_MSE_gamma_mean,train_MSE_gamma_std,{'r-','LineWidth',2},1);
shadedErrorBar(log10(gammas_range),test_MSE_gamma_mean,test_MSE_gamma_std,{'b-','LineWidth',2},1);
hold off
box on

xlabel('$log_{10}(\gamma)$','Interpreter','latex')
ylabel('$\bar{\epsilon}$', 'interpreter', 'latex')

xlim([-8 -3]);
ylim([0 0.6]);

set(gca, 'YTickMode','manual');
set(gca, 'YTick', 0:0.2:0.6);
set (gca, 'YTickLabel', {'0.0','0.2','0.4','0.6'},'FontSize',fsize)

set(gca, 'XTickMode','manual');
set(gca, 'XTick', -8:1:-3);
set (gca, 'XTickLabel', {'-8','-7','-6','-5','-4','-3'},'FontSize',fsize)

export_fig('./figures/KRR_validationcurve_gamma', '-a1','-q101','-eps','-painters'); % leave '-painters' can cause quality problems
close()


%% Lambda
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
shadedErrorBar(log10(lambdas_range(14:end)),train_MSE_lambda_mean(14:end),train_MSE_lambda_std(14:end),{'r-','LineWidth',2},1);
shadedErrorBar(log10(lambdas_range(14:end)),test_MSE_lambda_mean(14:end),test_MSE_lambda_std(14:end),{'b-','LineWidth',2},1);
hold off
box on

xlabel('$log_{10}(\lambda)$','Interpreter','latex')
ylabel('$\bar{\epsilon}$', 'interpreter', 'latex')

xlim([-7 -1]);
ylim([0 0.75]);

set(gca, 'YTickMode','manual');
set(gca, 'YTick', 0:0.25:0.75);
set (gca, 'YTickLabel', {'0.00','0.25','0.50','0.75'},'FontSize',fsize)

set(gca, 'XTickMode','manual');
set(gca, 'XTick', -7:1:-1);
set (gca, 'XTickLabel', {'-7','-6','-5','-4','-3','-2','-1'},'FontSize',fsize)

export_fig('./figures/KRR_validationcurve_lambda', '-a1','-q101','-eps','-painters'); % leave '-painters' can cause quality problems
close()