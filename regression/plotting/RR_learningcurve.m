load /data/ISOTROPIC/regression/RR_learningcurve.mat;
train_MSE_mean = mean(-train_MSE, 2);
train_MSE_std = std(-train_MSE, 0, 2);
test_MSE_mean = mean(-test_MSE, 2);
test_MSE_std = std(-test_MSE, 0, 2);


fsize=32;
fname='CMU Serif';

%% Ridge Regression
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
shadedErrorBar(train_sizes,train_MSE_mean,train_MSE_std,{'b-','LineWidth',3},1);
shadedErrorBar(train_sizes,test_MSE_mean,test_MSE_std,{'r-','LineWidth',3},1);
hold off
box on

xlabel('$Ntr$','Interpreter','latex')
ylabel('$\bar{\epsilon}$', 'interpreter', 'latex')

xlim([0 800]);
ylim([0 0.6]);

set(gca, 'YTickMode','manual');
set(gca, 'YTick', 0:0.2:0.6);
set (gca, 'YTickLabel', {'0.0','0.2','0.4','0.6'},'FontSize',fsize)

set(gca, 'XTickMode','manual');
set(gca, 'XTick', 0:200:800);
set (gca, 'XTickLabel', {'0','200','400','600','800'},'FontSize',fsize)


export_fig('./figures/RR_learningcurve', '-a1','-q101','-eps','-painters'); % leave '-painters' can cause quality problems
close()


%% Kernel Ridge Regression
load /data/ISOTROPIC/regression/KRR_learningcurve.mat;
train_MSE_mean = mean(-train_MSE, 2);
train_MSE_std = std(-train_MSE, 0, 2);
test_MSE_mean = mean(-test_MSE, 2);
test_MSE_std = std(-test_MSE, 0, 2);

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
shadedErrorBar(train_sizes,train_MSE_mean,train_MSE_std,{'b-','LineWidth',3},1);
shadedErrorBar(train_sizes,test_MSE_mean,test_MSE_std,{'r-','LineWidth',3},1);
hold off
box on

xlabel('$Ntr$','Interpreter','latex')
ylabel('$\bar{\epsilon}$', 'interpreter', 'latex')

xlim([0 800]);
ylim([0 0.6]);

set(gca, 'YTickMode','manual');
set(gca, 'YTick', 0:0.2:0.6);
set (gca, 'YTickLabel', {'0.0','0.2','0.4','0.6'},'FontSize',fsize)

set(gca, 'XTickMode','manual');
set(gca, 'XTick', 0:200:800);
set (gca, 'XTickLabel', {'0','200','400','600','800'},'FontSize',fsize)


export_fig('./figures/KRR_learningcurve', '-a1','-q101','-eps','-painters'); % leave '-painters' can cause quality problems
close()
