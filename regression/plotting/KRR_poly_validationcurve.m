clear all; clc; close all;
%% 
load /data/ISOTROPIC/regression/KRR_poly_validationcurve.mat;
train_MSE_degree_mean = mean(-train_MSE_degree, 2);
train_MSE_degree_std = std(-train_MSE_degree, 0, 2);
test_MSE_degree_mean = mean(-test_MSE_degree, 2);
test_MSE_degree_std = std(-test_MSE_degree, 0, 2);

num_degrees=size(degree_range,1);
%%
fsize=32;
fname='CMU Serif';

%% Gamma
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
shadedErrorBar(log10(degree_range),train_MSE_degree_mean,train_MSE_degree_std,{'b-','LineWidth',3},1);
shadedErrorBar(log10(degree_range),test_MSE_degree_mean,test_MSE_degree_std,{'r-','LineWidth',3},1);
hold off
box on

xlabel('$log_{10}(d)$','Interpreter','latex')
ylabel('$\bar{\epsilon}$', 'interpreter', 'latex')

xlim([0 1]);
ylim([0 0.8]);

set(gca, 'YTickMode','manual');
set(gca, 'YTick', 0:0.4:0.8);
set (gca, 'YTickLabel', {'0.0','0.4','0.8'},'FontSize',fsize)

set(gca, 'XTickMode','manual');
set(gca, 'XTick', 0:0.5:1);
set (gca, 'XTickLabel', {'0.0','0.5','1.0'},'FontSize',fsize)

export_fig('./figures/KRR_validationcurve_polydegree', '-a1','-q101','-eps','-painters'); % leave '-painters' can cause quality problems
close()