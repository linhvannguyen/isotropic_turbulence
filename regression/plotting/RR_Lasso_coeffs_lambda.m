clear all; close all; clc;

%% All constants
Nh = 96;
Nt = 37;
sspacing = 3;
tspacing = 4;

HTLS_sknots = 1:sspacing:Nh;
HTHS_sknots = 1:1:Nh;
LTHS_tknots = 0:tspacing:Nh;
Nl = numel(HTLS_sknots);
Ns = numel(LTHS_tknots);

Dh = Nh*Nh;
Dl = Nl*Nl;

N = Nt*Ns;

load /data/ISOTROPIC/regression/RR_Lasso_coeffs_lambda.mat;

%%
fsize=25;
fname='CMU Serif';

%% RR
h=figure;
set(h, 'Position', [200 200 1200 500]);
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


plot(log10(lambdas_RR(2:end)), coefs_RR_lambdas(2:end,:),'LineWidth',2)

xlabel('$log_{10}(\lambda)$','Interpreter','latex')
ylabel('$\beta_{ij}$', 'interpreter', 'latex')

xlim([-3 6]);
ylim([-0.2 1]);

set(gca, 'YTickMode','manual');
set(gca, 'YTick', 0:0.5:1);
set (gca, 'YTickLabel', {'0.0','0.5','1.0'},'FontSize',fsize)

set(gca, 'XTickMode','manual');
set(gca, 'XTick', -3:3:6);

box on

export_fig('./figures/RR_coeffs_lambda','-eps','-q101','-a4','-painters');
close()


%% RR
h=figure;
set(h, 'Position', [200 200 1200 500]);
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


plot(log10(lambdas_Lasso(2:end)), coefs_Lasso_lambdas(2:end,:),'LineWidth',2)

xlabel('$log_{10}(\lambda)$','Interpreter','latex')
ylabel('$\beta_{ij}$', 'interpreter', 'latex')

xlim([-3 6]);
ylim([-0.2 1]);

set(gca, 'YTickMode','manual');
set(gca, 'YTick', 0:0.5:1);
set (gca, 'YTickLabel', {'0.0','0.5','1.0'},'FontSize',fsize)

set(gca, 'XTickMode','manual');
set(gca, 'XTick', -3:3:6);

box on

export_fig('./figures/Lasso_coeffs_lambda','-eps','-q101','-a4','-painters');
close()