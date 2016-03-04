clear all; close all; clc;
path(path,'/home/nguyen/Documents/Codes/isotropic/DictionaryLearning/final/funcs')

%% INITIAL PARAMETERS
Nh=96; % Number of grid-point in HR field
scale=3;
Nl=Nh/scale;

size_l=8; 
patchsize_h = size_l*scale; % size of HR patches
patchsize_l = patchsize_h; % 4x4 LR patches

k_max=(Nh/2)*2/3;
k_cutoff = k_max/scale;

[gridz_h,gridy_h]=meshgrid(1:Nh+patchsize_h,1:Nh+patchsize_h);
gridz_l=gridz_h(2:scale:end,2:scale:end);
gridy_l=gridy_h(2:scale:end,2:scale:end);

%% Load 3D field
E_error_SRrec_features_ave=0; 
E_error_SRrec_ave=0;
E_error_LR_interp_ave=0; 
k_1D_error_SRrec=[0:Nh/2 -Nh/2+1:1:-1];
k_1D_error_LR_interp=[0:Nh/2 -Nh/2+1:1:-1];

E_HR_ave=0;
k_1D_HR=[0:Nh/2 -Nh/2+1:1:-1];


nc1 = netcdf('/data/ISOTROPIC/NS3D_DNS_FORCE_384/RECONSTRUCTINGPLANES_downsampled4_HR.nc','r');
nc2 = netcdf('/data/ISOTROPIC/DICTIONARYLEARNING/final/downsample4_3/joint_HRLR/SRrec/SR_patchesHRLR_patchsize08_all.nc','r');
nc3 = netcdf('/data/ISOTROPIC/DICTIONARYLEARNING/final/downsample4_3/joint_HRLRfeatures/SRrec/SR_patchesLRfeatures_patchsize08_all.nc','r');
ttds=1:nc1('Nt').itsDimsize;

for t=1:numel(ttds)
    t
    % Load data
    X_HR_org=nc1{'velocity_x_HR'}(ttds(t),:,:);
        
    % LR field
    X_LR=resize_nguyen(X_HR_org, 1/scale,'bicubic');
    X_LR_enlarged = [X_LR(Nl-size_l/2+1:Nl,Nl-size_l/2+1:Nl),X_LR(Nl-size_l/2+1:Nl,:),X_LR(Nl-size_l/2+1:Nl,1:size_l/2);...
        X_LR(:,Nl-size_l/2+1:Nl),X_LR,X_LR(:,1:size_l/2);...
        X_LR(1:size_l/2,Nl-size_l/2+1:Nl),X_LR(1:size_l/2,:),X_LR(1:size_l/2,1:size_l/2)];
    
    X_LR_interp=resize_nguyen(X_LR_enlarged, scale,'bicubic');
%     X_LR_interp=interp2(gridz_l,gridy_l,X_LR_enlarged,gridz_h,gridy_h,'spline');
    
    X_LR_interp=X_LR_interp(patchsize_h/2+1:end-patchsize_h/2,patchsize_h/2+1:end-patchsize_h/2);
    
    %% ERROR
    error_SRrec_features=X_HR_org-nc3{'X_HR_rec'}(ttds(t),:,:);
    error_SRrec=X_HR_org-nc2{'X_HR_rec'}(ttds(t),:,:);
    error_LR_interp=X_HR_org-X_LR_interp;

    %% 2D spectra of full resolution
    F_HR = fftn(X_HR_org);

    k_2D_HR=zeros(Nh,1);
    E_HR=zeros(Nh,1);

    for m=1:Nh
        for n=1:Nh
            temp=sqrt(k_1D_HR(m)^2+k_1D_HR(n)^2);
            k_2D_HR(round(temp)+1,1)=round(temp); % first wave number is zero
            E_HR(round(temp)+1,1) = E_HR(round(temp)+1,1) + abs(1/Nh^2*F_HR(m,n))^2;
        end
    end
    E_HR_ave=E_HR_ave+1/numel(ttds)*E_HR;

    %% 2D spectra of SR
    F_error_SRrec_features = fftn(error_SRrec_features);
    F_error_SRrec = fftn(error_SRrec);
    
    k_2D_errorSR=zeros(Nh,1);
    E_error_SRrec_features=zeros(Nh,1);
    E_error_SRrec=zeros(Nh,1);
    
    for m=1:Nh
        for n=1:Nh
            temp=sqrt(k_1D_error_SRrec(m)^2+k_1D_error_SRrec(n)^2);
            k_2D_errorSR(round(temp)+1,1)=round(temp); % first wave number is zero
            E_error_SRrec_features(round(temp)+1,1) = E_error_SRrec_features(round(temp)+1,1) + abs(1/Nh^2*F_error_SRrec_features(m,n))^2;
            E_error_SRrec(round(temp)+1,1) = E_error_SRrec(round(temp)+1,1) + abs(1/Nh^2*F_error_SRrec(m,n))^2;
        end
    end
    E_error_SRrec_features_ave=E_error_SRrec_features_ave+1/numel(ttds)*E_error_SRrec_features;
    E_error_SRrec_ave=E_error_SRrec_ave+1/numel(ttds)*E_error_SRrec;
    
    %% 2D spectra of interp LR 
    F_error_LR_interp = fftn(error_LR_interp);

    k_2D_errorLR_interp=zeros(Nh,1);
    E_error_LR_interp=zeros(Nh,1);
    for m=1:Nh
        for n=1:Nh
            temp=sqrt(k_1D_error_LR_interp(m)^2+k_1D_error_LR_interp(n)^2);
            k_2D_errorLR_interp(round(temp)+1,1)=round(temp); % first wave number is zero
            E_error_LR_interp(round(temp)+1,1) = E_error_LR_interp(round(temp)+1,1) + abs(1/Nh^2*F_error_LR_interp(m,n))^2;
        end
    end
    E_error_LR_interp_ave=E_error_LR_interp_ave+1/numel(ttds)*E_error_LR_interp;
end
close(nc1); close(nc3); 




%% Plot error spectra
fsize=25;
fname='CMU Serif';
set(0,'DefaultAxesFontName', fname)
set(0,'DefaultAxesFontSize', fsize)
set(0,'DefaultTextFontname', fname)
set(0,'DefaultTextFontSize', fsize)

fig1=figure(); 
set(gcf, 'Position', [400 100 1000 800]);
set(fig1, 'Color', 'w');

h1=loglog(k_2D_errorLR_interp(1:k_max),E_error_LR_interp_ave(1:k_max),'k-s','LineWidth',2); 
hold on;
h2=loglog(k_2D_errorSR(1:k_max),E_error_SRrec_ave(1:k_max),'r-s','LineWidth',2); 
h3=loglog(k_2D_errorSR(1:k_max),E_error_SRrec_features_ave(1:k_max),'b-s','LineWidth',2); 

plot([k_cutoff,k_cutoff] ,[8*10^-7,1.10^-2],'r--','LineWidth',3)

hold off
  
ylim([8*10^-7,10^-2]);
xlim([1,30]);

xlabel('k'); ylabel('E(k)');

leg=legend([h1 h2 h3],{'Interpolation','SR','SR with features'},'location','northwest');
set(leg,'FontSize',fsize-2);
legend boxoff

filename=strcat('./figures/spectra2d_error');
export_fig(filename,'-eps','-q101','-a4','-nocrop');


%% Plot normalized
fsize=25;
fname='CMU Serif';
set(0,'DefaultAxesFontName', fname)
set(0,'DefaultAxesFontSize', fsize)
set(0,'DefaultTextFontname', fname)
set(0,'DefaultTextFontSize', fsize)

nc1=netcdf('/data/ISOTROPIC/NS3D_DNS_FORCE_384/STAT.nc','r');
energy_spectra = nc1{'energy_spectra'}(:,:);
close(nc1)

fig1=figure(); 
set(gcf, 'Position', [400 100 1000 800]);
set(fig1, 'Color', 'w');

h1=loglog(k_2D_errorLR_interp(1:k_max),E_error_LR_interp_ave(1:k_max)./E_HR_ave(1:k_max),'k-s','LineWidth',2); 
hold on;
h2=loglog(k_2D_errorSR(1:k_max),E_error_SRrec_ave(1:k_max)./E_HR_ave(1:k_max),'r-s','LineWidth',2); 
h3=loglog(k_2D_errorSR(1:k_max),E_error_SRrec_features_ave(1:k_max)./E_HR_ave(1:k_max),'b-s','LineWidth',2); 

plot([k_cutoff,k_cutoff] ,[10^-6,2],'r--','LineWidth',3)

hold off
  
ylim([10^-6,2]);
xlim([0,40]);

xlabel('k'); %ylabel('E(k)');

leg=legend([h1 h2 h3],{'Interpolation','SR','SR with features'},'location','northwest');
set(leg,'FontSize',fsize-2);
legend boxoff

filename=strcat('./figures/spectra2d_error_normalized');
export_fig(filename,'-eps','-q101','-a4','-nocrop');




%% Plot error accumulated
cumsum_error_LR_interp=cumsum(E_error_LR_interp_ave(1:k_max))./sum(E_error_LR_interp_ave(1:k_max));
cumsum_error_SRrec=cumsum(E_error_SRrec_ave(1:k_max))./sum(E_error_LR_interp_ave(1:k_max));
cumsum_error_SRrec_features=cumsum(E_error_SRrec_features_ave(1:k_max))./sum(E_error_LR_interp_ave(1:k_max));

fsize=25;
fname='CMU Serif';
set(0,'DefaultAxesFontName', fname)
set(0,'DefaultAxesFontSize', fsize)
set(0,'DefaultTextFontname', fname)
set(0,'DefaultTextFontSize', fsize)

nc1=netcdf('/data/ISOTROPIC/NS3D_DNS_FORCE_384/STAT.nc','r');
energy_spectra = nc1{'energy_spectra'}(:,:);
close(nc1)

fig1=figure(); 
set(gcf, 'Position', [400 100 1000 800]);
set(fig1, 'Color', 'w');

h1=plot(k_2D_errorSR(1:k_max),cumsum_error_LR_interp,'k-*','LineWidth',2); 
hold on;
h2=plot(k_2D_errorSR(1:k_max),cumsum_error_SRrec,'r-s','LineWidth',2); 
h3=plot(k_2D_errorSR(1:k_max),cumsum_error_SRrec_features,'b-s','LineWidth',2); 

plot([k_cutoff,k_cutoff] ,[0,1],'r--','LineWidth',3)

hold off
  
ylim([0,1]);
xlim([0,33]);

xlabel('k'); % ylabel('E(k)');

leg=legend([h1 h2 h3],{'Interpolation','SR','SR with features'},'location','northwest');
set(leg,'FontSize',fsize-2);
legend boxoff

filename=strcat('./figures/spectra2d_error_normalized_interp');
export_fig(filename,'-eps','-q101','-a4','-nocrop');

