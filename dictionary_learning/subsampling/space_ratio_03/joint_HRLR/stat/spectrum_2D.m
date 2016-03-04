clear all; close all; clc;
path(path,'/home/nguyen/Documents/Codes/isotropic/DictionaryLearning/final/funcs')

%% INITIAL PARAMETERS
Nh=96; % Number of grid-point in HR field
scale=3;
Nl=Nh/scale;

size_l=8; 
patchsize_h = size_l*scale; % size of HR patches
patchsize_l = 8; % 4x4 LR patches

k_max=(Nh/2)*2/3;
k_cutoff = k_max/scale;

[gridz_h,gridy_h]=meshgrid(1:Nh+patchsize_h,1:Nh+patchsize_h);
gridz_l=gridz_h(2:scale:end,2:scale:end);
gridy_l=gridy_h(2:scale:end,2:scale:end);

%% Load 3D field
E_HR_ave=0; E_LR_interp_ave=0; E_SRrec_features_ave=0; E_SRrec_ave=0; 
k_1D_HR=[0:Nh/2 -Nh/2+1:1:-1];
k_1D_LR_interp=[0:Nh/2 -Nh/2+1:1:-1];
k_1D_LR=[0:Nh/(2*scale) -Nh/(2*scale)+1:1:-1];
k_1D_SRrec=[0:Nh/2 -Nh/2+1:1:-1];

nc1 = netcdf('/data/ISOTROPIC/NS3D_DNS_FORCE_384/RECONSTRUCTINGPLANES_downsampled4_HR.nc','r');
nc2 = netcdf('/data/ISOTROPIC/DICTIONARYLEARNING/final/downsample4_3/joint_HRLR/SRrec/SR_patchesHRLR_patchsize08_all.nc','r');
nc3 = netcdf('/data/ISOTROPIC/DICTIONARYLEARNING/final/downsample4_3/joint_HRLRfeatures/SRrec/SR_patchesLRfeatures_patchsize08_all.nc','r');

ttds=1:nc1('Nt').itsDimsize;

for t=1:numel(ttds)
    t
    X_HR_org=nc1{'velocity_x_HR'}(ttds(t),:,:);
    X_SRrec=nc2{'X_HR_rec'}(ttds(t),:,:);    
    X_SRrec_features=nc3{'X_HR_rec'}(ttds(t),:,:);
    
    % LR field
    X_LR=resize_nguyen(X_HR_org, 1/scale,'bicubic');
    X_LR_enlarged = [X_LR(Nl-patchsize_l/2+1:Nl,Nl-patchsize_l/2+1:Nl),X_LR(Nl-patchsize_l/2+1:Nl,:),X_LR(Nl-patchsize_l/2+1:Nl,1:patchsize_l/2);...
        X_LR(:,Nl-patchsize_l/2+1:Nl),X_LR,X_LR(:,1:patchsize_l/2);...
        X_LR(1:patchsize_l/2,Nl-patchsize_l/2+1:Nl),X_LR(1:patchsize_l/2,:),X_LR(1:patchsize_l/2,1:patchsize_l/2)];

    X_LR_interp=resize_nguyen(X_LR_enlarged, scale,'bicubic');

    X_LR_interp=X_LR_interp(patchsize_h/2+1:end-patchsize_h/2,patchsize_h/2+1:end-patchsize_h/2);

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
    
    %% 2D spectra of interp LR 
    F_LR_interp = fftn(X_LR_interp);

    k_2D_LR_interp=zeros(Nh,1);
    E_LR_interp=zeros(Nh,1);
    for m=1:Nh
        for n=1:Nh
            temp=sqrt(k_1D_LR_interp(m)^2+k_1D_LR_interp(n)^2);
            k_2D_LR_interp(round(temp)+1,1)=round(temp); % first wave number is zero
            E_LR_interp(round(temp)+1,1) = E_LR_interp(round(temp)+1,1) + abs(1/Nh^2*F_LR_interp(m,n))^2;
        end
    end
    E_LR_interp_ave=E_LR_interp_ave+1/numel(ttds)*E_LR_interp;
    
    %% 2D spectra of SR_rec
    F_SRrec_features = fftn(X_SRrec_features);

    k_2D_SRrec_features=zeros(Nh,1);
    E_SRrec_features=zeros(Nh,1);
    for m=1:Nh
        for n=1:Nh
            temp=sqrt(k_1D_SRrec(m)^2+k_1D_SRrec(n)^2);
            k_2D_SRrec_features(round(temp)+1,1)=round(temp); % first wave number is zero
            E_SRrec_features(round(temp)+1,1) = E_SRrec_features(round(temp)+1,1) + abs(1/Nh^2*F_SRrec_features(m,n))^2;
        end
    end
    E_SRrec_features_ave=E_SRrec_features_ave+1/numel(ttds)*E_SRrec_features;    
    
    
    %% 2D spectra of SR_rec
    F_SRrec = fftn(X_SRrec);

    k_2D_SRrec=zeros(Nh,1);
    E_SRrec=zeros(Nh,1);
    for m=1:Nh
        for n=1:Nh
            temp=sqrt(k_1D_SRrec(m)^2+k_1D_SRrec(n)^2);
            k_2D_SRrec(round(temp)+1,1)=round(temp); % first wave number is zero
            E_SRrec(round(temp)+1,1) = E_SRrec(round(temp)+1,1) + abs(1/Nh^2*F_SRrec(m,n))^2;
        end
    end
    E_SRrec_ave=E_SRrec_ave+1/numel(ttds)*E_SRrec;    
    
end
close(nc1); close(nc2) ; close(nc3)  

%% Plot
% -5/3 line
xline=[3,10];
yline=3*1e0*[exp((-5/3)*log(3)),exp((-5/3)*log(10))];

%
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
h1=loglog(k_2D_HR(2:k_max),3/2*E_HR_ave(2:k_max),'k-','LineWidth',2); 
hold on;
h2=loglog(k_2D_LR_interp(2:k_max),3/2*E_LR_interp_ave(2:k_max),'b-','LineWidth',2); 
h3=loglog(k_2D_LR_interp(2:k_max),3/2*E_SRrec_ave(2:k_max),'m-','LineWidth',2);
h4=loglog(k_2D_LR_interp(2:k_max),3/2*E_SRrec_features_ave(2:k_max),'r-','LineWidth',2);

% h4=loglog(k_2D_LR(2:11),3/2*E_LR_ave(2:11),'g-s','LineWidth',2); 

plot([k_cutoff,k_cutoff] ,[10^-8,1],'r--','LineWidth',3)
loglog(xline,yline,'k-','LineWidth',2)

hold off

ylim([10^-5,1]);
xlim([1,50]);

text(7,0.3,'-5/3','HorizontalAlignment','right','FontSize',fsize+2)

xlabel('k'); ylabel('E(k)');

leg=legend([h1 h2 h3 h4],{'Reference','Interpolation','SR','SR with features'},'location','northeast');
set(leg,'FontSize',fsize-2);
legend boxoff

filename=strcat('./figures/spectra2d');
export_fig(filename,'-eps','-q101','-a4','-nocrop');


%% 
energy_loss_smallscale=sum(E_HR_ave(k_cutoff+1:k_max))/sum(E_HR_ave(2:k_max)) % if using ideal filter (loss of small scale)
energy_loss_LR=(sum(E_HR_ave(2:k_max))-sum(E_LR_interp_ave(2:k_max)))/sum(E_HR_ave(2:k_max)) % if using bicubic interpolation
energy_loss_SRrec=(sum(E_HR_ave(2:k_max))-sum(E_SRrec_ave(2:k_max)))/sum(E_HR_ave(2:k_max))
energy_loss_SRrec_features=(sum(E_HR_ave(2:k_max))-sum(E_SRrec_features_ave(2:k_max)))/sum(E_HR_ave(2:k_max))