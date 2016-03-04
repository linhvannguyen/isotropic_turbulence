% clear all; close all; clc;
addpath('./funcs');

%% INITIAL PARAMETERS
filename_ref='/data/ISOTROPIC/data/data_downsampled4.nc';
nc = netcdf(filename_ref,'r');
Nt = nc('Nt').itsDimsize;
N_HR = nc('Nx').itsDimsize;
close(nc);

k_max=N_HR/2;
k_max_real=k_max*2/3; % ???

k_1D_HR=[0:k_max -k_max+1:1:-1];
[KZS,KYS] = meshgrid(k_1D_HR,k_1D_HR); 

%% 2D in space
% space_spacing_all=[2,3,4,6];
% 
% % Estimate loss of energy in 2D space
% ener_loss_space=zeros(numel(space_spacing_all),1);
% 
% nc=netcdf(filename_ref,'r');  
% for k=1:numel(space_spacing_all)  
%     space_spacing=space_spacing_all(k)
% 
%     u_HR_all=zeros(Nt,N_HR,N_HR,N_HR);
%     u_HR_fil_all=zeros(Nt,N_HR,N_HR,N_HR);
%     for t=1:Nt
%         for j=1:N_HR
%             u_HR=nc{'velocity_x'}(t,:,:,j);
%             u_HR_fil=filter_2D(u_HR, space_spacing, space_spacing);   
%             
%             u_HR_all(t,:,:,j)=u_HR;
%             u_HR_fil_all(t,:,:,j)=u_HR_fil;           
%         end
%         
%     end
%     RMS_ref=u_HR_all(:)'*u_HR_all(:);
%     RMS_fil=u_HR_fil_all(:)'*u_HR_fil_all(:);
%     ener_loss_space(k) = (RMS_ref-RMS_fil)/RMS_ref;
%     
%     clearvars u_HR_all u_HR_fil_all u_HR u_HR_fil;
% end
% close(nc)

%% 1D in time
time_spacing_all=[2,3,4,6, 8, 12];

% Estimate loss of energy in 1D streamwise
ener_loss_time=zeros(numel(time_spacing_all),1);
nc=netcdf(filename_ref,'r');
for k=1:numel(time_spacing_all) 
    time_spacing=time_spacing_all(k)

    u_HR_all=zeros(Nt,N_HR,N_HR,N_HR);
    u_HR_fil_all=zeros(Nt,N_HR,N_HR,N_HR);

    for t=1:Nt
        for m=1:N_HR
            for n=1:N_HR
                u_HR=nc{'velocity_x'}(t,m,n,:);
                u_HR_fil=filter_1D(u_HR, time_spacing);   
                u_HR_all(t,m,n,:)=u_HR;
                u_HR_fil_all(t,m,n,:)=u_HR_fil;
            end
        end
    end
    
    RMS_ref=u_HR_all(:)'*u_HR_all(:);
    RMS_fil=u_HR_fil_all(:)'*u_HR_fil_all(:);
    ener_loss_time(k) = (RMS_ref-RMS_fil)/RMS_ref;
    
    clearvars u_HR_all u_HR_fil_all u_HR u_HR_fil;   
end
close(nc) 

save energyloss_variousratios.mat space_spacing_all time_spacing_all ener_loss_space ener_loss_time 

%% plot sample
% velocity_x_interp=imresize(imresize(single(u_HR), 1/scale_factor_space,'cubic'), scale_factor_space,'cubic');
% fig1=figure(); imagesc(u_HR); caxis([-3,3]);
% fig2=figure(); imagesc(u_HR_fil); caxis([-3,3]);
% fig3=figure(); imagesc(velocity_x_interp); caxis([-3,3]);
% fig4=figure(); imagesc(u_HR-u_HR_fil); caxis([-3,3]);
% 
% fig5=figure();
% hold on
% plot(u_HR_fil,'r-')
% plot(u_HR,'k-')
% hold off;

