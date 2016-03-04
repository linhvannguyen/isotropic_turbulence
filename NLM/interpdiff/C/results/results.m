% clear all; close all; 
% clc;

%% INITIAL PARAMS
space_spacing=4;
time_spacing=4;

filename_ref='/data/ISOTROPIC/data/data_downsampled4.nc';
filename_NLM='/data/ISOTROPIC/NLM/interpdiff/NLmean_propag2dirs_sspacing3_tspacing4_sim12_acc12_neighbor2_tau0100.nc';
% filename_NLM='/data/ISOTROPIC/NLM/interpdiff/NLM_interpdiff_simpatch0_accumpatch0_searchbox0_tau0100.nc';
% filename_NLM='/data/ISOTROPIC/NLM/interpdiff/NLmean_propag_sspacing3_tspacing8_sim12_acc12_neighbor1_tau0100.nc';
% filename_NLM='/data/ISOTROPIC/NLM/interpdiff/NLmean_interpdiff_propag2dirs_simpatch0_accumpatch0_searchbox0_tau0100.nc';
nc = netcdf(filename_ref,'r');
Nt = nc('Nt').itsDimsize;
Nz = nc('Nz').itsDimsize;
Ny = nc('Ny').itsDimsize; 
Nx = nc('Nx').itsDimsize;
close(nc)

LTHS_idt=1:time_spacing:Nx;
HTHS_idy=1:Ny; % row indices of HTHS in space 
HTHS_idz=1:Nz; % column indices of HTHS in space
HTLS_idt=1:1:Nx; % indices of all HTHS snapshots in time
HTLS_idy= 1:space_spacing:Ny; % row indices of HTLS in space
HTLS_idz= 1:space_spacing:Nz; % column indices of HTLS in space

NRMSE=@(x_ref,x_est) sqrt(sum((x_est(:)-x_ref(:)).^2))/sqrt(sum(x_ref(:).^2));

%% Space interp
border=6*space_spacing;
left=border; right = border-(Ny-HTLS_idy(end));
bottom=border; top = border-(Nz-HTLS_idz(end));

nc1=netcdf(filename_ref,'r');
nc2=netcdf(filename_NLM,'r');

for t=10:10
    fprintf('Estimating block: %.2d-th snapshot.\n', t);    
%     for blockid=1:numel(LTHS_idt)-1
    for blockid=1:23
        t_PIV_prev = LTHS_idt(blockid);
        t_PIV_after = LTHS_idt(blockid+1);
               
        for pos_t=1:time_spacing-1           
            t_current = t_PIV_prev + pos_t;
            x_ref = nc1{'velocity_x'}(t,HTHS_idz,HTHS_idy,t_current);
            [x_interp, x_diff] = interp_border (x_ref, space_spacing, border, 1);

%             x_rec = x_interp+nc2{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
            x_rec = x_interp+nc2{'x_rec_all'}(t,t_current,HTHS_idz,HTHS_idy);

            NRMSE1 = NRMSE(x_ref, x_interp);
            NRMSE2 = NRMSE(x_ref, x_rec);
            fprintf(['\n Reconstruction with improvement of NRMSE is ', num2str(100*(1-NRMSE2/NRMSE1),'%.3f'), ' percents \n'])
        end 
    end 
end 
close(nc1); 

%%
% imagesc(weights1)   