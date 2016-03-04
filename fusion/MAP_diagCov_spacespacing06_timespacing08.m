clear all; close all; clc;

%% Define file locations (to load or to save)
addpath('./funcs/');
space_spacing=06; % subsampling ration in space
time_spacing=08; % subsampling ration in time (from 40Hz to 4Hz)

filename_ref='/data/ISOTROPIC/data/data_downsampled4.nc';
filename_interp_space=strcat('/data/ISOTROPIC/fusion/Uinterp_spatialspacing_',num2str(space_spacing,'%.2d'),'.nc');
filename_interp_time=strcat('/data/ISOTROPIC/fusion/Uinterp_timespacing_',num2str(time_spacing,'%.2d'),'.nc');
filename_fusion=strcat('/data/ISOTROPIC/fusion/FusedData_diagCov_timespacing_',num2str(time_spacing,'%.2d'),'_spacespacing_',num2str(space_spacing,'%.2d'),'.nc');

nc = netcdf(filename_ref,'r');
Nt = nc('Nt').itsDimsize;
Nx = nc('Nx').itsDimsize;
Ny = nc('Ny').itsDimsize;
Nz = nc('Nz').itsDimsize;
close(nc);

%% Define parameters
HTHS_idy=1:Ny; % row indices of HTHS in space 
HTHS_idz=1:Nz; % column indices of HTHS in space
HTLS_idt=1:1:Nx; % indices of all HTHS snapshots in time

HTLS_idy= 1:space_spacing:Ny; % row indices of HTLS in space
HTLS_idz= 1:space_spacing:Nz; % column indices of HTLS in space

LTHS_idt=1:time_spacing:Nx;

%% ============================ Interpolation step ========================
%  ========================================================================
%  ========================================================================

%% Space interp
[gridz_HS, gridy_HS] = meshgrid(1:Nz,1:Ny);
left=4*space_spacing; right = 4*space_spacing-(Ny-HTLS_idy(end));
bottom=4*space_spacing; top = 4*space_spacing-(Nz-HTLS_idz(end));
[gridz_HS_enlarged, gridy_HS_enlarged] = meshgrid(1:Nz+left+right,1:Ny+bottom+top);
gridz_LS_enlarged=gridz_HS_enlarged(1:space_spacing:end,1:space_spacing:end);
gridy_LS_enlarged=gridy_HS_enlarged(1:space_spacing:end,1:space_spacing:end);

fprintf('\n ********** START SPATIAL INTERPOLATION ************ \n');

nc1=netcdf(filename_ref,'r');

nc2 = netcdf(filename_interp_space,'clobber');
nc2('Nt')=0; 
nc2('Nx')=Nx;
nc2('Ny')=Ny;
nc2('Nz')=Nz;
nc2{'Uinterp'}=ncfloat('Nt','Nx','Ny','Nz');

for t=1:Nt
    fprintf('Interpolating %.4d-th snapshot.\n', t);    
    for i=1:Nx
        % interp
        asnap_LTHS=nc1{'velocity_x'}(t,HTHS_idz,HTHS_idy,i);
        asnap_LTHS_enlarged=enlarge_2D(asnap_LTHS,left, right, bottom, top);
        
        asnap_LTLS=asnap_LTHS_enlarged(1:space_spacing:end,1:space_spacing:end);
        asnap_LTHS_interp=interp2(gridz_LS_enlarged, gridy_LS_enlarged, asnap_LTLS, gridz_HS_enlarged, gridy_HS_enlarged,'spline');
               
        % save file
        nc2{'Uinterp'}(t,:,:,i)=asnap_LTHS_interp(bottom+1:bottom+Ny,left+1:left+Nz);
    end
end 
close(nc1); close(nc2);

fprintf('\n ********** COMPLETED SPATIAL INTERPOLATION ************ \n');

%% Time interp
left=4*time_spacing; right = 4*time_spacing-(Nx-LTHS_idt(end));

HTHS_idt_enlarged = 1:Nx+left+right;
LTHS_idt_enlarged = 1:time_spacing:numel(HTHS_idt_enlarged);

fprintf('\n ********** START TIME INTERPOLATION ************ \n');
nc1=netcdf(filename_ref,'r');

nc2 = netcdf(filename_interp_time,'clobber');
nc2('Nt')=0; 
nc2('Nx')=Nx;
nc2('Ny')=Ny;
nc2('Nz')=Nz;
nc2{'Uinterp'}=ncfloat('Nt','Nx','Ny','Nz');

for t=1:Nt
    fprintf('Interpolating %.4d-th snapshot.\n', t);    
    PIV_sampled=nc1{'velocity_x'}(t,:,:,:);
    PIV_sampled=cat(3, PIV_sampled(:,:,Nx-left+1:1:Nx), PIV_sampled, PIV_sampled(:,:,1:1:right));
    
    PIV_interp=interp1(LTHS_idt_enlarged, permute(PIV_sampled(:,:,LTHS_idt_enlarged),[3 1 2]),HTHS_idt_enlarged,'spline');
    nc2{'Uinterp'}(t,:,:,:)=permute(PIV_interp(left+1:left+Nx,:,:),[2,3,1]);
end

close(nc1); close(nc2);

fprintf('\n ********** COMPLETED TIME INTERPOLATION ************ \n');






%% ======================= Fusion parameters estimation ===================
%  ========================================================================
%  ========================================================================
fprintf('\n ********** START ESTIMATING FUSION PARAMS ************ \n');

%% Compute Cns 
start=tic();
nc1 = netcdf(filename_ref,'r');
nc2 = netcdf(filename_interp_space,'r');
hs=nc1{'velocity_x'}(:,HTHS_idz,HTHS_idy,LTHS_idt)-nc2{'Uinterp'}(:,HTHS_idz,HTHS_idy,LTHS_idt);
close(nc1); close(nc2);
var_hs=hs.^2; clearvars hs;

S=squeeze(sum(var_hs,1));
S=squeeze(sum(S,3));
S=S./(Nt*numel(LTHS_idt));

fprintf('Learning Cns: %.2f seconds.\n', toc(start)); 

% average
S_ave=zeros(Ny,Nz);
for i=1:space_spacing
    for j=1:space_spacing
        S_ave(i:space_spacing:Ny,j:space_spacing:Nz) = mean(mean(S(i:space_spacing:Ny,j:space_spacing:Nz)));
    end  
end

%% Compute Cnt
start=tic(); 
T_ave=zeros(time_spacing,1);

nc1 = netcdf(filename_ref,'r');
nc2 = netcdf(filename_interp_time,'r'); 

for i=1:time_spacing
    ht=nc1{'velocity_x'}(:,HTLS_idz,HTLS_idy,i:time_spacing:Nx)-nc2{'Uinterp'}(:,HTLS_idz,HTLS_idy,i:time_spacing:Nx);
    T_ave(i)=mean(ht(:).^2);
end
close(nc1); close(nc2);


fprintf('\n ********** FINISH ESTIMATING FUSION PARAMS ************ \n');
 

%% =============================== Fusion =================================
%  ========================================================================
%  ========================================================================
fprintf('\n ********** START FUSION STEP ************ \n');
nc1 = netcdf(filename_interp_space,'r');
nc2 = netcdf(filename_interp_time,'r');
nc = netcdf(filename_fusion,'clobber');
nc('Nt')=0;
nc('Nx')=Nx;
nc('Ny')=Ny;
nc('Nz')=Nz;
nc{'Zhat_all'}=ncfloat('Nt','Nz','Ny','Nx');

for t=1:Nt
    fprintf('Fusion %.2d-th block.\n', t);
    Zhat=zeros(Nz,Ny,Nx);
    for i=1:Nx
        pos_t=rem(i-1,time_spacing)+1;

        % estimate weights
        if(pos_t==1)
            Wt=ones(Nz,Ny);
            Ws=zeros(Nz,Ny);
        else
            Wt=S_ave./(S_ave+T_ave(pos_t));
            Ws=T_ave(pos_t)./(S_ave+T_ave(pos_t));
        end

        % fusion
        IsY =  nc1{'Uinterp'}(t,:,:,i);
        ItX =  nc2{'Uinterp'}(t,:,:,i);

        Zhat(:,:,i) = Wt.*ItX+Ws.*IsY;
    end 
    nc{'Zhat_all'}(t,:,:,:) = Zhat;
end
close(nc1); close(nc2); close(nc);

fprintf('\n ********** FINISH FUSION STEP ************ \n');