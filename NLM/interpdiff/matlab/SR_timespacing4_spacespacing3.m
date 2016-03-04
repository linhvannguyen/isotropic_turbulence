clear all; close all; 
clc;

%% INITIAL PARAMS
space_spacing=3;
time_spacing=4;

simpatch_haftsize = 6; % the half width of the patche
simpatch_fullsize = 2*simpatch_haftsize+1; % the full width.
accpatch_haftsize = simpatch_haftsize; 
accpatch_fullsize = 2*accpatch_haftsize + 1; 
neighbor_haftsize = 5; % the half width of the patch
neighbor_fullsize = 2*neighbor_haftsize+1; % the full width.
dim_simpatch=simpatch_fullsize*simpatch_fullsize;
dim_accpatch=accpatch_fullsize*accpatch_fullsize;
dim_neighbor=neighbor_fullsize*neighbor_fullsize;

acc_ids_1D = simpatch_haftsize+1-accpatch_haftsize:simpatch_haftsize+1+accpatch_haftsize;
acc_ids = reshape(1:dim_simpatch,[simpatch_fullsize,simpatch_fullsize]);
acc_ids = acc_ids(acc_ids_1D, acc_ids_1D);
acc_ids = acc_ids(:);

filename_ref='/data/ISOTROPIC/data/data_downsampled4.nc';
nc = netcdf(filename_ref,'r');
Nt = nc('Nt').itsDimsize;
Nz = nc('Nz').itsDimsize;
Ny = nc('Ny').itsDimsize; 
Nx = nc('Nx').itsDimsize;
close(nc)

dim_im=Ny*Nz;

tau = .1;

LTHS_idt=1:time_spacing:Nx;
HTHS_idy=1:Ny; % row indices of HTHS in space 
HTHS_idz=1:Nz; % column indices of HTHS in space
HTLS_idt=1:1:Nx; % indices of all HTHS snapshots in time
HTLS_idy= 1:space_spacing:Ny; % row indices of HTLS in space
HTLS_idz= 1:space_spacing:Nz; % column indices of HTLS in space

NRMSE=@(x_ref,x_est) sqrt(sum((x_est(:)-x_ref(:)).^2))/sqrt(sum(x_ref(:).^2));

%% Space interp
[gridz_HS, gridy_HS] = meshgrid(1:Nz,1:Ny);
left=4*space_spacing; right = 4*space_spacing-(Ny-HTLS_idy(end));
bottom=4*space_spacing; top = 4*space_spacing-(Nz-HTLS_idz(end));
[gridz_HS_enlarged, gridy_HS_enlarged] = meshgrid(1:Nz+left+right,1:Ny+bottom+top);
gridz_LS_enlarged=gridz_HS_enlarged(1:space_spacing:end,1:space_spacing:end);
gridy_LS_enlarged=gridy_HS_enlarged(1:space_spacing:end,1:space_spacing:end);

nc1=netcdf(filename_ref,'r');

tic()
[gridpatches_y, gridpatches_z] = grid_gen(int16(Ny),int16(Nz),int16(simpatch_haftsize),int16(neighbor_haftsize));
toc()
for t=1:1
    fprintf('Estimating block: %.2d-th snapshot.\n', t);    
%     for blockid=1:numel(LTHS_idt)-1
    for blockid=1:1
        t_PIV_prev = LTHS_idt(blockid);
        t_PIV_after = LTHS_idt(blockid+1);
        
        x_PIV_prev= nc1{'velocity_x'}(t,HTHS_idz,HTHS_idy,t_PIV_prev);
        x_PIV_interp_prev = interp_border (x_PIV_prev, space_spacing);
        x_PIV_diff_prev = x_PIV_prev - x_PIV_interp_prev;
        
        x_PIV_after= nc1{'velocity_x'}(t,HTHS_idz,HTHS_idy,t_PIV_after);
        x_PIV_interp_after = interp_border (x_PIV_after, space_spacing);
        x_PIV_diff_after = x_PIV_after - x_PIV_interp_after;
        
        for pos_t=1:time_spacing-1
            start=tic();
            
            t_current = t_PIV_prev + pos_t;
            x_current = nc1{'velocity_x'}(t,HTHS_idz,HTHS_idy,t_current);
            x_current_interp = interp_border(x_current, space_spacing);
          
            [x_rec1,weights1] = NLmean(x_current_interp, x_PIV_interp_prev, x_PIV_diff_prev, gridpatches_y, gridpatches_z, simpatch_haftsize,accpatch_haftsize, neighbor_haftsize, tau);
            [x_rec2,weights2] = NLmean(x_current_interp, x_PIV_interp_after, x_PIV_diff_after, gridpatches_y, gridpatches_z, simpatch_haftsize,accpatch_haftsize, neighbor_haftsize, tau);
            x_rec = x_current_interp + (x_rec1 + x_rec2)./(weights1+weights2);
%             x_rec = x_current_interp + x_rec1./weights1;
            
            NRMSE1 = NRMSE(x_current, x_current_interp);
            NRMSE2 = NRMSE(x_current, x_rec);
            fprintf(['\n Reconstruction in ',num2str(toc(start),'%.3f'),' seconds, with improvement of NRMSE is ', num2str(100*(1-NRMSE2/NRMSE1),'%.3f'), ' \n'])
        end
    end
end 
close(nc1); 

%%
% imagesc(weights1)   