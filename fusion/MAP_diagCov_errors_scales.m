clear all; close all; clc;

filename_ref='/data/ISOTROPIC/data/data_downsampled4.nc';

%%  TIME: 04; SPACE: 03
space_spacing=03; % subsampling ration in space
time_spacing=04; % subsampling ration in time (from 40Hz to 4Hz)

fprintf('compute errors TIME: 04  SPACE: 03\n');

filename_interp_space=strcat('/data/ISOTROPIC/fusion/Uinterp_spatialspacing_',num2str(space_spacing,'%.2d'),'.nc');
filename_interp_time=strcat('/data/ISOTROPIC/fusion/Uinterp_timespacing_',num2str(time_spacing,'%.2d'),'.nc');
filename_fusion=strcat('/data/ISOTROPIC/fusion/FusedData_diagCov_timespacing_',num2str(time_spacing,'%.2d'),'_spacespacing_',num2str(space_spacing,'%.2d'),'.nc');

nc1=netcdf(filename_ref,'r');
nc2 = netcdf(filename_interp_space,'r');
nc3=netcdf(filename_interp_time,'r');
nc4=netcdf(filename_fusion,'r');

x_nonoise=nc1{'velocity_x'}(:,:,:,:);
x_interp_spatial=nc2{'Uinterp'}(:,:,:,:);
x_interp_temporal=nc3{'Uinterp'}(:,:,:,:);
x_MAP=nc4{'Zhat_all'}(:,:,:,:);

close(nc1); close(nc2); close(nc3); close(nc4);

err_interp_spatial_time4_space3=sqrt(sum((x_nonoise(:)-x_interp_spatial(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
err_interp_temporal_time4_space3=sqrt(sum((x_nonoise(:)-x_interp_temporal(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
err_fusion_time4_space3=sqrt(sum((x_nonoise(:)-x_MAP(:)).^2))/sqrt(sum((x_nonoise(:)).^2));

%%  TIME: 06; SPACE: 04
space_spacing=4; % subsampling ration in space
time_spacing=6; % subsampling ration in time (from 40Hz to 4Hz)

fprintf('compute errors TIME: 06  SPACE: 04\n');

filename_interp_space=strcat('/data/ISOTROPIC/fusion/Uinterp_spatialspacing_',num2str(space_spacing,'%.2d'),'.nc');
filename_interp_time=strcat('/data/ISOTROPIC/fusion/Uinterp_timespacing_',num2str(time_spacing,'%.2d'),'.nc');
filename_fusion=strcat('/data/ISOTROPIC/fusion/FusedData_diagCov_timespacing_',num2str(time_spacing,'%.2d'),'_spacespacing_',num2str(space_spacing,'%.2d'),'.nc');

nc1=netcdf(filename_ref,'r');
nc2 = netcdf(filename_interp_space,'r');
nc3=netcdf(filename_interp_time,'r');
nc4=netcdf(filename_fusion,'r');

x_nonoise=nc1{'velocity_x'}(:,:,:,:);
x_interp_spatial=nc2{'Uinterp'}(:,:,:,:);
x_interp_temporal=nc3{'Uinterp'}(:,:,:,:);
x_MAP=nc4{'Zhat_all'}(:,:,:,:);

close(nc1); close(nc2); close(nc3); close(nc4);

err_interp_spatial_time6_space4=sqrt(sum((x_nonoise(:)-x_interp_spatial(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
err_interp_temporal_time6_space4=sqrt(sum((x_nonoise(:)-x_interp_temporal(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
err_fusion_time6_space4=sqrt(sum((x_nonoise(:)-x_MAP(:)).^2))/sqrt(sum((x_nonoise(:)).^2));


%%  TIME: 08; SPACE: 06
space_spacing=6; % subsampling ration in space
time_spacing=8; % subsampling ration in time (from 40Hz to 4Hz)

fprintf('compute errors TIME: 08  SPACE: 06\n');

filename_interp_space=strcat('/data/ISOTROPIC/fusion/Uinterp_spatialspacing_',num2str(space_spacing,'%.2d'),'.nc');
filename_interp_time=strcat('/data/ISOTROPIC/fusion/Uinterp_timespacing_',num2str(time_spacing,'%.2d'),'.nc');
filename_fusion=strcat('/data/ISOTROPIC/fusion/FusedData_diagCov_timespacing_',num2str(time_spacing,'%.2d'),'_spacespacing_',num2str(space_spacing,'%.2d'),'.nc');

nc1=netcdf(filename_ref,'r');
nc2 = netcdf(filename_interp_space,'r');
nc3=netcdf(filename_interp_time,'r');
nc4=netcdf(filename_fusion,'r');

x_nonoise=nc1{'velocity_x'}(:,:,:,:);
x_interp_spatial=nc2{'Uinterp'}(:,:,:,:);
x_interp_temporal=nc3{'Uinterp'}(:,:,:,:);
x_MAP=nc4{'Zhat_all'}(:,:,:,:);

close(nc1); close(nc2); close(nc3); close(nc4);

err_interp_spatial_time8_space6=sqrt(sum((x_nonoise(:)-x_interp_spatial(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
err_interp_temporal_time8_space6=sqrt(sum((x_nonoise(:)-x_interp_temporal(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
err_fusion_time8_space6=sqrt(sum((x_nonoise(:)-x_MAP(:)).^2))/sqrt(sum((x_nonoise(:)).^2));


%%  TIME: 08; SPACE: 03
space_spacing=03; % subsampling ration in space
time_spacing=08; % subsampling ration in time (from 40Hz to 4Hz)

fprintf('compute errors TIME: 08  SPACE: 03\n');

filename_interp_space=strcat('/data/ISOTROPIC/fusion/Uinterp_spatialspacing_',num2str(space_spacing,'%.2d'),'.nc');
filename_interp_time=strcat('/data/ISOTROPIC/fusion/Uinterp_timespacing_',num2str(time_spacing,'%.2d'),'.nc');
filename_fusion=strcat('/data/ISOTROPIC/fusion/FusedData_diagCov_timespacing_',num2str(time_spacing,'%.2d'),'_spacespacing_',num2str(space_spacing,'%.2d'),'.nc');

nc1=netcdf(filename_ref,'r');
nc2 = netcdf(filename_interp_space,'r');
nc3=netcdf(filename_interp_time,'r');
nc4=netcdf(filename_fusion,'r');

x_nonoise=nc1{'velocity_x'}(:,:,:,:);
x_interp_spatial=nc2{'Uinterp'}(:,:,:,:);
x_interp_temporal=nc3{'Uinterp'}(:,:,:,:);
x_MAP=nc4{'Zhat_all'}(:,:,:,:);

close(nc1); close(nc2); close(nc3); close(nc4);

err_interp_spatial_time8_space3=sqrt(sum((x_nonoise(:)-x_interp_spatial(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
err_interp_temporal_time8_space3=sqrt(sum((x_nonoise(:)-x_interp_temporal(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
err_fusion_time8_space3=sqrt(sum((x_nonoise(:)-x_MAP(:)).^2))/sqrt(sum((x_nonoise(:)).^2));


%%  TIME: 08; SPACE: 04
space_spacing=04; % subsampling ration in space
time_spacing=08; % subsampling ration in time (from 40Hz to 4Hz)

fprintf('compute errors TIME: 08  SPACE: 04\n');

filename_interp_space=strcat('/data/ISOTROPIC/fusion/Uinterp_spatialspacing_',num2str(space_spacing,'%.2d'),'.nc');
filename_interp_time=strcat('/data/ISOTROPIC/fusion/Uinterp_timespacing_',num2str(time_spacing,'%.2d'),'.nc');
filename_fusion=strcat('/data/ISOTROPIC/fusion/FusedData_diagCov_timespacing_',num2str(time_spacing,'%.2d'),'_spacespacing_',num2str(space_spacing,'%.2d'),'.nc');

nc1=netcdf(filename_ref,'r');
nc2 = netcdf(filename_interp_space,'r');
nc3=netcdf(filename_interp_time,'r');
nc4=netcdf(filename_fusion,'r');

x_nonoise=nc1{'velocity_x'}(:,:,:,:);
x_interp_spatial=nc2{'Uinterp'}(:,:,:,:);
x_interp_temporal=nc3{'Uinterp'}(:,:,:,:);
x_MAP=nc4{'Zhat_all'}(:,:,:,:);

close(nc1); close(nc2); close(nc3); close(nc4);

err_interp_spatial_time8_space4=sqrt(sum((x_nonoise(:)-x_interp_spatial(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
err_interp_temporal_time8_space4=sqrt(sum((x_nonoise(:)-x_interp_temporal(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
err_fusion_time8_space4=sqrt(sum((x_nonoise(:)-x_MAP(:)).^2))/sqrt(sum((x_nonoise(:)).^2));


%%  TIME: 04; SPACE: 06
space_spacing=6; % subsampling ration in space
time_spacing=4; % subsampling ration in time (from 40Hz to 4Hz)

fprintf('compute errors TIME: 04  SPACE: 06\n');

filename_interp_space=strcat('/data/ISOTROPIC/fusion/Uinterp_spatialspacing_',num2str(space_spacing,'%.2d'),'.nc');
filename_interp_time=strcat('/data/ISOTROPIC/fusion/Uinterp_timespacing_',num2str(time_spacing,'%.2d'),'.nc');
filename_fusion=strcat('/data/ISOTROPIC/fusion/FusedData_diagCov_timespacing_',num2str(time_spacing,'%.2d'),'_spacespacing_',num2str(space_spacing,'%.2d'),'.nc');

nc1=netcdf(filename_ref,'r');
nc2 = netcdf(filename_interp_space,'r');
nc3=netcdf(filename_interp_time,'r');
nc4=netcdf(filename_fusion,'r');

x_nonoise=nc1{'velocity_x'}(:,:,:,:);
x_interp_spatial=nc2{'Uinterp'}(:,:,:,:);
x_interp_temporal=nc3{'Uinterp'}(:,:,:,:);
x_MAP=nc4{'Zhat_all'}(:,:,:,:);

close(nc1); close(nc2); close(nc3); close(nc4);

err_interp_spatial_time4_space6=sqrt(sum((x_nonoise(:)-x_interp_spatial(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
err_interp_temporal_time4_space6=sqrt(sum((x_nonoise(:)-x_interp_temporal(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
err_fusion_time4_space6=sqrt(sum((x_nonoise(:)-x_MAP(:)).^2))/sqrt(sum((x_nonoise(:)).^2));


%%  TIME: 06; SPACE: 06
space_spacing=6; % subsampling ration in space
time_spacing=6; % subsampling ration in time (from 40Hz to 4Hz)

fprintf('compute errors TIME: 06  SPACE: 06\n');

filename_interp_space=strcat('/data/ISOTROPIC/fusion/Uinterp_spatialspacing_',num2str(space_spacing,'%.2d'),'.nc');
filename_interp_time=strcat('/data/ISOTROPIC/fusion/Uinterp_timespacing_',num2str(time_spacing,'%.2d'),'.nc');
filename_fusion=strcat('/data/ISOTROPIC/fusion/FusedData_diagCov_timespacing_',num2str(time_spacing,'%.2d'),'_spacespacing_',num2str(space_spacing,'%.2d'),'.nc');

nc1=netcdf(filename_ref,'r');
nc2 = netcdf(filename_interp_space,'r');
nc3=netcdf(filename_interp_time,'r');
nc4=netcdf(filename_fusion,'r');

x_nonoise=nc1{'velocity_x'}(:,:,:,:);
x_interp_spatial=nc2{'Uinterp'}(:,:,:,:);
x_interp_temporal=nc3{'Uinterp'}(:,:,:,:);
x_MAP=nc4{'Zhat_all'}(:,:,:,:);

close(nc1); close(nc2); close(nc3); close(nc4);

err_interp_spatial_time6_space6=sqrt(sum((x_nonoise(:)-x_interp_spatial(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
err_interp_temporal_time6_space6=sqrt(sum((x_nonoise(:)-x_interp_temporal(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
err_fusion_time6_space6=sqrt(sum((x_nonoise(:)-x_MAP(:)).^2))/sqrt(sum((x_nonoise(:)).^2));