%% LOAD TO NETCDF
clear all; close all; clc;
Nh=96; % Number of grid-point in HR field

%% PLANES TO LEARN THE DICTIONARIES
nc = netcdf('/data/ISOTROPIC/data/data_downsampled4.nc','r');
U_HR_3D=nc{'velocity_x'}(1,:,:,:);
close(nc)
%%
ids=2:3:Nh;
[gridz, gridy, gridx] = meshgrid(0:1/(Nh-1):1, 0:1/(Nh-1):1, 0:1/(Nh-1):1);

%% 3D isosurface
fig1=figure();
axis equal;
box on;
set(fig1,'color','w')

fv = isosurface (gridz,gridy,gridx,U_HR_3D,0);
p = patch(fv);
set(p,'FaceColor','red','EdgeColor','none');
camlight;
lighting gouraud;
xlabel('x');
ylabel('y');
zlabel('z');

view(45,45)
camlight(-45,45)

xlim([0,1]);
ylim([0,1]);
zlim([0,1]);

% print('./figures/isosurf_u_HR_3D','-dpng')

%% PLOT 2D plane
u_HR_2D=squeeze(U_HR_3D(:,:,1));

fig2=figure(); set(fig2,'color','w'); imagesc(u_HR_2D); caxis([-4,4]); axis off; axis equal; 
export_fig('./figures/u_HR_2D','-eps','-q101','-a4'); %close();
