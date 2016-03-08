%%
load wind %load built-in vector field dataset 
[cu,cv,cw] = curl(x,y,z,u,v,w); %calculate vorticity of the vector field 
div = divergence(x,y,z,u,v,w); %calculate divergence of the vector field 
vtkwrite('wind.vtk', 'structured_grid',x,y,z, 'vectors','vector_field',u,v,w, 'vectors','vorticity',cu,cv,cw, 'scalars','divergence',div); 

%%
ncload('/data/ISOTROPIC/data/FIELD-020.nc');
Nh=size(velocity_x,1);
[gridz, gridy, gridx] = meshgrid(0:1/(Nh-1):1, 0:1/(Nh-1):1, 0:1/(Nh-1):1);

[cz,cy,cx] = curl(gridz,gridy,gridx,velocity_z, velocity_y, velocity_x); %calculate vorticity of the vector field 

vtkwrite('/data/ISOTROPIC/data/field.vtk', 'structured_grid', gridz, gridy, gridx, 'vectors','vector_field', velocity_z, velocity_y, velocity_x, 'vectors','vorticity',cz,cy,cx); 