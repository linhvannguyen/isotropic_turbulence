function x_interp = interp_border (x_ref, space_spacing)
%interp_border interpolate, taking care of border by periodicity
[Ny,Nz] = size(x_ref);
HTLS_idy= 1:space_spacing:Ny; % row indices of HTLS in space
HTLS_idz= 1:space_spacing:Nz; % column indices of HTLS in space

left=4*space_spacing; right = 4*space_spacing-(Ny-HTLS_idy(end));
bottom=4*space_spacing; top = 4*space_spacing-(Nz-HTLS_idz(end));

[gridz_HS_enlarged, gridy_HS_enlarged] = meshgrid(1:Nz+left+right,1:Ny+bottom+top);
gridz_LS_enlarged=gridz_HS_enlarged(1:space_spacing:end,1:space_spacing:end);
gridy_LS_enlarged=gridy_HS_enlarged(1:space_spacing:end,1:space_spacing:end);

x_ref_enlarged=enlarge_2D(x_ref,left, right, bottom, top);
x_LR_enlarged = x_ref_enlarged(1:space_spacing:end,1:space_spacing:end);
x_interp_enlarged=interp2(gridz_LS_enlarged, gridy_LS_enlarged, x_LR_enlarged, gridz_HS_enlarged, gridy_HS_enlarged,'spline');
x_interp=x_interp_enlarged(bottom+1:bottom+Ny,left+1:left+Nz);

