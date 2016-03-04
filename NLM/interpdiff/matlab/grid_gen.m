function [gridy, gridz] = grid_gen(Ny,Nz,simpatch_haftsize,neighbor_haftsize)
%grid_gen generate grids of size [Ny,Nz,dim_neighbor,dim_simpatch]
    dim_neighbor = (2*neighbor_haftsize+1)^2;
    dim_simpatch = (2*simpatch_haftsize+1)^2;

    % grids of reference pixels
    [gridz,gridy] = meshgrid(1:Nz,1:Ny);

    % neighboring offsets
    [dz,dy] = meshgrid(-neighbor_haftsize:neighbor_haftsize,-neighbor_haftsize:neighbor_haftsize);
    dz = reshape(dz(:), [1 1 dim_neighbor]);
    dy = reshape(dy(:), [1 1 dim_neighbor]);
    gridz = repmat(gridz, [1 1 dim_neighbor]) + repmat(dz, [Ny Nz 1]);
    gridy = repmat(gridy, [1 1 dim_neighbor]) + repmat(dy, [Ny Nz 1]);

    % similarity patch offsets
    [dz,dy] = meshgrid(-simpatch_haftsize:simpatch_haftsize,-simpatch_haftsize:simpatch_haftsize);
    dz = reshape(dz(:), [1 1 1 dim_simpatch]);
    dy = reshape(dy(:), [1 1 1 dim_simpatch]);
    gridz = repmat(gridz, [1 1 1 dim_simpatch]) + repmat(dz, [Ny Nz dim_neighbor 1]);
    gridy = repmat(gridy, [1 1 1 dim_simpatch]) + repmat(dy, [Ny Nz dim_neighbor 1]);

    % handle boundary condition by periodicity
    gridy = border_per(gridy,1,Ny); 
    gridz = border_per(gridz,1,Nz); 
end

function ids = border_per (ids, minid, maxid) 
%border_per handle boundary by periodicity
    ids(ids<minid) = maxid + ids(ids<minid);
    ids(ids>maxid) = ids(ids>maxid) - maxid;
end
