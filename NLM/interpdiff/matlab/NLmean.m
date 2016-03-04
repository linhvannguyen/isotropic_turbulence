function [x_rec,weights] = NLmean(x_ref, x_moving, x_fusion, gridy, gridz, simpatch_haftsize, accpatch_haftsize, neighbor_haftsize, tau)
%NLmean: non-local similarity-based reconstruction: x_fusion is
    % Prepare some constants
    [Ny,Nz] = size(x_ref);
    dim_neighbor = (2*neighbor_haftsize+1)^2;
    dim_simpatch = (2*simpatch_haftsize+1)^2;
    
    acc_ids_1D = simpatch_haftsize+1-accpatch_haftsize:simpatch_haftsize+1+accpatch_haftsize;
    acc_ids = reshape(1:dim_simpatch,[2*simpatch_haftsize+1,2*simpatch_haftsize+1]);
    acc_ids = acc_ids(acc_ids_1D, acc_ids_1D);
    acc_ids = acc_ids(:);
 
    % Load reference patches of size (n,n,w1,w1). 
    ref_patches = extract_patch(x_ref, squeeze(gridy(:,:,ceil(dim_neighbor/2),:,:)), squeeze(gridz(:,:,ceil(dim_neighbor/2),:,:)));
    ref_patches = reshape(ref_patches,[Ny Nz 1 dim_simpatch]);
    ref_patches = repmat(ref_patches, [1 1 dim_neighbor  1]);

    % Load moving patches of size (n,n,m,m,w1,w1). 
    moving_patches = extract_patch(x_moving, gridy, gridz); 

    % Compute distance
    D = (ref_patches - moving_patches).^2;
    D = mean(D,4); 
    K=exp(-D/(2.0*tau*tau));

    % fusion of information
    x_rec=zeros(Ny,Nz); 
    weights=zeros(Ny,Nz);
    

%     patch = @(i,j) extract_patch(x_fusion, squeeze(gridy(i,j,:,acc_ids)), squeeze(gridz(i,j,:,acc_ids)));
%     row = @(i,j) squeeze(gridy(i,j,ceil(dim_neighbor/2),acc_ids));
%     col = @(i,j) squeeze(gridz(i,j,ceil(dim_neighbor/2),acc_ids));
%     NLM = @(i,j) put_patch(x_rec, weights, patch(i,j), row(i,j), col(i,j),squeeze(K(i,j,:)));
%     [Z,Y] = meshgrid(1:Nz,1:Ny);
%     [x_rec, weights] = arrayfun(@(i,j)NLM(i,j), Y,Z, 'UniformOutput', false);

    for i=1:Ny
        for j=1:Nz
            patches=extract_patch(x_fusion, squeeze(gridy(i,j,:,acc_ids)), squeeze(gridz(i,j,:,acc_ids)));
            w = squeeze(K(i,j,:));
            rows_ref = squeeze(gridy(i,j,ceil(dim_neighbor/2),acc_ids));
            cols_ref = squeeze(gridz(i,j,ceil(dim_neighbor/2),acc_ids));
            [x_rec, weights]=put_patch(x_rec, weights, patches, rows_ref, cols_ref,w);
        end
    end 
end


function [x_2D, W_2D]=put_patch(x_2D, W_2D, patches, rows_ref, cols_ref, w)
%put_patch accumulate moving patches of size [m1xm2xn1xn2] into a 2D image at
% the positions rows_ref [n1xn2] and cols_ref[n1xn2]
    patches=reshape(patches,[numel(w) numel(rows_ref)]);
    temp = patches.*repmat(w,[1 numel(rows_ref)]);
    x_2D(rows_ref+(cols_ref-1)*size(x_2D,1)) = x_2D(rows_ref+(cols_ref-1)*size(x_2D,1)) + sum(temp,1)';
    W_2D(rows_ref+(cols_ref-1)*size(x_2D,1)) = W_2D(rows_ref+(cols_ref-1)*size(x_2D,1)) + sum(w(:));
end


function patches=extract_patch(x_2D, rows, cols)
    patches=x_2D(rows + (cols-1)*size(x_2D,1));
end