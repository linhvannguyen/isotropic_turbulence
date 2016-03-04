function [x_2D, W_2D]=put_patch(x_2D, W_2D, patches, rows_ref, cols_ref, w)
%put_patch accumulate moving patches of size [m1xm2xn1xn2] into a 2D image at
% the positions rows_ref [n1xn2] and cols_ref[n1xn2]
patches=reshape(patches,[size(w) size(rows_ref)]);
temp = patches.*repmat(w,[1 1 size(rows_ref)]);
x_2D(rows_ref, cols_ref) = x_2D(rows_ref, cols_ref) + sum(temp(:));
W_2D(rows_ref, cols_ref) = W_2D(rows_ref, cols_ref) + sum(w(:));