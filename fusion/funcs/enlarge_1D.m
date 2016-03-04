function x_1D_expand = enlarge_1D(x_1D, left, right)
%EXPAND_1D expand a 1D signal by [left,right]

N=numel(x_1D);

gridx = -left+1:N+right;

gridx(gridx<1) = N+gridx(gridx<1); 
gridx(gridx>N) = gridx(gridx>N) - N; 

x_1D_expand = x_1D(gridx);

end