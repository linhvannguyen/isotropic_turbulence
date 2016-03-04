function x_2D_expand = enlarge_2D(x_2D, left, right, bottom, top)
%EXPAND_2D expand a 2D field by [left,right,bottom,top]

[nrow,ncol]=size(x_2D);

[gridy,gridx] = meshgrid(-left+1:ncol+right,-bottom+1:nrow+top);

gridx(gridx<1) = ncol+gridx(gridx<1); 
gridy(gridy<1) = nrow+gridy(gridy<1);
gridx(gridx>nrow) = gridx(gridx>nrow) - nrow; 
gridy(gridy>ncol) = gridy(gridy>ncol) - ncol;

x_2D_expand = x_2D(gridx + (gridy-1)*nrow);
end