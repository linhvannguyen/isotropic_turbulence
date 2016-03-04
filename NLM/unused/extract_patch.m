function patches=extract_patch(x_2D, rows, cols)
    patches=x_2D(rows + (cols-1)*size(x_2D,1));
