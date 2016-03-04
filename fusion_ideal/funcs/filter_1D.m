function x_1D_fil = filter_1D(x_1D, scale)

N = numel(x_1D);
ks = (mod(1/2 + (0:(N-1))/N, 1) - 1/2); 
kc = max(ks)/scale;
LPF = (ks.^2 < kc^2); 

F = fft(x_1D);
x_1D_fil=ifft(LPF.*F);

end