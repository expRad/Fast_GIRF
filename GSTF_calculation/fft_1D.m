function out = fft_1D(in,dim)
% One-dimensional Fourier transformation
% of the input array in along the dimension dim

out = fftshift(fft(fftshift(in,dim),[],dim),dim);

end