function out = ifft_1D(in,dim)
% One-dimensional inverse Fourier transformation
% of the input array in along the dimension dim

out = ifftshift(ifft(ifftshift(in,dim),[],dim),dim);

end