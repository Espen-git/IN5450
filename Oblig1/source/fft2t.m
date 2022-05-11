 function [out] = fft2t(in, NFFT)
%
%  Function for doing a 2-d fft and shifting everything 
%  appropriately.
%
%  Based on fft2s by Hayden John Callow

if ~exist('NFFT','var') || all(size(in)==NFFT)
   out = fftshift(fft2(ifftshift(in))); 
else
   [m, n] = size(in);
   tmp = complex( zeros(NFFT, class(in)));
%   tmp(1:m,1:n) = in;
   startind = floor((NFFT - m)/2)
   tmp(startind+(1:m),startind+(1:n)) = in;
   out = fftshift(fft2(ifftshift(tmp))); 
end
