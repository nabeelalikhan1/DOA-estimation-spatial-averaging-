function PQ = paddedsize(AB, CD, PARAM)
%PADDEDSIZE Computes padded sizes useful for FFT-based filtering.
% PQ = PADDEDSIZE(AB), where AB is a two-element size vector,
% computes the two-element size vector PQ = 2*AB.
%
% PQ = PADDEDSIZE(AB, 'PWR2') computes the vector PQ such that
% PQ(1) = PQ(2) = 2^nextpow2(2*m), where m is MAX(AB).
%
% PQ = PADDEDSIZE(AB, CD), where AB and CD are two-element size
% vectors, computes the two-element size vector PQ. The elements
% of PQ are the smallest even integers greater than or equal to
% AB + CD -1.
%
% PQ = PADDEDSIZE(AB, CD, 'PWR2') computes the vector PQ such that
% PQ(1) = PQ(2) = 2^nextpow2(2*m), where m is MAX([AB CD]).
if nargin == 1
 PQ = 2*AB;
elseif nargin == 2 & ~ischar(CD)
 PQ = AB + CD - 1;
 PQ = 2 * ceil(PQ / 2);
elseif nargin == 2
 m = max(AB); % Maximum dimension.
 % Find power-of-2 at least twice m.
 P = 2^nextpow2(2*m);
 PQ = [P, P];
elseif nargin == 3
 m = max([AB CD]); %Maximum dimension.
 P = 2^nextpow2(2*m);
 PQ = [P, P];
else
 error('Wrong number of inputs.')
end
function [U, V] = dftuv(M, N)
%DFTUV Computes meshgrid frequency matrices.
% [U, V] = DFTUV(M, N) computes meshgrid frequency matrices U and
% V. U and V are useful for computing frequency-domain filter
% functions that can be used with DFTFILT. U and V are both M-by-N.
% Set up range of variables.
u = 0:(M-1);
v = 0:(N-1);
% Compute the indices for use in meshgrid
idx = find(u > M/2);
u(idx) = u(idx) - M;
idy = find(v > N/2);
v(idy) = v(idy) - N;
% Compute the meshgrid arrays
[V, U] = meshgrid(v, u);
function H = lpfilter(type, M, N, D0, n)
%LPFILTER Computes frequency domain lowpass filters
% H = LPFILTER(TYPE, M, N, D0, n) creates the transfer function of
% a lowpass filter, H, of the specified TYPE and size (M-by-N). To
% view the filter as an image or mesh plot, it should be centered
% using H = fftshift(H).
%
% Valid values for TYPE, D0, and n are:
%
% 'ideal' Ideal lowpass filter with cutoff frequency D0. n need
% not be supplied. D0 must be positive
%
% 'btw' Butterworth lowpass filter of order n, and cutoff D0.
% The default value for n is 1.0. D0 must be positive.
%
% 'gaussian' Gaussian lowpass filter with cutoff (standard deviation)
% D0. n need not be supplied. D0 must be positive.
% Use function dftuv to set up the meshgrid arrays needed for
% computing the required distances.
[U, V] = dftuv(M, N);
% Compute the distances D(U, V).
D = sqrt(U.^2 + V.^2);
% Begin fiter computations.
switch type
case 'ideal'
 H = double(D <=D0);
case 'btw'
 if nargin == 4
 n = 1;
 end
 H = 1./(1 + (D./D0).^(2*n));
case 'gaussian'
 H = exp(-(D.^2)./(2*(D0^2)));
otherwise
 error('Unknown filter type.')
end
function H = hpfilter(type, M, N, D0, n)
%HPFILTER Computes frequency domain highpass filters
% H = HPFILTER(TYPE, M, N, D0, n) creates the transfer function of
% a highpass filter, H, of the specified TYPE and size (M-by-N).
% Valid values for TYPE, D0, and n are:
%
% 'ideal' Ideal highpass filter with cutoff frequency D0. n
% need not be supplied. D0 must be positive
%
% 'btw' Butterworth highpass filter of order n, and cutoff D0.
% The default value for n is 1.0. D0 must be positive.
%
% 'gaussian' Gaussian highpass filter with cutoff (standard deviation)
% D0. n need not be supplied. D0 must be positive.
%
% The transfer function Hhp of a highpass filter is 1 - Hlp,
% where Hlp is the transfer function of the corresponding lowpass
% filter. Thus, we can use function lpfilter to generate highpass
% filters.
if nargin == 4
 n = 1; % Default value of n.
end
% Generate highpass filter.
Hlp = lpfilter(type, M, N, D0, n);
H = 1 - Hlp;
function g = dftfilt(f, H)
%DFTFILT Performs frequency domain filtering.
% G = DFTFILT(F, H) filters F in the frequency domain using the
% filter transfer function H. The output, G, is the filtered
% image, which has the same size as F. DFTFILT automatically pads
% F to be the same size as H. Function PADDEDSIZE can be used to
% determine an appropriate size for H.
%
% DFTFILT assumes that F is real and that H is a real, uncentered
% circularly-symmetric filter function.
% Copyright 2002-2004 R. C. Gonzalez, R. E. Woods, & S. L. Eddins
% Digital Image Processing Using MATLAB, Prentice-Hall, 2004
% $Revision: 1.5 $ $Date: 2003/08/25 14:28:22 $
% Obtain the FFT of the padded input.
F = fft2(f, size(H, 1), size(H, 2));
% Perform filtering.
g = real(ifft2(H.*F));
% Crop to original size.
g = g(1:size(f, 1), 1:size(f, 2));