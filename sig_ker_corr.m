function sig_ker_m = sig_ker_corr(x,m)
%  Author:
%     Samir Ouelha (post-doc of Prof. Boualem Boashash)
%

%Please cite the following reference:

%   [1] B.Boashash and S. Ouelha
%       “An Improved design of high-resolution Quadratic Time-Frequency
%       distributions for the analysis of non-stationary multicomponent 
%       signals using directional compact kernels”
%       , IEEE Transactions on Signal Processing, 2016
%
% This work was supported by ARC and QNRF grants.
% A correlation function to aid in the creation
% of the signal kernel
%

N = length(x);
z_nmm = [zeros(1,N+m) x zeros(1,N-m)];
z_npm = [zeros(1,N-m) x zeros(1,N+m)];
ker_nm = z_npm.*conj(z_nmm);
sig_ker_m = ker_nm(N+1:2*N);
