function [amb, tfrep] = wvd1(x)
%

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

% Wigner Ville Distribution
%  No windowing or time-resolution variablility

N = length(x);
analytic_sig_ker = signal_kernal(x);
tfrep = real(1./N.*fft(ifftshift(analytic_sig_ker,1), N, 1));
amb = fftshift(1./N.*fft(analytic_sig_ker, N, 2),2);
