function analytic_sig_ker = signal_kernal(x)
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


%
%    Both real and analytic versions.
%
%   K(n,m) = z(n+m)z*(n-m)
%
% where, z() is the analytic associate of real input signal s()
%    n is the sampled or discrete time
%    m is the sample or disrete lag
%



N = length(x);

if isreal(x)
    if mod(length(x),2) == 0
        true_X = fft(x);
        analytic_X = [true_X(1) 2.*true_X(2:N/2) true_X(N/2+1) zeros(1,N/2-1)];
        analytic_x = ifft(analytic_X);
    else
        true_X = fft(x);
        analytic_X = [true_X(1) 2.*true_X(2:ceil(N/2)) zeros(1,floor(N/2))];
        analytic_x = ifft(analytic_X);
    end
else
    analytic_x=x;
end

analytic_sig_ker = zeros(N,N);
for m = -round(N/2-1):1:round(N/2-1);
    analytic_sig_ker(m+round(N/2)+1,:) = sig_ker_corr(analytic_x,m);
end


