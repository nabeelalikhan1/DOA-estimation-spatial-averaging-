function P = music(x, n, m, lamda, d, theta)
% Syntax : P = music(x, n, m, lamda, d, theta)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <INPUTs>
% x        : Input signals (m x N).
% n        : Number of sources.
% m        : Number of sensors.
% lamda    : Wavelength of the received signal.
% d        : ULA elements spacing.
% theta    : Solution space for the 'MUSIC'.
%
% <OUTPUTs>
% P        : Estimated Spectrum for 'MUSIC'.
theta   = theta/180*pi;
N       = length(x); 
theta_N = length(theta);
Rxx     = (1/N)*(x*x');   % Covariance Matrix
[vv,~]  = svd(Rxx);       % Find the eigenvalues and eigenvectors of Rxx
NN      = vv(:,n+1:m);    % Estimate/Selection of noise subspace

%% MUSIC Main
P = zeros(1,theta_N);
for ii = 1:theta_N
    a_theta = exp(-1j*2*pi*(d/lamda)*cos(theta(ii))*(0:m-1));
    P_temp  = conj(a_theta)*(NN*NN')*a_theta.';
    P(ii) = abs(1/P_temp);
end
P = P/max(P);
end