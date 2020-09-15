%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 Boualem Boashash
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%
% Authors: Prof. Boualem Boashash        (boualem.boashash@gmail.com)
%          Dr. Abdeldjalil Aissa-El-Bey  (abdeldjalil.aissaelbey@telecom-bretagne.eu)
%          RA: Md.F.A
%
% The following references should be cited whenever this script is used:
% [1] B. Boashash, A. Aissa-El-Bey, Multisensor Time-Frequency Signal Processing:
%     A tutorial review with illustrations in selected application areas, Digital
%     Signal Processing, In Press.
% [2] B. Boashash, A. Aissa-El-Bey, M. F. Al-Sa'd, Multisensor time-frequency
%     signal processing software Matlab package: An analysis tool for multichannel
%     non-stationary data , SoftwareX, In Press.
%
% Last Modification: 25-04-2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 DOA estimation using TF MUSIC algorithm
%
% Syntax : P = tf_music(Ds, n, m, lamda, d, theta)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <INPUTs>
% Ds       : Selected MTFD Matrix (m x m).
% n        : Number of sources.
% m        : Number of sensors.
% lamda    : Wavelength of the received signal.
% d        : ULA elements spacing.
% theta    : Solution space for the 'MUSIC'.
%
% <OUTPUTs>
% P        : Estimated Spectrum for 'MUSIC'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function P = tf_music(Ds, n, m, lamda, d, theta)
theta   = theta/180*pi;
theta_N = length(theta);
[vv,~]  = svd(Ds);        % Find the eigenvalues and eigenvectors of STFD
UN_TF   = vv(:,n+1:m);    % Estimate/Selection of noise subspace

%% TF MUSIC Main
P = zeros(1,theta_N);
for ii = 1:theta_N
    a_theta = exp(-1i*2*pi*(d/lamda)*sin(theta(ii))*(0:m-1));

    PP_TF   = conj(a_theta)*(UN_TF*UN_TF')*a_theta.';
    P(ii)   = abs(1/PP_TF);
end
P = P/max(P);
end