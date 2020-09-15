function Pmusic = TMMUSIC(Rxx, lambda, m, n, d, theta)

[Ntime,Vtime]  = svd(Rxx);              % Find the eigenvalues and eigenvectors of Rxx
NN             = Ntime(:,n+1:m);        % Estimate/Selection of noise subspace
Pmusic         = zeros(1,length(theta));
for ii=1:length(theta)
    a_theta=zeros(1,m);
    for jj  =  0:m-1
        a_theta(1+jj)=exp(-1i*2*jj*pi*d*sin(theta(ii)/180*pi)/lambda);
    end

    PP      =  a_theta*(NN*NN')*a_theta';

    Pmusic(ii)  =  abs(1/ PP);
end
Pmusic = Pmusic/max(Pmusic);

%Pmusic   =    10*log10(Pmusic/max(Pmusic)); %Spatial spectrum function
