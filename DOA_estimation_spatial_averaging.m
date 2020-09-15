close all;
clear all;
N_sensors=10;
n=0:127;

%addpath('D:\D\win64_bin\win64_bin');
addpath('D:\tfsa_5-5\windows\win64_bin');
%crossing componentsi8

s1=exp(2*pi*1i*(0.05*n+0.45*n.^3/(128*128*3)));
s3=exp(2*pi*1i*(0.11*n+0.45*n.^3/(128*128*3)));%+exp(2*pi*1i*(0.45*n+0*0.1*n.^2/(2*128)-0*0.5*n.^3/(128*128*3)));
%s3=exp(2*pi*1i*(0.45*n-0.45*n.^3/(128*128*3)));
%crossing components LFM sources
%s2=exp(2*pi*1i*(0.1*n+0.2*n.^2/(2*128)));
%s3=exp(2*pi*1i*(0.3*n-0.2*n.^2/(2*128)));

s = [(s1.')  (s3.') ];

% s1=exp(2*pi*1i*(0.05*n+0*0.2*n.^2/(2*128)));

n_sources=2;
N_C=2;
s_orig=s;

% set mixing matrix A
%theta = [15,30,50]*pi/180;   % sensor separation angles in radians
theta = [-4,4]*pi/180;   % sensor separation angles in radians

A = exp(1j*pi*[0:N_sensors-1].'*sin(theta));  % mixing matrix A
theta9=round(theta *180/pi);
% generate noise
SNR=-10;
%SNR=0;

N_sim=200;
N_sim=1000;

for iii=1:N_sim
X = A*s.';                             % mixed source

sigma = 10^(-SNR/20);
w = sigma*(randn(N_sensors,length(n)) + 1j*(randn(N_sensors,length(n))))/sqrt(2); % noise

X=X+w;

%D = mtfd(X, 'CKD',0.05,0.05);
D = mtfd(X, 'CKD',1,0.25,0.25);

I_avg=zeros(length(s),length(s));
for ii=1:N_sensors
    I=D{ii,ii};
    I_avg=I+I_avg;
end

I_avg=real(I_avg)/N_sensors;
D_avg=I_avg;
for ii=1:N_sensors
    DD{ii,ii}=I_avg;
end

I_avg1=I_avg;
I_avg1(I_avg1<0)=0;
%figure;imagesc(real(I_avg1));
lag=N_sensors-1;
DD=D;
for jj=1:lag
    I_avg=zeros(length(s),length(s));
    
    for ii=1:N_sensors-jj
        I1=D{ii,ii+jj};
        I_avg=I_avg+I1;
    end
    I_avg=I_avg/(N_sensors-jj);
    for ii=1:N_sensors-jj
        DD{ii,ii+jj}=I_avg;
        DD{ii+jj,ii}=conj(I_avg);
        
    end
end
%%% DOA Estimation TF MUSIC

perc=0.4;
thr = perc*max(max(D_avg));
Tr = (D_avg) >= thr;
[F_trace, ~] = find(Tr);
n_p = length(F_trace);
D_s = zeros(N_sensors, N_sensors);
for m1 = 1:N_sensors
    for m2 = 1:N_sensors
        D_s(m1,m2) = (1/n_p).*sum(sum(D{m1,m2}.*Tr));
    end
end
theta1=-20:1:20;

P= tf_music(D_s, n_sources, N_sensors, 2,1, theta1);
 P_CKD(iii,:)=P;%/max(P);

%%% DOA Estimation New spatial averaging

Tr = (D_avg) >= thr;
[F_trace, ~] = find(Tr);
n_p = length(F_trace);
D_s = zeros(N_sensors, N_sensors);
for m1 = 1:N_sensors
    for m2 = 1:N_sensors
        D_s(m1,m2) = (1/n_p).*sum(sum(DD{m1,m2}.*Tr));
    end
end
P=tf_music(D_s, n_sources, N_sensors, 2,1, theta1);
P_spatial(iii,:) = P;%/max(P);

[D_avg,DD,~]=SADTFD_new(X,2,30,length(X)/2);
%%%%%%%%%%%%%%%%%%%%%%%%


% 
 perc=0.4;
 thr = perc*max(max(D_avg));
 Tr = (D_avg) >= thr;
 [F_trace, ~] = find(Tr);
 n_p = length(F_trace);
 D_s = zeros(N_sensors, N_sensors);
 for m1 = 1:N_sensors
     for m2 = 1:N_sensors
         D_s(m1,m2) = (1/n_p).*sum(sum(DD{m1,m2}.*Tr));
     end
 end
% 
% %%% DOA Estimation
 P_tf_music_SADTFD(iii,:)  = tf_music(D_s, n_sources, N_sensors, 2,1, theta1);

end
P_CKD=mean(P_CKD);
P_CKD=P_CKD/max(P_CKD);
P_spatial=mean(P_spatial);
P_spatial=P_spatial/max(P_spatial);
P_tf_music_SADTFD=mean(P_tf_music_SADTFD);
P_tf_music_SADTFD=P_tf_music_SADTFD/max(P_tf_music_SADTFD);
A=zeros(size(theta1))+1;
A(theta1==-4)=-3;
A(theta1==4)=-3;
figure;plot(theta1,log10(P_CKD),'linewidth',3);
hold on;plot(theta1,log10(P_spatial),'r','linewidth',3);
 hold on;plot(theta1,log10(P_tf_music_SADTFD),'k:','linewidth',3);
 hold on;stem(theta1,A);
xlabel('Spatial Spectrum','FontSize',24);
ylabel('Amplitude (dB)','FontSize',24);
legend('STFD based TF MUSIC','Averaged STFD based TF MUSIC','Both averaged and directionally smoothed STFD');
axis([min(theta1) max(theta1)  -3 0],'FontSize',24);
%%% DOA Estimation New spatial averaging+SADTFD

% 
% [D_avg,D1,~]=SADTFD_new(X,2,30,length(X)/2);
% D_avg(D_avg<0)=0;
% perc=0.4;
% thr = perc*max(max(D_avg));
% Tr = abs(D_avg) >= thr;
% [F_trace, ~] = find(Tr);
% n_p = length(F_trace);
% D_s = zeros(N_sensors, N_sensors);
% for m1 = 1:N_sensors
%     for m2 = 1:N_sensors
%         D_s(m1,m2) = (1/n_p).*sum(sum(D1{m1,m2}.*Tr));
%     end
% end
% theta1=-90:0.5:90;
% 
% %%% DOA Estimation
% P_tf_music_SADTFD = tf_music(D_s, n_sources, N_sensors, 2,1, theta1);
% hold on;plot(theta1,P_tf_music_SADTFD,'k:');
