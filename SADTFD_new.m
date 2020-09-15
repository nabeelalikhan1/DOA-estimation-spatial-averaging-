function [Inew,D,orient]= SADTFD_new(S,alpha1,alpha2,N)

%length

%alpha=0.25;
R=3;
%D = mtfd(S, 'WVD');%,0.04,0.04);
D = mtfd(S, 'CKD',1,0.3,0.3);

[MM,~]=size(D);
D_avg = zeros(length(S), length(S));
for mm = 1:MM, D_avg = D{mm,mm} + D_avg; end
DD=D;
D_avg=real(D_avg)/MM;

%D=DD;
M=180;


[X,Y]=meshgrid(-1:2/N:1,-1:2/N:1);

for kk=0:M/R-1
    angle=pi*kk*R/M;
    
    X1=X*cos(angle)-Y*sin(angle);
    Y1=X*sin(angle)+Y*cos(angle);
    
    
    A=exp((-1/2)*(((alpha1*X1).^2)+(alpha2*Y1).^2));
    A=A.*(1-alpha2*alpha2*Y1.^2);
    
    A=A/sum(sum(abs(A)));
    BB{kk+1} = A;
    
end

PQ = paddedsize(size(D_avg));
F = fft2(double(D_avg), PQ(1), PQ(2));

FIabs = fft2(double(abs(D_avg)), PQ(1), PQ(2));

jjjj=0;
for ii=0:M/R-1
    B = BB{ii+1};
    
    H = fft2(double(B), PQ(1), PQ(2));
    H1 = fft2(double(B), PQ(1), PQ(2));
    
    F_fH = H.*F;
    ffi = real(ifft2(F_fH));
    II3(:,:,jjjj+1)=(ffi(round(length(B)/2):end-round(length(B)/2), round(length(B)/2):end-round(length(B)/2)));
    %F_fH = H.*FIabs;
    
    F_fH = H1.*FIabs;
    
    ffi =( abs(ifft2(F_fH))).^2;
    
    II(:,:,jjjj+1)=(ffi(round(length(B)/2):end-round(length(B)/2), round(length(B)/2):end-round(length(B)/2)));
    %  end
    jjjj=jjjj+1;
end

Inew=zeros(size(ffi(1:end/2, 1:end/2)));
[M1,N1]=size(ffi(1:end/2, 1:end/2));
[b,a]=max(II,[],3);
for m=1:M1
    for n=1:N1
        Inew(m,n)=II3(m,n,a(m,n));
        orient(m,n)=a(m,n);
    end
end
%Inew(Inew<0)=0;
for ii=1:MM
    DD{ii,ii}=Inew;
end



lag=MM-1;
for jj=1:lag
    I_avg=zeros(length(S),length(S));
    
    for ii=1:MM-jj
        I1=D{ii,ii+jj};
        I_avg=I_avg+I1;
    end
    I_avg=I_avg/(MM-jj);
    F = fft2(double(I_avg), PQ(1), PQ(2));
    
    jjjj=0;
    for ii=0:M/R-1
        B = BB{ii+1};
        
        H = fft2(double(B), PQ(1), PQ(2));
        F_fH = H.*F;
        ffi = (ifft2(F_fH));
        II3(:,:,jjjj+1)=(ffi(round(length(B)/2):end-round(length(B)/2), round(length(B)/2):end-round(length(B)/2)));
        jjjj=jjjj+1;
    end
    for m=1:M1
        for n=1:N1
            IS(m,n)=II3(m,n,a(m,n));
        end
    end
    
    %IS(Inew<=0)=0;

    
    for ii=1:MM-jj
        DD{ii,ii+jj}=IS;
        DD{ii+jj,ii}=conj(IS);
        
    end
end
%Inew(Inew<0)=0;
for ii=1:MM
    DD{ii,ii}=Inew;
end

D=DD;
end




%end