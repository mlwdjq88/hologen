clc
clear all
lambda=632.8e-6;
R=500e-3;
N=1000;
[x,y]=meshgrid(linspace(-R,R,N));
ph=zeros(N);
Ns=N/10;
for i=1:10
    for j=1:10
        ph0=zeros(Ns);
        k=Ns-(j)*8;
        ph0(1:k,1:k)=i*2*pi/10;
        ph(1+Ns*(i-1):Ns*i,1+Ns*(j-1):Ns*j)=ph0;
    end
end
E=exp(1i*2*pi/lambda*x*sin(0.2/180*pi)+1i*ph);
A=ones(N);
%   E=exp(-1i*2*pi/lambda*sqrt(xc.^2+yc.^2+f.^2));
 
%E=exp(1i*2*pi/lambda*sqrt(x.^2+y.^2+f.^2));
% E=exp(1i*2*pi/lambda*xc)+...
%     exp(-1i*2*pi/lambda*xc);
% ph=angle(2+E);
 I=A.*abs(1+E);
 I=I.^2;
 th=2;
 I(I<th)=0;
 I(I>th)=1;
 figure(1),
imshow(I,[])