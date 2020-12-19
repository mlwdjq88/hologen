clc
clear all
for p=1:1
    for q=1:1
lambda=632.8e-6;
delta=0.1570;
f=20;
R=0.75;
% Rx=-500e-3+p*4e-3;
% Ry=500e-3+q*4e-3;
N=2000;
[x,y]=meshgrid(linspace(-R,R,N*2+1),linspace(-R,R,N*2+1));
[TH,r] = cart2pol(x,y) ;
% s=sqrt(linspace(0,R^2,1000));
% ss=[-fliplr(s(2:end)) s];
% [xs,ys]=meshgrid((ss));
rs=sqrt(r*R);
xs=rs.*cos(TH);
ys=rs.*sin(TH);
s=1;
switch s
    case 1
        rc=r;
        xc=x;
        yc=y;
    case 2
        rc=rs;
        xc=xs;
        yc=ys;
end
incidentAngle=6/180*pi;
 offset=f*tan(incidentAngle);
% yc=yc+offset;
E=exp(-1i*2*pi/lambda*sqrt((xc-delta/2).^2+(yc-offset).^2+f.^2))+...
     exp(-1i*2*pi/lambda*sqrt((xc+delta/2).^2+(yc-offset).^2+f.^2))+...
     exp(-1i*2*pi/lambda*sqrt((xc).^2+(yc-offset-delta/2).^2+f.^2))+...
     exp(-1i*2*pi/lambda*sqrt((xc).^2+(yc-offset+delta/2).^2+f.^2));
 
A=ones(size(E));
% A(rc>R)=0;
%   E=exp(-1i*2*pi/lambda*sqrt(xc.^2+yc.^2+f.^2));
 
%E=exp(1i*2*pi/lambda*sqrt(x.^2+y.^2+f.^2));
% E=exp(1i*2*pi/lambda*xc)+...
%     exp(-1i*2*pi/lambda*xc);
% ph=angle(2+E);
 I=A.*abs(4*exp(1i*2*pi/lambda*yc*sin(incidentAngle))+E);
 I=I.^2;
 th=32;
 I(I<th)=0;
 I(I>th)=1;
 I(xc.^2+(yc).^2>R^2)=0;
   figure(1),imshow(I,[])
 [B,L,Ns] = bwboundaries(I);
%  C=B;
 figure(2)
 for i=1:Ns
%      B{i}((B{i}(:,1)==1)|(B{i}(:,1)==2*N+1)|(B{i}(:,2)==1)|(B{i}(:,2)==2*N+1),:)=[];
     C{i}=B{i};
     %
     for j=1:length(B{i})
    C{i}(j,1)=xs((B{i}(j,1)),(B{i}(j,2)));
    C{i}(j,2)=ys((B{i}(j,1)),(B{i}(j,2)));
     end
%      C{i}(abs(C{i}(:,1))>0.003)=[];
     plot(C{i}(:,1),C{i}(:,2)),hold on
 end
 pause(0.5)
    end
end
% L = bwlabel(I);
% imshow(L,[])