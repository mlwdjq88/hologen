clc
clear all
lambda=13.5e-6;
delta=0.02;
f=3;
R=30e-3;
N=1000;
[x,y]=meshgrid(linspace(-R,R,N*2+1));
[TH,r] = cart2pol(x,y) ;
% s=sqrt(linspace(0,R^2,1000));
% ss=[-fliplr(s(2:end)) s];
% [xs,ys]=meshgrid((ss));
rs=sqrt(r*R);
xs=rs.*cos(TH);
ys=rs.*sin(TH);
s=2;
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
E=exp(1i*2*pi/lambda*sqrt((xc-delta/2).^2+yc.^2+f.^2))+...
     exp(1i*2*pi/lambda*sqrt((xc+delta/2).^2+yc.^2+f.^2));
A=ones(2*N+1);
A(rc>R)=0;
%   E=exp(-1i*2*pi/lambda*sqrt(xc.^2+yc.^2+f.^2));
 
%E=exp(1i*2*pi/lambda*sqrt(x.^2+y.^2+f.^2));
% E=exp(1i*2*pi/lambda*xc)+...
%     exp(-1i*2*pi/lambda*xc);
% ph=angle(2+E);
 I=A.*abs(2+E);
 I=I.^2;
 th=max(I(:))/2;
 I(I<th)=0;
 I(I>th)=1;
 figure(1),
imshow(I,[])
figure(2), c=contour(I);
c(:,(((c(1,:)-N).^2+(c(2,:)-N).^2)>(N+1).^2))=[];
c(:,c(1,:)==0|c(2,:)==0)=[];
% cr=(c(1,:).^2+c(2,:).^2)/R;
% cth=atan2(c(1,:),c(2,:));
% c0(1,:)=cr.*cos(cth);
% c0(2,:)=cr.*sin(cth);
% c0(1,:)=xs(round(c(2,:)*(2*N)+c(1,:)));
% c0(2,:)=ys(round(c(2,:)*(2*N)+c(1,:)));
c0=c;
for i=1:length(c)
    if round(c(1,i))>0&&round(c(2,i))>0
    c0(1,i)=xs(round(c(2,i)),round(c(1,i)));
    c0(2,i)=ys(round(c(2,i)),round(c(1,i)));
    end
end
c1=c0;
c2=c0;
c1(1,:)=c1(1,:)+delta/2;
c2(1,:)=c2(1,:)-delta/2;
c3(1,:)=[c1(1,:),c2(1,:)];
c3(2,:)=[c1(2,:),c2(2,:)];

%select coordinates
c3(:,c3(1,:).^2+c3(2,:).^2>R.^2)=[];
% ss=sqrt((c3(1,:)+delta/2).^2+c3(2,:).^2+f^2)-sqrt((c3(1,:)-delta/2).^2+c3(2,:).^2+f^2);

ss=exp(-1i*2*pi/lambda*sqrt((c3(1,:)-delta/2).^2+c3(2,:).^2+f.^2))+...
     exp(-1i*2*pi/lambda*sqrt((c3(1,:)+delta/2).^2+c3(2,:).^2+f.^2));
 Is=abs(2+ss);
  c3(:,Is<2)=[];
% Is(Is<2)=0;
% Is(Is>2)=1;
%c3(:,abs(sqrt((c3(1,:)+delta/2).^2+c3(2,:).^2+f^2)-sqrt((c3(1,:)-delta/2).^2+c3(2,:).^2+f^2))<0.00001)=[];
figure(3),plot(c1(1,:),c1(2,:),'.');
% hold on
% plot(c2(1,:),c2(2,:),'.');




        





