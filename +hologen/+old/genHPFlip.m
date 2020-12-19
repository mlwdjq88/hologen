function genHPFlip()
%  mex CgetPolygonCoords.cpp
%   mex CgetPolygonCoords_polar.cpp
%    mex CdownSamplingUsingRealCoords.cpp
% mex CgetPolygonCoordsUsingPolarCoords.cpp
 mex CgetPolygonCoordsUsingFlippedAdjacentCoords.cpp
% mex CreadCor.cpp
%% parameter setting
db=10000000; %set unit to anstrom
lambda=13.5e-6;
delta=0.008;
f=3;
NA=0.0875;
R=f*tan(asin(NA));
Nx=1000;
Ny=1000;
% dire=1;
AccuCtrl=0.001;
[polyPre,polyForm,polyPost]=polyDef();
Mx=10;
My=1;
tic,
for p=1:Mx
    for q=1:My
        getCoords(lambda,delta,f,R,db,Nx,Ny,p,q,AccuCtrl,polyPre,polyForm,polyPost);
    end
end
toc
%% generate GDS
CreadCor(Mx,My);
delete('*.cor');
winopen('test.gds');


function T=getCoords(lambda,delta,f,R,db,Nx,Ny,p,q,AccuCtrl,polyPre,polyForm,polyPost)
lambda=lambda*db;
delta=delta*db;
f=f*db;
R=R*db;
sx=db*(p-1)*5e-3;
sy=db*(q-1)*5e-3;
dx=db*5e-3;
dy=db*5e-3;
T=lambda/sin(atan(delta/f/2));
Rx=sx;
Rx(2)=Rx+dx;
Ry=sy;
Ry(2)=Ry+dy;
[x,y]=meshgrid(linspace(Rx(1),Rx(2),Nx*2+1),linspace(Ry(1),Ry(2),Ny*2+1));
[TH,r] = cart2pol(x,y) ;
rs=sqrt(r*R);
xs=rs.*cos(TH);
ys=rs.*sin(TH);
% xs(xs<10e-5)=0;
% ys(ys<10e-5)=0;
B=getBoundaries(xs,ys,delta,f,lambda);
Bs=CdownSamplingUsingRealCoords(B,lambda,delta,f,AccuCtrl);
C=splitConcave(Bs);
CgetPolygonCoordsUsingFlippedAdjacentCoords(C,polyPre,polyForm,polyPost,Ny,p,q);



function B=getBoundaries(xs,ys,delta,f,lambda)
[sr,~]=size(xs);
I=getIntensity(xs,ys,delta,f,lambda);
I0=I;
th=8;
deltas=0.1;
I(I<th-deltas)=0;
I(I>=th-deltas)=1;
[L,num]=bwlabel(I);
L0=L;
L(I0>th+deltas)=0;
L(1,1)=L0(1,1);
L(1,end)=L0(1,end);
L(end,end)=L0(end,end);
L(end,1)=L0(end,1);
B=cell(num,1);
k=1;
for i=1:num
    [y,x]=find(L==i);
    xr=xs((x-1)*sr+y);
    yr=ys((x-1)*sr+y);
    [ths,rs]=cart2pol(xr,yr);
    if ~all(ths<0)
    ths(ths<-pi/2)=ths(ths<-pi/2)+2*pi;
    end
    th0=mean(ths);
    r0=mean(rs);
    [x0,y0]=pol2cart(th0,r0);
%     [thi,ri]=cart2pol(xr-x0,yr-y0);
     [xr,yr]=getAccurateCoords(xr,yr,x0,y0,delta,f,lambda,th,0.0001);
    
%     [xr,yr]=pol2cart(thi,ri);
%     xr=xr+x0;
%     yr=yr+y0;
    [ths,rs]=cart2pol(xr,yr);
    if ~all(ths<0)
    ths(ths<-pi/2)=ths(ths<-pi/2)+2*pi;
    end
    
    
    
    [thmin,Imin]=min(ths);
    [thmax,Imax]=max(ths);
    rmin=rs(Imin);
    rmax=rs(Imax);
    sides=calSides(ths,rs,thmin,rmin,thmax,rmax);
    index=sides<0;
    th1=ths(index);
    xr1=xr(index);
    yr1=yr(index);
    th2=ths(~index);
    xr2=xr(~index);
    yr2=yr(~index);  
    [~,I1]=sort(th1);
    [~,I2]=sort(th2,'descend');
    xr1=xr1(I1);
    xr2=xr2(I2);
    yr1=yr1(I1);
    yr2=yr2(I2);
    xr=[xr1;xr2];
    yr=[yr1;yr2];
%       plot(xr,yr,'.');hold on;
%     Ic=getIntensity(x0,y0,delta,f,lambda);
%     [thi,~]=cart2pol(xr-x0,yr-y0);
%     [~,Is]=sort(thi);
    B{k}=[xr,yr];
    k=k+1;
%      if Ic<th
         
%          plot(x0,y0,'.')
%     end
end

function sides=calSides(x,y,x1,y1,x2,y2)
sides=y*(x2-x1)-x*(y2-y1)-y1*x2+x1*y2;


function [xr,yr]=getAccurateCoords(xr,yr,x0,y0,delta,f,lambda,th,accu)
 N=length(xr);
[thi,ri]=cart2pol(xr-x0,yr-y0);
for i=1:N
    Ia=getIntensity(xr(i),yr(i),delta,f,lambda)-th;
    if abs(Ia)<0.1
            a=ri(i);
            b=a*(Ia+th)/th;
            [xs,ys]=pol2cart(thi(i),b);
    Ib=getIntensity(xs+x0,ys+y0,delta,f,lambda)-th;
    if abs(Ia)>abs(Ib)
        ri(i)=b;
    end
    Im=min(abs(Ia),abs(Ib));
    k=0;
    while Im>accu&&abs(Ia-Ib)>0.0000001&&k<100
            if abs(Ia)>abs(Ib)
                    a=b-(b-a)/(Ib-Ia)*Ib;
                    [xs,ys]=pol2cart(thi(i),a);
                    if abs(getIntensity(xs+x0,ys+y0,delta,f,lambda)-th)<abs(Ia)
                        Ia=abs(getIntensity(xs+x0,ys+y0,delta,f,lambda)-th);
                        ri(i)=a;
                    else
                        break;
                    end
            else
                    b=a-(b-a)/(Ib-Ia)*Ia;
                    [xs,ys]=pol2cart(thi(i),b);
                    if abs(getIntensity(xs+x0,ys+y0,delta,f,lambda)-th)<abs(Ib)
                        Ib=abs(getIntensity(xs+x0,ys+y0,delta,f,lambda)-th);
                        ri(i)=b;
                    else
                        break
                    end
            end
                Im=min(abs(Ia),abs(Ib));
                k=k+1;
    end
    end
end
[xr,yr]=pol2cart(thi,ri);
xr=xr+x0;
yr=yr+y0;


function C=splitConcave(B)
C=cell(2,1);
N=length(B);
k=1;
for i=1:N
    num=length(B{i});
    angle=zeros(num,1);
    for j=1:num
        if j==1
            angle(1)=calAngle(B{i}(end,1),B{i}(end,2),B{i}(1,1),B{i}(1,2),B{i}(2,1),B{i}(2,2));
        elseif j==num
            angle(num)=calAngle(B{i}(num-1,1),B{i}(num-1,2),B{i}(num,1),B{i}(num,2),B{i}(1,1),B{i}(1,2));
        else
            angle(j)=calAngle(B{i}(j-1,1),B{i}(j-1,2),B{i}(j,1),B{i}(j,2),B{i}(j+1,1),B{i}(j+1,2));
        end
    end
%      sum(angle)
    if sum(angle)>0
        sign=1;
    else
        sign=-1;
    end
     convex=find(angle*sign<0);
%      if i==79
%          s=1;
%      end
     while ~isempty(convex)&&length(angle)>2
%          if i==61
%          figure(1),plot(B{i}(:,1),B{i}(:,2))
%          end
         B{i}=circshift(B{i},1-convex(1),1);
         angle=circshift(angle,1-convex(1));
         if angle(2)*sign>0
             C{k}=B{i}(1:3,:);
             k=k+1;
%              plot(C{k-1}(:,1),C{k-1}(:,2)),hold on,pause(0.01)
             B{i}(2,:)=[];
             angle(1)=calAngle(B{i}(end,1),B{i}(end,2),B{i}(1,1),B{i}(1,2),B{i}(2,1),B{i}(2,2));
             angle(2)=[];
             angle(2)=calAngle(B{i}(1,1),B{i}(1,2),B{i}(2,1),B{i}(2,2),B{i}(3,1),B{i}(3,2));
         elseif angle(end)*sign>0
             C{k}(1,:)=B{i}(1,:);
             C{k}(2:3,:)=B{i}(end-1:end,:);
%              plot(C{k-1}(:,1),C{k-1}(:,2)),hold on,pause(0.01)
             k=k+1;
             B{i}(end,:)=[];
             angle(1)=calAngle(B{i}(end,1),B{i}(end,2),B{i}(1,1),B{i}(1,2),B{i}(2,1),B{i}(2,2));
             angle(end)=[];
             angle(end)=calAngle(B{i}(end-1,1),B{i}(end-1,2),B{i}(end,1),B{i}(end,2),B{i}(1,1),B{i}(1,2));
         else
             if all(angle*sign<=0)
                 sign=sign*-1;
             end
             B{i}=circshift(B{i},-1,1);
             angle=circshift(angle,-1);
         end
          convex=find(angle*sign<0); 
     end
     C{k}=B{i};
     k=k+1;
%       plot(C{k-1}(:,1),C{k-1}(:,2)),hold on,pause(0.01)
end

function angle=calAngle(x0,y0,x1,y1,x2,y2)
A=[x1-x0,y1-y0];
B=[x2-x1,y2-y1];
angle=atan2d(det([A;B]),dot(A,B));
% angle=acos(A*B/norm(A,2)/norm(B,2));


function I=getIntensity(xs,ys,delta,f,lambda)
E=exp(1i*2*pi/lambda*sqrt((xs-delta/2).^2+ys.^2+f.^2))+...
    exp(1i*2*pi/lambda*sqrt((xs+delta/2).^2+ys.^2+f.^2));
A=ones(size(E));
I=A.*abs(2+E);
I=I.^2;