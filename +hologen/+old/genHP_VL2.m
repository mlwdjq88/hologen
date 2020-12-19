function genHP_VL2()
%  mex CgetPolygonCoords.cpp
%   mex CgetPolygonCoords_polar.cpp
%     mex CdownSamplingUsingRealCoords.cpp
% mex CgetPolygonCoordsUsingPolarCoords.cpp
%  mex CgetPolygonCoordsUsingAdjacentCoords.cpp
% mex CgetPolygonCoordsUsingFlippedAdjacentCoords.cpp
  % mex CgetBoundariesFromLabel.cpp
% mex CsplitConcave.cpp
% mex CreadCor.cpp
%% parameter setting
%current unit: mm
% Ts=[0.02,0.03,0.06,0.1,0.2;0.07,0.1,0.2,0.4,10];
% global stopsign
fs=[15,20,35,50,75];
Ts(1,:)=linspace(0.05,0.6,9);
Ts(2,:)=linspace(0.07,0.8,9);
Ts(3,:)=linspace(0.12,1.5,9);
Ts(4,:)=linspace(0.16,2,9);
Ts(5,:)=linspace(0.25,3,9);
Ts(1,10)=6;
Ts(2,10)=8;
Ts(3,10)=15;
Ts(4,10)=20;
Ts(5,10)=50;
for i=1:1
    for j=10:10
%         for k=1:5
% i=1
% j=5
% k=5
i=2;
j=6;
%             f=5+(i-1)*15;
%             incidentAngle=6+2*(j-1)+10*(2-i);
%             T=Ts(i,k);
% k=5
% f=5+(i-1)*15;
% incidentAngle=6+2*(j-1)+10*(2-i);
% T=Ts(i,k);
f=fs(i);
incidentAngle=6;
T=Ts(i,j);
db=10000000; %set unit to anstrom
lambda=632.8e-6;
delta=2*f*tan(asin(lambda/T));
% NA=0.06;
% R=f*tan(asin(NA));
Nx=500;
Ny=500;
ringSampling=50;
ringSamplingNum=Nx/ringSampling;
% dire=1;
DownSamplingAccuCtrl=0.001;
CoordsAccuCtrl=0.0001;
filename=[num2str(i),num2str(j),'F',num2str(f),'_T',num2str(T),'.gds'];
incidentAngle=incidentAngle/180*pi;
R=tan(incidentAngle)*f/2;
% R=3.8;
%  R=R/2;
RingNum=(sqrt(f^2+R^2)-f)/lambda;
divnum=2*floor(RingNum/ringSamplingNum)+1;
dx=2*R/divnum;
dy=2*R/divnum;
% offset=f*tan(incidentAngle);
% R=offset+R;
[polyPre,polyForm,polyPost]=polyDef();
Mx=divnum;
My=divnum;
fliped=divnum/2+1;
offsetx=8*(i-1)*db;
offsety=8*(j-1)*db;
tic
parfor p=1:Mx
    for q=1:My
%         if p==5&&q==6
%             stopsign=1;
%         else
%             stopsign=0;
%         end
        getCoords(lambda,delta,f,R,db,Nx,Ny,p,q,dx,dy,offsetx,offsety,fliped,incidentAngle,DownSamplingAccuCtrl,CoordsAccuCtrl,polyPre,polyForm,polyPost);
    end
end
toc
%% generate GDS
CreadCor(Mx,My,filename);
delete('*.cor');
    end
end
% end
  winopen(filename);


function B=getCoords(lambda,delta,f,R,db,Nx,Ny,p,q,dx,dy,offsetx,offsety,fliped,incidentAngle,AccuCtrl,CoordsAccuCtrl,polyPre,polyForm,polyPost)
lambda=lambda*db;
delta=delta*db;
f=f*db;
% offset=f*tan(incidentAngle);
R=R*db;
T=lambda/sin(atan(delta/f/2));
% fliped=1.5;
dx=db*dx;
dy=db*dy;
sx=(p-fliped)*dx;
sy=(q-fliped)*dy;
Rx=sx;
Rx(2)=Rx+dx;
Ry=sy;
Ry(2)=Ry+dy;
%%find two non-zero points
[thn,rn] = cart2pol(Rx,0);
rsn=sqrt(rn*R);
xn=rsn.*cos(thn);
for i=1:2
    if getIntensity(xn(i),0,delta,f,lambda,incidentAngle)>4
        temp=linspace(xn(i)-T/8,xn(i)+T/8,50);
        tempI=getIntensity(temp,0,delta,f,lambda,incidentAngle);
        tempR=temp(tempI<4);
        if isempty(tempR)
            break;
        end
        [~,Index]=min(abs(tempR-xn(i)));
        Rx(i)=(tempR(Index(1))/cos(thn(i))).^2/R*cos(thn(i));
    end
end
[x,y]=meshgrid(linspace(Rx(1),Rx(2),Nx*2+1),linspace(Ry(1),Ry(2),Ny*2+1));
[TH,r] = cart2pol(x,y);
rs=sqrt(r*R);
xs=rs.*cos(TH);
ys=rs.*sin(TH);
B=getBoundaries(xs,ys,delta,f,lambda,incidentAngle,CoordsAccuCtrl);
% for i=1:length(B)
%     figure(10),plot(B{i}(:,1),B{i}(:,2)),hold on;
% end
% pause;
Bs=CdownSamplingUsingRealCoords(B,lambda,delta,f,AccuCtrl);
% if p==7&&q==7
% for i=1:length(B)
%     figure(10),plot(Bs{i}(:,1),Bs{i}(:,2)),hold on;
% end
% end
C=CsplitConcave(Bs);
if fliped~=0
    CgetPolygonCoordsUsingAdjacentCoords(C,polyPre,polyForm,polyPost,Ny,p,q,offsetx,offsety,R,0.24*db);
else
    CgetPolygonCoordsUsingFlippedAdjacentCoords(C,polyPre,polyForm,polyPost,Ny,p,q);
end


function B=getBoundaries(xs,ys,delta,f,lambda,incidentAngle,CoordsAccuCtrl)
I=getIntensity(xs,ys,delta,f,lambda,incidentAngle);
I0=I;
th=8;
deltas=0.1;
mask=zeros(size(I));
mask(2:end-1,2:end-1)=1;
I(I<th-3*deltas&mask==0)=0;
I(I>=th-3*deltas&mask==0)=1;
I(I<th-deltas&mask==1)=0;
I(I>=th-deltas&mask==1)=1;
[L,~]=bwlabel(I);
L0=L;
L(I0>th+deltas&mask==1)=0;
L(I0>th+3*deltas&mask==0)=0;
L(1,1)=L0(1,1);
L(1,end)=L0(1,end);
L(end,end)=L0(end,end);
L(end,1)=L0(end,1);
B=getBoundariesFromLabel(L,xs,ys,delta,f,lambda,th,incidentAngle,CoordsAccuCtrl);
% tic
%     B=getBoundariesFromLabels(L,xs,ys,delta,f,lambda,th,CoordsAccuCtrl);
%     B=getBoundariesFromLabels_mex(L,xs,ys,delta,f,lambda,th,CoordsAccuCtrl);
%     s=1
   

function B=getBoundariesFromLabel(L,xs,ys,delta,f,lambda,th,incidentAngle,CoordsAccuCtrl)
% global stopsign
[sr,sc]=size(xs);
T=lambda/sin(atan(delta/f/2));
num=max(L(:));
if num>1
B=cell(num,1);
end
k=1;
for i=1:num
    [y,x]=find(L==i);
    bx=find(x==1|x==sc);
    by=find((x~=1&x~=sc)&(y==1|y==sr));
    xr=xs((x-1)*sr+y);
    yr=ys((x-1)*sr+y);
    [ths,rs]=cart2pol(xr,yr);
    if ~all(ths<0)
        thsN=ths(ths<0);
        thsP=ths(ths>=0);
        if(2*pi+min(thsN)-max(thsP)< min(thsP)-max(thsN))
        ths(ths<0)=ths(ths<0)+2*pi;
        end
    end
    th0=mean(ths);
    r0=mean(rs);
    [x0,y0]=pol2cart(th0,r0);
     
    [xr,yr]=getAccurateCoords(xr,yr,x0,y0,bx,by,delta,f,lambda,th,incidentAngle,CoordsAccuCtrl);
%     if i==19
%            plot(xr,yr,'.'),hold on,pause(0.01)
%     end
    if length(xr)<3
        continue;
    end
    [ths,rs]=cart2pol(xr,yr);
    tempt=0;
    tempr=0;
    angM=false;
    if ~all(ths<0)
        %         ths(ths<-pi/2)=ths(ths<-pi/2)+2*pi;
        thsN=ths(ths<0);
        thsP=ths(ths>=0);
        angM=(2*pi+min(thsN)-max(thsP)< min(thsP)-max(thsN));
        if angM
            ths(ths<0)=ths(ths<0)+2*pi;
            if min(thsP)-max(thsN)<20*pi/180
                tempt=min(thsP);
                tempr=rs(ths==min(thsP));
            end
        elseif (2*pi+min(thsN)-max(thsP))<20*pi/180
            tempt=min(thsN);
            tempr=rs(ths==min(thsN));
        end
%          (2*pi+min(thsN)-max(thsP))<90*pi/180

%         [xr,yr]=pol2cart(tempt,tempr);
    end
    dire=1;
    [thmin,~]=min(ths);
    [thmax,~]=max(ths);
    if std(rs(ths<=thmin+(thmax-thmin)/4))<std(rs(ths>thmax-(thmax-thmin)/4))
        dire=0;
        if tempt~=0&&tempr~=0
            if angM&&(min(thsP)-max(thsN)<20*pi/180)
                tempt=max(thsN);
                tempr=rs(ths==(max(thsN)+2*pi));
            elseif (2*pi+min(thsN)-max(thsP))<20*pi/180
                tempt=max(thsP);
                tempr=rs(ths==max(thsP));
            end
        end
    end
    
    l=r0*(thmax-thmin);
    if max(rs)-min(rs)>2.5*lambda
        div=T/5;
    else
        div=T/10;
    end
    if div>2*r0*pi/180
        div=2*r0*pi/180;
    end
%     if i==61
%         l
%         div
%     end
%      rd=(max(rs)-min(rs));
    if l>div
        Ns=ceil(l/div);
        for m=1:Ns
            if dire==1
                indexm=ths<=thmin+(thmax-thmin)*m/Ns;
            else
                indexm=ths>=thmax-(thmax-thmin)*m/Ns;
            end
            xm=xr(indexm);
%              lw=(max(rs(ths<=thmin+(thmax-thmin)*m/Ns))-min(rs(ths<=thmin+(thmax-thmin)*m/Ns)))/rd;
            if length(xm)<5&&m~=Ns||length(xm)<3%||(lw<0.1)
                continue;
            end
            xr(indexm)=[];
            ym=yr(indexm);
            yr(indexm)=[];
            ths(indexm)=[];
            [thm,rm]=cart2pol(xm,ym);
            if angM
                %                 thm(thm<-pi/2)=thm(thm<-pi/2)+2*pi;
%                 thsN=thm(thm<0);
%                 thsP=thm(thm>=0);
%                 if(2*pi+min(thsN)-max(thsP)< min(thsP)-max(thsN))
                    thm(thm<0)=thm(thm<0)+2*pi;
%                 end
            end
            %             if i==104&&m==4
            %                 s=1;
            %             end
%             if i==61
%                 s=1
%             end
            [xm,ym]=sortCoods(xm,ym,thm,rm,1);
            [thm,rm]=cart2pol(xm,ym);
%             
            if angM
% %                 thm(thm<-pi/2)=thm(thm<-pi/2)+2*pi;
%                 thsN=thm(thm<0);
%                 thsP=thm(thm>=0);
%                 if(2*pi+min(thsN)-max(thsP)< min(thsP)-max(thsN))
                    thm(thm<0)=thm(thm<0)+2*pi;
%                 end
            end
            if dire==1
                [~,Im]=max(thm);
            else
                [~,Im]=min(thm);
            end
            switch Im
                case 1
%                     if abs(rm(1)-rm(2))>abs(rm(Im)-rm(end))
                    if getIntensity((xm(Im)+xm(Im+1))/2,(ym(Im)+ym(Im+1))/2,delta,f,lambda,incidentAngle)>getIntensity((xm(Im)+xm(end))/2,(ym(Im)+ym(end))/2,delta,f,lambda,incidentAngle)
                        xr(end+1,1)=xm(1);
                        xr(end+1,1)=xm(2);
                        yr(end+1,1)=ym(1);
                        yr(end+1,1)=ym(2);
                        ths(end+1,1)=thm(1);
                        ths(end+1,1)=thm(2);
                    else
                        xr(end+1,1)=xm(1);
                        xr(end+1,1)=xm(end);
                        yr(end+1,1)=ym(1);
                        yr(end+1,1)=ym(end);
                        ths(end+1,1)=thm(1);
                        ths(end+1,1)=thm(end);
                    end
                case length(thm)
%                     if abs(rm(end)-rm(1))>abs(rm(end)-rm(end-1))
                    if getIntensity((xm(Im)+xm(1))/2,(ym(Im)+ym(1))/2,delta,f,lambda,incidentAngle)>getIntensity((xm(Im)+xm(Im-1))/2,(ym(Im)+ym(Im-1))/2,delta,f,lambda,incidentAngle)
                        xr(end+1,1)=xm(end);
                        xr(end+1,1)=xm(1);
                        yr(end+1,1)=ym(end);
                        yr(end+1,1)=ym(1);
                        ths(end+1,1)=thm(end);
                        ths(end+1,1)=thm(1);
                    else
                        xr(end+1,1)=xm(end);
                        xr(end+1,1)=xm(end-1);
                        yr(end+1,1)=ym(end);
                        yr(end+1,1)=ym(end-1);
                        ths(end+1,1)=thm(end);
                        ths(end+1,1)=thm(end-1);
                    end
                otherwise
%                     if abs(rm(Im)-rm(Im+1))>abs(rm(Im)-rm(Im-1))
                    if getIntensity((xm(Im)+xm(Im+1))/2,(ym(Im)+ym(Im+1))/2,delta,f,lambda,incidentAngle)>getIntensity((xm(Im)+xm(Im-1))/2,(ym(Im)+ym(Im-1))/2,delta,f,lambda,incidentAngle)
                        xr(end+1,1)=xm(Im);
                        xr(end+1,1)=xm(Im+1);
                        yr(end+1,1)=ym(Im);
                        yr(end+1,1)=ym(Im+1);
                        ths(end+1,1)=thm(Im);
                        ths(end+1,1)=thm(Im+1);
                    else
                        xr(end+1,1)=xm(Im);
                        xr(end+1,1)=xm(Im-1);
                        yr(end+1,1)=ym(Im);
                        yr(end+1,1)=ym(Im-1);
                        ths(end+1,1)=thm(Im);
                        ths(end+1,1)=thm(Im-1);
                    end
            end
            B{k}=[xm,ym];
            k=k+1;
            if m==Ns&&(tempt~=0||tempr~=0)
                 [xt,yt]=pol2cart(tempt,tempr);
                 xr(end+1)=xt;
                 yr(end+1)=yt;
                 B{k}=[xr,yr];
                 k=k+1;
            end
%             if stopsign==1
%                      plot(xm,ym),hold on;pause(0.01)
%             end
        end
    else
        [xr,yr]=sortCoods(xr,yr,ths,rs,1);
        B{k}=[xr,yr];
        k=k+1;
%         if stopsign==1
%              plot(xr,yr),hold on;pause(0.01)
%         end
    end
end
for n=num:-1:k
    B(n)=[];
end


function [xr,yr]=sortCoods(xr,yr,ths,rs,method)
if method==1
    [thmin,Imin]=min(ths);
    [thmax,Imax]=max(ths);
    rmin=rs(Imin);
    rmax=rs(Imax);
    ravg=mean(rs);
    if (rmin-ravg)*(rmax-ravg)>0&&...
            ((abs(rmin-ravg)/(ravg-min(rs))>1/3&&abs(rmax-ravg)/(ravg-min(rs))>1/3)||...
            (abs(rmin-ravg)/(max(rs)-ravg)>1/3&&abs(rmax-ravg)/(max(rs)-ravg)>1/3))
        rmin=ravg;
    end
    sides=calSides(ths,rs,thmin,rmin,thmax,rmax);
    index=sides<0;
    sc=false;
    if mean(sides(~index))~=0
        s=-mean(sides(index))/mean(sides(~index));
        if (s>1)
            s=1/s;
        end
        if s<0.08
            sc=true;
        end
    else
        sc=true;
    end
    if sum(index)<=1||sum(index)>=length(index)-3||sc
    [xr,yr]=sortCoods(xr,yr,ths,rs,0);
    return;
    end
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
else
    x0=mean(xr);
    y0=mean(yr);
    [th,~]=cart2pol(xr-x0,yr-y0);
    [~,Index]=sort(th);
    xr=xr(Index);
    yr=yr(Index);
%     plot(xr,yr),pause(2)
end


function sides=calSides(x,y,x1,y1,x2,y2)
sides=y*(x2-x1)-x*(y2-y1)-y1*x2+x1*y2;


function [xr,yr]=getAccurateCoords(xr,yr,x0,y0,bx,by,delta,f,lambda,th,incidentAngle,accu)
N=length(xr);
[thi,ri]=cart2pol(xr-x0,yr-y0);
for i=1:length(bx)
    [thib,rib]=cart2pol(0,yr(bx(i))-y0);
    thi(bx(i))=thib;
    ri(bx(i))=rib;
    xb(i)=xr(bx(i));
end
for i=1:length(by)
    [thib,rib]=cart2pol(xr(by(i))-x0,0);
    thi(by(i))=thib;
    ri(by(i))=rib;
    yb(i)=yr(by(i));
end
rimin=min(ri);
rimax=max(ri);
r0=ri;
for i=1:N
    Ia=getIntensity(xr(i),yr(i),delta,f,lambda,incidentAngle)-th;
    if abs(Ia)<0.1
        a=ri(i);
        b=a*(Ia/2+th)/th;
        [xs,ys]=pol2cart(thi(i),b);
        Ib=getIntensity(xs+x0,ys+y0,delta,f,lambda,incidentAngle)-th;
        if abs(Ia)>abs(Ib)
            ri(i)=b;
        end
        Im=min(abs(Ia),abs(Ib));
        k=0;
        while Im>accu&&abs(Ia-Ib)>0.0000001&&k<100
            if abs(Ia)>abs(Ib)
                s=(b-a)/(Ib-Ia)*Ib;
                a=b-s;
                [xs,ys]=pol2cart(thi(i),a);
                while abs(getIntensity(xs+x0,ys+y0,delta,f,lambda,incidentAngle)-th)>abs(Ia)&&a~=b
                    s=s/2;
                    a=b-s;
                    [xs,ys]=pol2cart(thi(i),a);
                end
                Ia=abs(getIntensity(xs+x0,ys+y0,delta,f,lambda,incidentAngle)-th);
                ri(i)=a;
            else
                s=(b-a)/(Ib-Ia)*Ia;
                b=a-s;
                [xs,ys]=pol2cart(thi(i),b);
                while abs(getIntensity(xs+x0,ys+y0,delta,f,lambda,incidentAngle)-th)>abs(Ib)&&a~=b
                    s=s/2;
                    b=a-s;
                    [xs,ys]=pol2cart(thi(i),b);
                end
                Ib=abs(getIntensity(xs+x0,ys+y0,delta,f,lambda,incidentAngle)-th);
                ri(i)=b;
            end
            Im=min(abs(Ia),abs(Ib));
            k=k+1;
        end
        %         Im>accu
        %             abs(Ia-Ib)>0.0000001
        %             k<100
    end
end
Index=ri>1.1*rimax|ri<0.9*rimin|abs(ri-r0)./r0>1/4;
% ri(Index)=[];
% thi(Index)=[];
[xr,yr]=pol2cart(thi,ri);
xr=xr+x0;
yr=yr+y0;
for i=1:length(bx)
    xr(bx(i))=xb(i);
end
for i=1:length(by)
    yr(by(i))=yb(i);
end
xr(Index)=[];
yr(Index)=[];
%  plot(xr,yr,'.'),hold on;pause(0.5)

% plot(xr,yr,'.'),hold on,
% plot(x0,y0,'.'),hold on,
% pause(0.01)


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
    %     sum(angle)
    if sum(angle)>0
        sign=1;
    else
        sign=-1;
    end
    convex=find(angle*sign<0);
    while ~isempty(convex)&&length(angle)>2
        
%         if k==434
%         s=1;
%     end
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
y=det([A;B]);
% y=(x1-x0).*(y2-y1)-(x2-x1).*(y1-y0);
x=dot(A,B);
% x=(x1-x0).*(x2-x1)+(y1-y0).*(y2-y1);
angle=atan2(y,x);
% angle=acos(A*B/norm(A,2)/norm(B,2));

function [polyPre,polyForm,polyPost]=polyDef()
polyPre   = [0, 4, 8, 0, 0, 6, 13, 2, 0, 1, 0, 6, 14, 2, 0, 0];
polyPost    = [0, 4, 17, 0];
polyForm  = [0, 44, 16, 3];
str=dec2hex(polyPre,2);
str=str';
str=str(:);
polyPres(1,:)=str(1:8)';
polyPres(2,:)=str(9:16)';
polyPres(3,:)=str(17:24)';
polyPres(4,:)=str(25:32)';
polyPre=hex2dec(polyPres);
str=dec2hex(polyForm);
str=str';
str=str(:);
polyForm=str';
polyForm=hex2dec(polyForm);
str=dec2hex(polyPost);
str=str';
str=str(:);
polyPost=str';
polyPost=hex2dec(polyPost);

function I=getIntensity(xs,ys,delta,f,lambda,incidentAngle)
% E=exp(1i*2*pi/lambda*sqrt((xs-delta/2).^2+ys.^2+f.^2))+...
%     exp(1i*2*pi/lambda*sqrt((xs+delta/2).^2+ys.^2+f.^2));
% A=ones(size(E));
% I=A.*abs(2+E);
% I=I.^2;
PI=3.14159265;
T=lambda/sin(atan(delta/f/2));
offset=f*tan(incidentAngle);
delta2=f*(tan(asin(sin(incidentAngle)+lambda/T)))-offset;
delta1=offset-f*(tan(asin(sin(incidentAngle)-lambda/T)));
% ys=ys+offset;
I=(cos(2*PI/lambda*sqrt((xs-delta/2).^2+(ys-offset).^2+f.^2))+...
    cos(2*PI/lambda*sqrt((xs+delta/2).^2+(ys-offset).^2+f.^2))+...
    2*cos(2*PI*ys/lambda*sin(incidentAngle))).^2+...
    (-sin(2*PI/lambda*sqrt((xs-delta/2).^2+(ys-offset).^2+f.^2))+...
    -sin(2*PI/lambda*sqrt((xs+delta/2).^2+(ys-offset).^2+f.^2))+...
    2*sin(2*PI*ys/lambda*sin(incidentAngle))).^2;
% I=(cos(2*PI/lambda*sqrt((ys-delta1).^2+xs.^2+f.^2))+...
%     cos(2*PI/lambda*sqrt((ys+delta2).^2+xs.^2+f.^2))+...
%     2*cos(2*PI*xs/lambda*sin(incidentAngle))).^2+...
%     (sin(2*PI/lambda*sqrt((ys-delta1).^2+xs.^2+f.^2))+...
%     sin(2*PI/lambda*sqrt((ys+delta2).^2+xs.^2+f.^2))+...
%     2*sin(2*PI*xs/lambda*sin(incidentAngle))).^2;