function genHL_QWLSI_onAxis_EUV()
%%this version can generate polygon having more than 4 sides.

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
% current unit: mm
cAppPath = fileparts(mfilename('fullpath'));
f=3; % lens focal length
incidentAngle=6; % offset angle
T=0.54; % grating pitch
db=10000000; %set unit to anstrom
lambda=13.5e-6; % wavelength
delta=2*f*tan(asin(lambda/T));
NA=0.0875; 
R=f*tan(asin(NA)); % lens radius
Nx=500;
Ny=500;
obscuration =0; % central obscuration
printResolution = 1e-6; % minimum printable feature 
ringSampling=50;
ringSamplingNum=Nx/ringSampling;
% dire=1;
DownSamplingAccuCtrl=0.001; % ratio
CoordsAccuCtrl=0.0001; % intensity
filenamestr=['F',num2str(f),'_T',num2str(T),'_wl',num2str(lambda)];
filename=[filenamestr,'.gds'];
filename = fullfile(cAppPath, '..', '..', 'Data','gds',filename);
outputFile=hologen.utils.OpenGDS(filename,filenamestr);
incidentAngle=incidentAngle/180*pi;
% R=tan(incidentAngle)*f/2;
%  R=R/2;
RingNum=(sqrt(f^2+R^2)-f)/lambda;
divnum=2*floor(RingNum/ringSamplingNum)+1;
dx=2*R/divnum;
dy=2*R/divnum;
% offset=f*tan(incidentAngle);
% R=offset+R;
Mx=divnum;
My=divnum;
fliped=divnum/2+1;
offsetx=8*0*db; % offset of the lens center
offsety=8*0*db;
th=32; % threshold for generating hololens
tic
%  Mx=3;
%  My=5;

ceout=cell(Mx,My);
% parfor_progress(Mx);
parfor p=1:Mx
%     parfor_progress;
    for q=1:My
        ceout{p,q}=getCoords(lambda,delta,f,R,db,Nx,Ny,p,q,Mx,My,dx,dy,fliped,incidentAngle,th);
    end
end
% parfor_progress(0);
toc
Ixy=[];
uxy=[];
Ixy0=[];
sns=0;
for p=1:Mx
    for q=1:My
        Ixy0=[Ixy0,ceout{p,q}.Ixy0];
        Ixy=[Ixy,ceout{p,q}.Ixy];
        if ~isempty(ceout{p,q}.uxy)
            ceout{p,q}.uxy(:,4)=ceout{p,q}.uxy(:,4)+sns;
        end
        uxy=[uxy;ceout{p,q}.uxy];
        sns=length(Ixy);
    end
end

while ~isempty(uxy)
    t=uxy(1,4);
    temp=find(uxy(:,1)==uxy(1,1)&uxy(:,2)==uxy(1,2)&uxy(:,3)==uxy(1,3)&uxy(:,4)~=t);
    if isempty(temp)
        temp=find(uxy(:,1)==uxy(1,1)&uxy(:,2)==uxy(1,2)&uxy(:,3)==uxy(1,3)&uxy(:,4)==t);
        if isempty(temp)
            break;
        end
        %         if length(temp)==1
        %            uxy(:,:)=uxy([end,1:end-1],:);
        %            continue;
        %         end
        uxy(temp,:)=[];
        temp3=find(uxy(:,4)==t);
        if isempty(temp3)
            Ixy0(end+1)=Ixy(t);
        end
        continue;
    end
    k=uxy(temp,4);
    Ixy(t).xr=[Ixy(t).xr;Ixy(k).xr];
    Ixy(t).yr=[Ixy(t).yr;Ixy(k).yr];
    uxy(temp,:)=[];
    uxy(1,:)=[];
    temp2=find(uxy(:,4)==k);
    if ~isempty(temp2)
        uxy(temp2,4)=t;
    end
    temp3=find(uxy(:,4)==t);
    if isempty(temp3)
        Ixy0(end+1)=Ixy(t);
    end
end

lambda=lambda*db;
printResolution = printResolution*db;
delta=delta*db;
f=f*db;
R=R*db;
Rb=obscuration*db;
for t=length(Ixy0):-1:1
    if ~all((Ixy0(t).xr).^2+(Ixy0(t).yr).^2>Rb.^2&(Ixy0(t).xr).^2+(Ixy0(t).yr).^2<=(R).^2)
        Ixy0(t)=[];
    end
end
%    Ixy0=Ixy0(28);
len=length(Ixy0);
% parfor_progress(len);
for q=1:len
%     parfor_progress;
    cxy=[Ixy0(q).xr,Ixy0(q).yr];
    cxy=unique(cxy,'rows');
    cn=length(cxy);
    k=k+1;
    cx=cxy(:,1);
    cy=cxy(:,2);
%     if q~=70
%         continue;
%     end
    B=hologen.utils.CgetBoundariesFromLabelQWLSI(cx,cy,cn,delta,f,lambda,th,incidentAngle,CoordsAccuCtrl);
      
%     B=getBoundariesFromLabel(Ixy0s,delta,f,lambda,th,incidentAngle,CoordsAccuCtrl);
%       Bs =B;
    Bs=hologen.utils.CdownSamplingUsingRealCoords(B,lambda,delta,f,DownSamplingAccuCtrl);
%     q,
%     plot(Bs{q}(:,1),Bs{q}(:,2)),hold on;
% figure(2),hold off;
    for si=length(Bs):-1:1
%                 figure(2),plot(Bs{si}(:,1),Bs{si}(:,2));hold on
     
        while 1 % remove glitch
            num=length(Bs{si});
            angle=zeros(num,1);
            for js=1:num
                if js==1
                    angle(1)=calAngle(Bs{si}(end,1),Bs{si}(end,2),Bs{si}(1,1),Bs{si}(1,2),Bs{si}(2,1),Bs{si}(2,2));
                elseif js==num
                    angle(num)=calAngle(Bs{si}(num-1,1),Bs{si}(num-1,2),Bs{si}(num,1),Bs{si}(num,2),Bs{si}(1,1),Bs{si}(1,2));
                else
                    angle(js)=calAngle(Bs{si}(js-1,1),Bs{si}(js-1,2),Bs{si}(js,1),Bs{si}(js,2),Bs{si}(js+1,1),Bs{si}(js+1,2));
                end
            end
            da=diff([angle;angle(1)]);
            index=abs(da)>5; % remove glitch, if this value is too small, it may remove correct shapes
            if sum(abs(angle)>3.1)>1
                index=index|abs(angle)>3.1;
            end
            Bs{si}(index,:)=[];
            [~,rss]=cart2pol(Bs{si}(:,1),Bs{si}(:,2));
            if max(rss)-min(rss)<printResolution
                Bs{si}=[];
            end
            if sum(index)==0||length(Bs{si})<3
                break;
            end
        end
    end
    % Bs=CdownSamplingUsingRealCoords(Bs,lambda,delta,f,DownSamplingAccuCtrl);
    
    
    
    %% generate boundary
    for ns=1:length(Bs)
        xy=Bs{ns};
        if isempty(xy)||size(xy,1)<3
            continue;
        end
        if (xy(1,1)^2+xy(1,2)^2)>R^2||(xy(1,1)^2+xy(1,2)^2)<Rb^2
            continue;
        end
        xy(end+1,:)=xy(1,:);
        xy(:,1)=xy(:,1)+offsetx;
        xy(:,2)=xy(:,2)+offsety;
        Np=size(xy,1)-1;%多边形顶点数
        hologen.utils.CreateBoundary(outputFile,xy',Np);
    end
end
parfor_progress(0);

%% finish GDS
hologen.utils.CloseGDS(outputFile);
% end
winopen(filename);


function stout=getCoords(lambda,delta,f,R,db,Nx,Ny,p,q,Mx,My,dx,dy,fliped,incidentAngle,th)
lambda=lambda*db;
delta=delta*db;
f=f*db;
% offset=f*tan(incidentAngle);
R=R*db;
% T=lambda/sin(atan(delta/f/2));
% fliped=1.5;
dx=db*dx;
dy=db*dy;
sx=(p-fliped)*dx;
sy=(q-fliped)*dy;
Rx=sx;
Rx(2)=Rx+dx;
Ry=sy;
Ry(2)=Ry+dy;
[x,y]=meshgrid(linspace(Rx(1),Rx(2),Nx*2+1),linspace(Ry(1),Ry(2),Ny*2+1));
[TH,r] = cart2pol(x,y);
rs=sqrt(r*R);
xs=rs.*cos(TH);
ys=rs.*sin(TH);
stout=getBoundaries(xs,ys,delta,f,lambda,incidentAngle,p,q,Mx,My,th);




function stout=getBoundaries(xs,ys,delta,f,lambda,incidentAngle,p,q,Mx,My,th)
% global Ixy0 Ixy uxy sns
Ixy=[];
Ixy0=[];
uxy=[];
sns=1;
% xynum=length(Ixy);
% if xynum==0
% Ixy(1).xr=[];
% Ixy(1).yr=[];
% end
I=getIntensity(xs,ys,delta,f,lambda,incidentAngle);
[sr,sc]=size(xs);
I0=I;
deltas=0.1;
% mask=zeros(size(I));
% mask(2:end-1,2:end-1)=1;
% I(I<th-3*deltas&mask==0)=0;
% I(I>=th-3*deltas&mask==0)=1;
% if (p==9||p==10)&&q==6
%     s=1
% end
I(I<th-deltas)=0;
I(I>=th-deltas)=1;
[L,num]=bwlabel(I);
L0=L;
L0(I0>th+deltas)=0;
for i=1:num
    %up
    l=size(uxy,1);
    temp=find(L(1,:)==i);
    if ~isempty(temp)&&q~=1
        uxy(end+1,1:4)=[p,q-0.5,temp(1),sns];
        if (temp(end)-temp(1))~=length(temp)-1
            df=find(diff(temp)~=1);
            df=df(1)+1;
            uxy(end+1,1:4)=[p,q-0.5,temp(df),sns];
        end
        %     else
        %         xu(1)=0;
        %         yu(1)=0;
    end
    %down
    temp=find(L(end,:)==i);
    if ~isempty(temp)&&q~=My
        uxy(end+1,1:4)=[p,q+0.5,temp(1),sns];
        if (temp(end)-temp(1))~=length(temp)-1
            df=find(diff(temp)~=1);
            df=df(1)+1;
            uxy(end+1,1:4)=[p,q+0.5,temp(df),sns];
        end
        %     else
        %         xu(2)=0;
        %         yu(2)=0;
    end
    %left
    temp=find(L(:,1)==i);
    if ~isempty(temp)&&p~=1
        uxy(end+1,1:4)=[p-0.5,q,temp(1),sns];
        if (temp(end)-temp(1))~=length(temp)-1
            df=find(diff(temp)~=1);
            df=df(1)+1;
            uxy(end+1,1:4)=[p-0.5,q,temp(df),sns];
        end
        %     else
        %         xu(3)=0;
        %         yu(3)=0;
    end
    %right
    temp=find(L(:,end)==i);
    if ~isempty(temp)&&p~=Mx
        uxy(end+1,1:4)=[p+0.5,q,temp(1),sns];
        if (temp(end)-temp(1))~=length(temp)-1
            df=find(diff(temp)~=1);
            df=df(1)+1;
            uxy(end+1,1:4)=[p+0.5,q,temp(df),sns];
        end
        %     else
        %         xu(4)=0;
        %         yu(4)=0;
    end
    if l==size(uxy,1)
        [y,x]=find(L0==i);
        Ixy0(end+1).xr=xs((x-1)*sr+y);
        Ixy0(end).yr=ys((x-1)*sr+y);
    else
        sns=sns+1;
        [y,x]=find(L0==i);
        Ixy(end+1).xr=xs((x-1)*sr+y);
        Ixy(end).yr=ys((x-1)*sr+y);
    end
    %     sns(xynum+i,1)=sn;
    %     uxy(xynum+i,1:16)=[xu,yu];
end
stout.Ixy=Ixy;
stout.Ixy0=Ixy0;
stout.uxy=uxy;
stout.sns=sns;



function B=getBoundariesFromLabel(Ixy0,delta,f,lambda,th,incidentAngle,CoordsAccuCtrl)
% global stopsign
% [sr,sc]=size(xs);
T=lambda/sin(atan(delta/f/2));
num=length(Ixy0);
if num>1
    B=cell(num,1);
end
k=1;
for i=1:num
    %     [y,x]=find(L==i);
    %     bx=find(x==1|x==sc);
    %     by=find((x~=1&x~=sc)&(y==1|y==sr));
    xr=Ixy0(i).xr;
    yr=Ixy0(i).yr;
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
    x0=mean(xr);
    y0=mean(yr);
    %     [x0,y0]=pol2cart(th0,r0);
    %      th0
    %      r0
    [xr,yr]=getAccurateCoords(xr,yr,ths,rs,x0,y0,delta,f,lambda,th,incidentAngle,CoordsAccuCtrl);
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
    if div>3*r0*pi/180
        div=3*r0*pi/180;
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
            if m==1&&(tempt~=0||tempr~=0)
                Im=find(thm==tempt|thm==tempt+2*pi);
                if isempty(Im)
                    tempt=0;
                    tempr=0;
                else
                    switch Im
                        case 1
                            %                     if abs(rm(1)-rm(2))>abs(rm(Im)-rm(end))
                            if getIntensity((xm(Im)+xm(Im+1))/2,(ym(Im)+ym(Im+1))/2,delta,f,lambda,incidentAngle)>getIntensity((xm(Im)+xm(end))/2,(ym(Im)+ym(end))/2,delta,f,lambda,incidentAngle)
                                tempt(1)=xm(1);
                                tempt(2)=xm(2);
                                tempr(1)=ym(1);
                                tempr(2)=ym(2);
                            else
                                tempt(1)=xm(1);
                                tempt(2)=xm(end);
                                tempr(1)=ym(1);
                                tempr(2)=ym(end);
                            end
                        case length(thm)
                            %                     if abs(rm(end)-rm(1))>abs(rm(end)-rm(end-1))
                            if getIntensity((xm(Im)+xm(1))/2,(ym(Im)+ym(1))/2,delta,f,lambda,incidentAngle)>getIntensity((xm(Im)+xm(Im-1))/2,(ym(Im)+ym(Im-1))/2,delta,f,lambda,incidentAngle)
                                tempt(1)=xm(end);
                                tempt(2)=xm(1);
                                tempr(1)=ym(end);
                                tempr(2)=ym(1);
                            else
                                tempt(1)=xm(end);
                                tempt(2)=xm(end-1);
                                tempr(1)=ym(end);
                                tempr(2)=ym(end-1);
                            end
                        otherwise
                            %                     if abs(rm(Im)-rm(Im+1))>abs(rm(Im)-rm(Im-1))
                            if getIntensity((xm(Im)+xm(Im+1))/2,(ym(Im)+ym(Im+1))/2,delta,f,lambda,incidentAngle)>getIntensity((xm(Im)+xm(Im-1))/2,(ym(Im)+ym(Im-1))/2,delta,f,lambda,incidentAngle)
                                tempt(1)=xm(Im);
                                tempt(2)=xm(Im+1);
                                tempr(1)=ym(Im);
                                tempr(2)=ym(Im+1);
                            else
                                tempt(1)=xm(Im);
                                tempt(2)=xm(Im-1);
                                tempr(1)=ym(Im);
                                tempr(2)=ym(Im-1);
                            end
                    end
                end
            end
            
            B{k,1}=[xm,ym];
            k=k+1;
            if m==Ns&&(tempt(1)~=0||tempr(1)~=0)
                %                  [xt,yt]=pol2cart(tempt,tempr);
                xr(end+1:end+2)=tempt;
                yr(end+1:end+2)=tempr;
                [ths,rs]=cart2pol(xr,yr);
                [xr,yr]=sortCoods(xr,yr,ths,rs,0);
                B{k}=[xr,yr];
                k=k+1;
            end
            %             if stopsign==1
            %                      plot(xm,ym),hold on;pause(0.01)
            %             end
        end
    else
        [xr,yr]=sortCoods(xr,yr,ths,rs,1);
        B{k,1}=[xr,yr];
        k=k+1;
        %         if stopsign==1
        %             plot(xr,yr),hold on;pause(0.01)
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


function [xr,yr]=getAccurateCoords(xr,yr,ths,rs,x0,y0,delta,f,lambda,th,incidentAngle,accu)
N=length(xr);
[thi,ri]=cart2pol(xr-x0,yr-y0);
% [thi,ri]=cart2pol(xr-x0,yr-y0);
% for i=1:length(bx)
%     [thib,rib]=cart2pol(0,yr(bx(i))-y0);
%     thi(bx(i))=thib;
%     ri(bx(i))=rib;
%     xb(i)=xr(bx(i));
% end
% for i=1:length(by)
%     [thib,rib]=cart2pol(xr(by(i))-x0,0);
%     thi(by(i))=thib;
%     ri(by(i))=rib;
%     yb(i)=yr(by(i));
% end
% rimin=min(ri);
% rimax=max(ri);
rsmin=min(rs);
rsmax=max(rs);
% r0=ri;
xr0=xr;
yr0=yr;
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
[xr,yr]=pol2cart(thi,ri);
xr=xr+x0;
yr=yr+y0;
rs=sqrt(xr.^2+yr.^2);

% Index=abs(getIntensity((xr+xr0)/2,(yr+yr0)/2,delta,f,lambda,incidentAngle)-th)>abs(getIntensity(xr0,yr0,delta,f,lambda,incidentAngle)-th);
% xr(Index)=xr0(Index);
% yr(Index)=yr0(Index);
Index=rs>(1.1*rsmax-0.1*rsmin)|rs<(1.1*rsmin-0.1*rsmax)|...
    abs(getIntensity((xr+xr0)/2,(yr+yr0)/2,delta,f,lambda,incidentAngle)-th)>abs(getIntensity(xr0,yr0,delta,f,lambda,incidentAngle)-th);
% Index=ri>(1.1*rsmax-0.1*rsmin)|ri<(1.1*rsmin-0.1*rsmax)|abs(ri-r0)./(rimax-rimin)>1/20|abs(ri-r0)./r0>1/4;
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
PI=3.14159265;
offset=f*tan(incidentAngle);
I=(cos(2*PI/lambda*sqrt((xs-delta/2).^2+(ys-offset).^2+f.^2))+...
    cos(2*PI/lambda*sqrt((xs+delta/2).^2+(ys-offset).^2+f.^2))+...
    cos(2*PI/lambda*sqrt((xs).^2+(ys-offset-delta/2).^2+f.^2))+...
    cos(2*PI/lambda*sqrt((xs).^2+(ys-offset+delta/2).^2+f.^2))+...
    4*cos(2*PI*ys/lambda*sin(incidentAngle))).^2+...
    (-sin(2*PI/lambda*sqrt((xs-delta/2).^2+(ys-offset).^2+f.^2))+...
    -sin(2*PI/lambda*sqrt((xs+delta/2).^2+(ys-offset).^2+f.^2))+...
    -sin(2*PI/lambda*sqrt((xs).^2+(ys-offset-delta/2).^2+f.^2))+...
    -sin(2*PI/lambda*sqrt((xs).^2+(ys-offset+delta/2).^2+f.^2))+...
    4*sin(2*PI*ys/lambda*sin(incidentAngle))).^2;
