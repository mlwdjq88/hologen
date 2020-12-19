function B=getBoundariesFromLabels(L,xs,ys,delta,f,lambda,th,CoordsAccuCtrl)
[sr,~]=size(xs);
T=lambda/sin(atan(delta/f/2));
num=max(L(:));
B=cell(1000,1);
for i=1:1000
    B{i}=0;
end
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
    %     plot(xr,yr,'.'),hold on,pause(0.01)
    [xr,yr]=getAccurateCoords(xr,yr,x0,y0,delta,f,lambda,th,CoordsAccuCtrl);
    if length(xr)<3
        continue;
    end
    [ths,rs]=cart2pol(xr,yr);
    if ~all(ths<0)
        ths(ths<-pi/2)=ths(ths<-pi/2)+2*pi;
    end
    [thmin,~]=min(ths);
    [thmax,~]=max(ths);
    l=r0*(thmax-thmin);
    if max(rs)-min(rs)>2.5*lambda
        div=T/5;
    else
        div=T/10;
    end
    if l>div
        Ns=ceil(l/div);
        for m=1:Ns
            xm=xr(ths<=thmin+(thmax-thmin)*m/Ns);
            if length(xm)<5
                continue;
            end
            xr(ths<=thmin+(thmax-thmin)*m/Ns)=[];
            ym=yr(ths<=thmin+(thmax-thmin)*m/Ns);
            yr(ths<=thmin+(thmax-thmin)*m/Ns)=[];
            ths(ths<=thmin+(thmax-thmin)*m/Ns)=[];
            [thm,rm]=cart2pol(xm,ym);
            if ~all(thm<0)
                thm(thm<-pi/2)=thm(thm<-pi/2)+2*pi;
            end
            %             if i==104&&m==4
            %                 s=1;
            %             end
            [xm,ym]=sortCoods(xm,ym,thm,rm,1);
            [thm,rm]=cart2pol(xm,ym);
            if ~all(thm<0)
                thm(thm<-pi/2)=thm(thm<-pi/2)+2*pi;
            end
            [~,Im]=max(thm);
            rN=length(xr);
            xrs=zeros(length(xr)+2,1);
            xrs(1:end-2,1)=xr;
            yrs=zeros(length(yr)+2,1);
            yrs(1:end-2,1)=yr;
            thss=zeros(length(ths)+2,1);
            thss(1:end-2,1)=ths;
            xr=xrs;
            yr=yrs;
            ths=thss;
            
            switch Im
                case 1
                    if abs(rm(1)-rm(2))>abs(rm(Im)-rm(end))
                        xr(rN+1,1)=xm(1);
                        xr(rN+2,1)=xm(2);
                        yr(rN+1,1)=ym(1);
                        yr(rN+2,1)=ym(2);
                        ths(rN+1,1)=thm(1);
                        ths(rN+2,1)=thm(2);
                    else
                        xr(rN+1,1)=xm(1);
                        xr(rN+2,1)=xm(end);
                        yr(rN+1,1)=ym(1);
                        yr(rN+2,1)=ym(end);
                        ths(rN+1,1)=thm(1);
                        ths(rN+2,1)=thm(end);
                    end
                case length(thm)
                    if abs(rm(end)-rm(1))>abs(rm(end)-rm(end-1))
                        xr(rN+1,1)=xm(end);
                        xr(rN+2,1)=xm(1);
                        yr(rN+1,1)=ym(end);
                        yr(rN+2,1)=ym(1);
                        ths(rN+1,1)=thm(end);
                        ths(rN+2,1)=thm(1);
                    else
                        xr(rN+1,1)=xm(end);
                        xr(rN+2,1)=xm(end-1);
                        yr(rN+1,1)=ym(end);
                        yr(rN+2,1)=ym(end-1);
                        ths(rN+1,1)=thm(end);
                        ths(rN+2,1)=thm(end-1);
                    end
                otherwise
                    if abs(rm(Im)-rm(Im+1))>abs(rm(Im)-rm(Im-1))
                        xr(rN+1,1)=xm(Im);
                        xr(rN+2,1)=xm(Im+1);
                        yr(rN+1,1)=ym(Im);
                        yr(rN+2,1)=ym(Im+1);
                        ths(rN+1,1)=thm(Im);
                        ths(rN+2,1)=thm(Im+1);
                    else
                        xr(rN+1,1)=xm(Im);
                        xr(rN+2,1)=xm(Im-1);
                        yr(rN+1,1)=ym(Im);
                        yr(rN+2,1)=ym(Im-1);
                        ths(rN+1,1)=thm(Im);
                        ths(rN+2,1)=thm(Im-1);
                    end
            end
            B{k}=[xm,ym];
            k=k+1;
            %                         plot(xm,ym),hold on;pause(0.01)
        end
    else
        [xr,yr]=sortCoods(xr,yr,ths,rs,1);
        B{k}=[xr,yr];
        k=k+1;
        %                  plot(xr,yr),hold on;pause(0.01)
    end
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
end


function sides=calSides(x,y,x1,y1,x2,y2)
sides=y*(x2-x1)-x*(y2-y1)-y1*x2+x1*y2;


function [xr,yr]=getAccurateCoords(xr,yr,x0,y0,delta,f,lambda,th,accu)
N=length(xr);
[thi,ri]=cart2pol(xr-x0,yr-y0);
rimin=min(ri);
rimax=max(ri);
r0=ri;
for i=1:N
    Ia=getIntensity(xr(i),yr(i),delta,f,lambda)-th;
    if abs(Ia)<0.1
        a=ri(i);
        b=a*(Ia/2+th)/th;
        [xs,ys]=pol2cart(thi(i),b);
        Ib=getIntensity(xs+x0,ys+y0,delta,f,lambda)-th;
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
                while abs(getIntensity(xs+x0,ys+y0,delta,f,lambda)-th)>abs(Ia)&&a~=b
                    s=s/2;
                    a=b-s;
                    [xs,ys]=pol2cart(thi(i),a);
                end
                Ia=abs(getIntensity(xs+x0,ys+y0,delta,f,lambda)-th);
                ri(i)=a;
            else
                s=(b-a)/(Ib-Ia)*Ia;
                b=a-s;
                [xs,ys]=pol2cart(thi(i),b);
                while abs(getIntensity(xs+x0,ys+y0,delta,f,lambda)-th)>abs(Ib)&&a~=b
                    s=s/2;
                    b=a-s;
                    [xs,ys]=pol2cart(thi(i),b);
                end
                Ib=abs(getIntensity(xs+x0,ys+y0,delta,f,lambda)-th);
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
Index=ri>1.1*rimax|ri<0.9*rimin|abs(ri-r0)./r0>1/2;
ri(Index)=[];
thi(Index)=[];
[xr,yr]=pol2cart(thi,ri);
xr=xr+x0;
yr=yr+y0;
% plot(xr,yr,'.'),hold on,
% plot(x0,y0,'.'),hold on,
% pause(0.01)



function I=getIntensity(xs,ys,delta,f,lambda)
% E=exp(1i*2*pi/lambda*sqrt((xs-delta/2).^2+ys.^2+f.^2))+...
%     exp(1i*2*pi/lambda*sqrt((xs+delta/2).^2+ys.^2+f.^2));
% A=ones(size(E));
% I=A.*abs(2+E);
% I=I.^2;
PI=3.14159265;
I=(cos(2*PI/lambda*sqrt((xs-delta/2).^2+ys.^2+f.^2))+cos(2*PI/lambda*sqrt((xs+delta/2).^2+ys.^2+f.^2))+2).^2+...
    (sin(2*PI/lambda*sqrt((xs-delta/2).^2+ys.^2+f.^2))+sin(2*PI/lambda*sqrt((xs+delta/2).^2+ys.^2+f.^2))).^2;