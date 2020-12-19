clc
clear all
% pitch=5:10:95;
% thita=6:2:24;
% pitch=5:10:95;
%unit: um
lambda=632.8e-3;
db=10000;
for m=1:5
    for n=1:5
        clear x y
thita=6;
type=n*100+(m-1)*500;
filename=['A',num2str(thita),'_','S',num2str(type),'.gds'];
offsetx=42000+(n-1)*8000;
offsety=(m-1)*8000;
T=lambda/asin(thita/180*pi);
pitch=ceil(type/T)*T;
np=floor(5000/type);
if np>10
    np=10;
end
psy=0;
k=1;
for p=50/100*type/np:50/100*type/np:50/100*type
    psx=0;
    for offset = 1/np:1/np:1
        psx=psx+offset*T;
        if offset>0.5
            x(k,1)=pitch*(offset*np-1);
            x(k,2)=pitch*(offset*np-1)+(offset-0.5)*T;
            x(k,3)=pitch*(offset*np-1)+(offset-0.5)*T;
            x(k,4)=pitch*(offset*np-1);
            x(k,5)=x(k,1);
            y(k,1)=psy;
            y(k,2)=psy;
            y(k,3)=psy+p;
            y(k,4)=psy+p;
            y(k,5)=y(k,1);
%             plot(x(k,:),y(k,:),'b'),hold on,pause(0.01);
            k=k+1;
        end
            while mod(psx,pitch)+T/2<p
                x(k,1)=psx;
                x(k,2)=psx+T/2;
                x(k,3)=psx+T/2;
                x(k,4)=psx;
                x(k,5)=x(k,1);
                y(k,1)=psy;
                y(k,2)=psy;
                y(k,3)=psy+p;
                y(k,4)=psy+p;
                y(k,5)=y(k,1);
%                 plot(x(k,:),y(k,:),'b'),hold on,pause(0.01);
                k=k+1;
                psx=psx+T;
            end
            if mod(psx,pitch)<p
                x(k,1)=psx;
                x(k,2)=(offset*np-1)*pitch+p;
                x(k,3)=(offset*np-1)*pitch+p;
                x(k,4)=psx;
                x(k,5)=x(k,1);
                y(k,1)=psy;
                y(k,2)=psy;
                y(k,3)=psy+p;
                y(k,4)=psy+p;
                y(k,5)=y(k,1);
%                 plot(x(k,:),y(k,:),'b'),hold on,pause(0.01);
                k=k+1;
            end
            if mod(p,T)/T<0.5
                x(k,1)=(offset*np-1)*pitch+p;
                x(k,2)=(offset*np-1)*pitch+p+T/2-mod(p,T);
                x(k,3)=(offset*np-1)*pitch+p+T/2-mod(p,T);
                x(k,4)=(offset*np-1)*pitch+p;
                x(k,5)=x(k,1);
                y(k,1)=psy;
                y(k,2)=psy;
                y(k,3)=psy+p;
                y(k,4)=psy+p;
                y(k,5)=y(k,1);
%                 plot(x(k,:),y(k,:),'b'),hold on,pause(0.01);
                k=k+1;
            end
            psx=(offset*np-1)*pitch+p-mod(p,T)+T;
            while psx+T/2<(offset*np)*pitch
                x(k,1)=psx;
                x(k,2)=psx+T/2;
                x(k,3)=psx+T/2;
                x(k,4)=psx;
                x(k,5)=x(k,1);
                y(k,1)=psy;
                y(k,2)=psy;
                y(k,3)=psy+p;
                y(k,4)=psy+p;
                y(k,5)=y(k,1);
%                 plot(x(k,:),y(k,:),'b'),hold on,pause(0.01);
                k=k+1;
                psx=psx+T;
            end
%             if psx/pitch<offset*10
%                 x(k,1)=psx;
%                 x(k,2)=(offset*10)*pitch;
%                 x(k,3)=(offset*10)*pitch;
%                 x(k,4)=psx;
%                 x(k,5)=x(k,1);
%                 y(k,1)=psy;
%                 y(k,2)=psy;
%                 y(k,3)=psy+p;
%                 y(k,4)=psy+p;
%                 y(k,5)=y(k,1);
%                 plot(x(k,:),y(k,:),'b'),hold on,pause(0.05);
%                 k=k+1;
%             end
            psx=offset*np*pitch;
    end
     psx=0;
            psy=psy+p;
            while psx+T/2<np*pitch
                x(k,1)=psx;
                x(k,2)=psx+T/2;
                x(k,3)=psx+T/2;
                x(k,4)=psx;
                x(k,5)=x(k,1);
                y(k,1)=psy;
                y(k,2)=psy;
                y(k,3)=psy+pitch-p;
                y(k,4)=psy+pitch-p;
                y(k,5)=y(k,1);
%                 plot(x(k,:),y(k,:),'b'),hold on,pause(0.01);
                k=k+1;
                psx=psx+T;
            end
            psy=psy+pitch-p;
end
psx=0;
psy=0;
while psx>-pitch
    psx=psx-T;
    x(k,1)=psx;
    x(k,2)=psx+T/2;
    x(k,3)=psx+T/2;
    x(k,4)=psx;
    x(k,5)=x(k,1);
    y(k,1)=-pitch;
    y(k,2)=-pitch;
    y(k,3)=pitch*np;
    y(k,4)=pitch*np;
    y(k,5)=y(k,1);    
    k=k+1;
end
 while psx+T/2<np*pitch
    x(k,1)=psx;
    x(k,2)=psx+T/2;
    x(k,3)=psx+T/2;
    x(k,4)=psx;
    x(k,5)=x(k,1);
    y(k,1)=-pitch;
    y(k,2)=-pitch;
    y(k,3)=0;
    y(k,4)=0;
    y(k,5)=y(k,1);    
    k=k+1;
    psx=psx+T;
 end
centershift=-(np-1)/2*pitch;
offsetx=offsetx+centershift;
offsety=offsety+centershift;
maskGDS(filename,db*x,db*y,db*offsetx,db*offsety);
    end
end
winopen(filename)