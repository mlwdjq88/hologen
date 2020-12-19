clc
clear all
% pitch=5:10:95;
% thita=6:2:24;
% pitch=5:10:95;
%unit: um
lambda=632.8e-3;
db=10000;
for m=6:10
    for n=1:5
        clear x y
        thita=6;
        type=500;
        
        offsetx=42000+(n-1)*8000;
        offsety=(m-1)*8000;
        T=lambda/asin(thita/180*pi);
        pitch=ceil(type/T)*T;
        np=1;
        psy=0;
        k=1;
        p=20+(m-6)*50;
        psx=0;
        offset = n/6;
        filename=['S',num2str(p),'_','O',num2str(offset),'.gds'];
        psx=psx+offset*T;
        if offset>0.5
            x(k,1)=pitch*(0);
            x(k,2)=pitch*(0)+(offset-0.5)*T;
            x(k,3)=pitch*(0)+(offset-0.5)*T;
            x(k,4)=pitch*(0);
            x(k,5)=x(k,1);
            y(k,1)=psy;
            y(k,2)=psy;
            y(k,3)=psy+pitch;
            y(k,4)=psy+pitch;
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
            y(k,3)=psy+pitch;
            y(k,4)=psy+pitch;
            y(k,5)=y(k,1);
            %                 plot(x(k,:),y(k,:),'b'),hold on,pause(0.01);
            k=k+1;
            psx=psx+T;
        end
        if mod(psx,pitch)<p
            x(k,1)=psx;
            x(k,2)=(0)*pitch+p;
            x(k,3)=(0)*pitch+p;
            x(k,4)=psx;
            x(k,5)=x(k,1);
            y(k,1)=psy;
            y(k,2)=psy;
            y(k,3)=psy+pitch;
            y(k,4)=psy+pitch;
            y(k,5)=y(k,1);
            %                 plot(x(k,:),y(k,:),'b'),hold on,pause(0.01);
            k=k+1;
        end
        if mod(p,T)/T<0.5
            x(k,1)=(0)*pitch+p;
            x(k,2)=(0)*pitch+p+T/2-mod(p,T);
            x(k,3)=(0)*pitch+p+T/2-mod(p,T);
            x(k,4)=(0)*pitch+p;
            x(k,5)=x(k,1);
            y(k,1)=psy;
            y(k,2)=psy;
            y(k,3)=psy+pitch;
            y(k,4)=psy+pitch;
            y(k,5)=y(k,1);
            %                 plot(x(k,:),y(k,:),'b'),hold on,pause(0.01);
            k=k+1;
        end
        psx=(0)*pitch+p-mod(p,T)+T;
        while psx+T/2<(1)*pitch
            x(k,1)=psx;
            x(k,2)=psx+T/2;
            x(k,3)=psx+T/2;
            x(k,4)=psx;
            x(k,5)=x(k,1);
            y(k,1)=psy;
            y(k,2)=psy;
            y(k,3)=psy+pitch;
            y(k,4)=psy+pitch;
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
        psx=-1*pitch;
        
        
        
        while psx+T/2<0
            x(k,1)=psx;
            x(k,2)=psx+T/2;
            x(k,3)=psx+T/2;
            x(k,4)=psx;
            x(k,5)=x(k,1);
            y(k,1)=psy;
            y(k,2)=psy;
            y(k,3)=psy+pitch;
            y(k,4)=psy+pitch;
            y(k,5)=y(k,1);
            %                 plot(x(k,:),y(k,:),'b'),hold on,pause(0.01);
            k=k+1;
            psx=psx+T;
        end
        % psx=0;
        % psy=0;
        % while psx>-pitch
        %     psx=psx-T;
        %     x(k,1)=psx;
        %     x(k,2)=psx+T/2;
        %     x(k,3)=psx+T/2;
        %     x(k,4)=psx;
        %     x(k,5)=x(k,1);
        %     y(k,1)=-pitch;
        %     y(k,2)=-pitch;
        %     y(k,3)=pitch*np;
        %     y(k,4)=pitch*np;
        %     y(k,5)=y(k,1);
        %     k=k+1;
        % end
        %  while psx+T/2<np*pitch
        %     x(k,1)=psx;
        %     x(k,2)=psx+T/2;
        %     x(k,3)=psx+T/2;
        %     x(k,4)=psx;
        %     x(k,5)=x(k,1);
        %     y(k,1)=-pitch;
        %     y(k,2)=-pitch;
        %     y(k,3)=0;
        %     y(k,4)=0;
        %     y(k,5)=y(k,1);
        %     k=k+1;
        %     psx=psx+T;
        %  end
        
%         offsetx=offsetx;
        offsety=offsety-0.5*pitch;
        maskGDS(filename,db*x,db*y,db*offsetx,db*offsety);
        % end
    end
end
winopen(filename)