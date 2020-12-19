function stout=getCoords(lambda,delta,f,R,subR,xOffset,yOffset,db,Nx,Ny,p,q,Mx,My,dx,dy,fliped,incidentAngle,th,hologram)
lambda=lambda*db;
delta=delta*db;
f=f*db;
% offset=f*tan(incidentAngle);
R=R*db;
subR=subR*db;
xOffset = xOffset*db;
yOffset = yOffset*db;
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
[x,y] = meshgrid(Rx,Ry);
[TH,r] = cart2pol(x,y);
rs=sqrt(r*R);
xs=rs.*cos(TH);
ys=rs.*sin(TH);
if all((xs-xOffset).^2+(ys-yOffset).^2>subR^2)%|(xs-xOffset).^2+(ys-yOffset).^2<(0.99*subR)^2|xs>0)
    stout =[];
    return;
end
[x,y]=meshgrid(linspace(Rx(1),Rx(2),Nx*2+1),linspace(Ry(1),Ry(2),Ny*2+1));
[TH,r] = cart2pol(x,y);
rs=sqrt(r*R);
xs=rs.*cos(TH);
ys=rs.*sin(TH);
stout=getBoundaries(xs,ys,delta,f,lambda,incidentAngle,p,q,Mx,My,th,hologram);



function stout=getBoundaries(xs,ys,delta,f,lambda,incidentAngle,p,q,Mx,My,th,hologram)

Ixy=[];
Ixy0=[];
uxy=[];
sns=1;
switch hologram
    case 1
        % QWLSI
        I = hologen.utils.getIntensity_QWLSI(xs,ys,delta,f,lambda,incidentAngle);
    case 2
        % LSI
        I = hologen.utils.getIntensity_LSI(xs,ys,delta,f,lambda,incidentAngle);
end
[sr,sc]=size(xs);
I0=I;
deltas=0.1;

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
    
end
stout.Ixy=Ixy;
stout.Ixy0=Ixy0;
stout.uxy=uxy;
stout.sns=sns;
