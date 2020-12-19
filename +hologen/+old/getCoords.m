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
tic
B=getBoundaries(xs,ys,delta,f,lambda);
toc
Bs=CdownSampling(B,lambda,delta,f,xs,ys,AccuCtrl);
% getPolygonCoordss(Bs,dire,p,q,Ny,polyPre,polyForm,polyPost);
CgetPolygonCoords_polar(Bs,polyPre,polyForm,polyPost,Ny,p,q);
% CgetPolygonCoords(Bs,dire,polyPre,polyForm,polyPost,Ny,p,q);