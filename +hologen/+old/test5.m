clc
clear all
%  mex CgetPolygonCoords.cpp
%   mex CgetPolygonCoords_polar.cpp
%    mex CdownSampling.cpp
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
AccuCtrl=0.005;
[polyPre,polyForm,polyPost]=polyDef();
Mx=1;
My=1;
% tic,
for p=1:Mx
    for q=1:My
        getCoords(lambda,delta,f,R,db,Nx,Ny,p,q,AccuCtrl,polyPre,polyForm,polyPost);
    end
end
% toc
%% generate GDS
CreadCor(Mx,My);
delete('*.cor');
winopen('test.gds')