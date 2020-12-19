function getPolygonCoords(B,dire,p,q,Ny,xs,ys,polyPre,polyForm,polyPost)
% polyPre   = [0, 4, 8, 0, 0, 6, 13, 2, 0, 1, 0, 6, 14, 2, 0, 0];
% polyPost    = [0, 4, 17, 0];
% polyForm  = [0, 44, 16, 3];
% str=dec2hex(polyPre,2);
% str=str';
% str=str(:);
% polyPres(1,:)=str(1:8)';
% polyPres(2,:)=str(9:16)';
% polyPres(3,:)=str(17:24)';
% polyPres(4,:)=str(25:32)';
% polyPre=hex2dec(polyPres);
% str=dec2hex(polyForm);
% str=str';
% str=str(:);
% polyForm=str';
% polyForm=hex2dec(polyForm);
% str=dec2hex(polyPost);
% str=str';
% str=str(:);
% polyPost=str';
% polyPost=hex2dec(polyPost);
filename=strcat('test',num2str(p),'&',num2str(q),'.cor');
% fid=fopen(filename,'wb');
% filename='test.gds';
allcoords=[];
outputFile=fopen(filename,'wb');
    Ns=length(B);
       for i=1:Ns
            C=polygons(B{i},dire);
            LC=size(C,2);
            if LC>0
                coords=zeros(16,LC);
                coords(1:5,:)=repmat([polyPre;polyForm],1,LC);
                coords(6:2:14,:)=xs((2*Ny+1)*(C(2:2:10,:)-1)+C(1:2:10,:));
                coords(7:2:15,:)=ys((2*Ny+1)*(C(2:2:10,:)-1)+C(1:2:10,:));
                coords(16,:)=repmat(polyPost,1,LC);
                coords2=coords;
                coords2(6:2:14,:)=-coords2(6:2:14,:);
                coords3=coords2;
                coords3(7:2:15,:)=-coords3(7:2:15,:);
                coords4=coords3;
                coords4(6:2:14,:)=-coords4(6:2:14,:);
                allcoords=[allcoords;coords(:);coords2(:);coords3(:);coords4(:)];
%                 plot(coords(6:2:14,:),coords(7:2:15,:)),hold on
%                 error=CalError(coords(6:2:14,:),coords(7:2:15,:),delta,f,lambda,th);
%                 errors=[errors;error(:)];
                
            end
       end
       fwrite(outputFile,allcoords, 'int32','b');
        fclose(outputFile);
%                  et=toc*(num+1-p)*(num+1-q)/60
%           errors(errors>3*mean(errors))=[];
%     errors(errors>mean(errors)+3*std(errors))=[];
%       er(p,q)=mean(errors);
%       plot(errors(:),'.');hold on