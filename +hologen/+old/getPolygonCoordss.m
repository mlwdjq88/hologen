function getPolygonCoordss(B,dire,p,q,Ny,polyPre,polyForm,polyPost)
filename=strcat('test',num2str(p),'&',num2str(q),'.cor');
allcoords=[];
outputFile=fopen(filename,'wb');
    Ns=length(B);
       for i=1:Ns
            C=polygons_polar(B{i});
            LC=size(C,2);
            if LC>0
                coords=zeros(16,LC);
                coords(1:5,:)=repmat([polyPre;polyForm],1,LC);
                coords(6:2:14,:)=C(1:2:10,:);
                coords(7:2:15,:)=C(2:2:10,:);
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