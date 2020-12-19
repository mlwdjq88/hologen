function B=downSampling(B,lambda,delta,f,xs,ys)
N=length(B);
for i=1:N
    Rs=xs(B{i}(:,1),B{i}(:,2)).^2+ys(B{i}(:,1),B{i}(:,2)).^2;
    Rs=Rs(:);
    K=mean(sqrt(Rs+f^2+delta^2/4)/lambda-pi/4);
    lv=(floor(K)+pi/4)^2*lambda^2-f^2-delta^2/4;
    if lv>=0
        localT=sqrt((ceil(K)+pi/4)^2*lambda^2-f^2-delta^2/4)-sqrt(lv);
    else
        localT=sqrt((ceil(K)+pi/4)^2*lambda^2-f^2-delta^2/4);
    end
    k=1;
    while k<length(B{i})-1
    X1=xs(B{i}(k,1),B{i}(k,2));
    Y1=ys(B{i}(k,1),B{i}(k,2));
        for j=length(B{i}):-1:k+2
            X2=xs(B{i}(j,1),B{i}(j,2));
            Y2=ys(B{i}(j,1),B{i}(j,2));
            X0=xs(B{i}(k+1:j-1,1),B{i}(k+1:j-1,2));
            Y0=ys(B{i}(k+1:j-1,1),B{i}(k+1:j-1,2));
            d=calDistance(X1,Y1,X2,Y2,X0,Y0);
            if max(d(:))<=localT/1000
                B{i}(k+1:j-1,:)=NaN;
                k=j;
                break;
            end
        end
        if j==k+2
            break;
        end
    end   
     B{i}(isnan(B{i}(:,1)),:)=[];
end

function d=calDistance(X1,Y1,X2,Y2,X0,Y0)
A = Y2 - Y1;
B = X1 - X2;
C = X2*Y1 - X1*Y2;
d=abs(A*X0+B*Y0+C)/sqrt(A^2+B^2);