function B=getBoundaries(xs,ys,delta,f,lambda)
% [L,num]=bwlabel(I);
I=getIntensity(xs,ys,delta,f,lambda);
I0=I;
th=8;
delta=0.1;
I(I<th-delta)=0;
I(I>=th-delta)=1;
[L,num]=bwlabel(I);
L0=L;
L0(I0>th+delta)=0;
L(2:end-1,2:end-1)=L0(2:end-1,2:end-1);

% I(I~=0)=1;
B=cell(num,1);
for i=1:num
    [y,x]=find(L==i);
    x0=mean(x);
    y0=mean(y);

    [th,~]=cart2pol(x-x0,y-y0);
    [~,Is]=sort(th);
%     Bx=xs((x(Is)-1)*4001+y(Is));
%     By=ys((x(Is)-1)*4001+y(Is));
    Bx=y(Is);
    By=x(Is);
    B{i}=[Bx,By];
%     plot(Bx,By)
end
%  figure(5),mesh(L)
% I(I<th)=0;
% I(I>th)=1;

% for i=1:num
%     
% end