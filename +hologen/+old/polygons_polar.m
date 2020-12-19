function D=polygons_polar(A)
% [TH,r] = cart2pol(A(1:end-1,1),A(1:end-1,2));
% [THs,I]=sort(TH);
% rs=r(I);
L=length(A)-1;
if L<3
    D=[];
    return;
elseif mod(L,2)==0
    A(end,:)=[];
end
num=ceil((L-2)/2);
D=zeros(10,num);
% [X,Y]=pol2cart(THs,rs);
% D(1,:)=A(1,1);
% D(2,:)=A(1,2);
% D(9,:)=A(1,1);
% D(10,:)=A(1,2);
for i=1:num
    D(1,i)=A(i,1);
    D(2,i)=A(i,2);
    D(3,i)=A(end-i+1,1);
    D(4,i)=A(end-i+1,2);
    D(5,i)=A(end-i,1);
    D(6,i)=A(end-i,2);
    D(7,i)=A(i+1,1);
    D(8,i)=A(i+1,2);
    D(9,i)=A(i,1);
    D(10,i)=A(i,2);
end

% 
% [B,~] = sort(A(:,dire));
% B2=unique(B);
% LB=length(B);
% LB2=length(B2);
% if LB==LB2
%     LD=LB2-2;
% else
%     LD=LB2-1;
% end
% D=zeros(10,LD);
% if LD<=0
%     return;
% end
% k=1;
% t=1;
% for i=1:LB2
%     I1=find(A(:,dire)==B2(i));
%     L=length(I1);
%     if L==1
%         D(2*t-1:2*t,k)=A(I1(1),:)';
%         t=t+1;
%     else
%         temp=A(I1,3-dire);
%         [~,I2]=sort(temp);
%         if t==3&&D(3-dire,k)>D(5-dire,k)
%             D(2*t-1:2*t,k)=A(I1(I2(1)),:)';
%             D(2*t+1:2*t+2,k)=A(I1(I2(end)),:)';
%         else
%             D(2*t-1:2*t,k)=A(I1(I2(end)),:)';
%             D(2*t+1:2*t+2,k)=A(I1(I2(1)),:)';
%         end
%         t=t+2;
%     end
%     if t>3
%         D(2*t-1:2*t,k)=D(1:2,k);
%         if t==4
%             D(2*t+1:2*t+2,k)=D(1:2,k);
%         end
%         if i~=LB2
% %              D(1:4,k+1)=D(2*t-5:2*t-2,k);
%              ts=CalSlope(D(2*t-5:2*t,k));
%              if ts==-1
%                  D(1:4,k+1)=D(2*t-5:2*t-2,k);
%              else
%                  D(1:4,k+1)=D(2*t-3:2*t,k);
%              end
%         end
% %         plot(D(1:2:end,k),D(2:2:end,k)),hold on;
% %         pause(1);
%         k=k+1;
%         t=3;
%     end
% end
% if k-1<LD
%     D(:,k:end)=[];
% end
% 
% function t=CalSlope(xy)
% x=xy(1:2:end);
% y=xy(2:2:end);
% if abs(atan2(y(1)-y(2),x(1)-x(2)))>abs(atan2(y(3)-y(2),x(3)-x(2)))
%     t=1;
% else
%     t=-1;
% end
