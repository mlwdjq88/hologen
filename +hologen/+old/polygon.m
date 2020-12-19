function D=polygon(A,dire)
% N=length(A);
[B,~] = sort(A(:,dire));
B2=unique(B);
% C=A(I,:);
k=1;
t=1;
for i=1:length(B2)
    I1=find(A(:,dire)==B2(i));
    L=length(I1);
    if L==1
        %         case 1
        D{k}(t,:)=A(I1,:);
        t=t+1;
        %             if t>3
        %                 D{k+1}(1:2,:)=D{k}(t-2:t-1,:);
        %                 k=k+1;
        %                 t=3;
        %             end
        %         case 2
    else     %             D{k}(t:t+1,:)=A(I1,:);
        temp=A(I1,3-dire);
        [~,I2]=sort(temp);
            if t==3&&D{k}(1,3-dire)>D{k}(2,3-dire)
                D{k}(t,:)=A(I1(I2(1)),:);
                D{k}(t+1,:)=A(I1(I2(end)),:);
            else
                D{k}(t,:)=A(I1(I2(end)),:);
                D{k}(t+1,:)=A(I1(I2(1)),:);
            end
        t=t+2;
        %             if t>3
        %                 D{k+1}(1:2,:)=D{k}(t-2:t-1,:);
        %                 k=k+1;
        %                 t=3;
        %             end
        %         otherwise
        %             temp=A(I1,3-dire);
        %             [~,I2]=sort(temp);
        %             if t==3&&D{k}(1,1)>D{k}(2,1)
        %                 D{k}(t,:)=A(I1(I2(1)),:);
        %                 D{k}(t+1,:)=A(I1(I2(end)),:);
        %             else
        %                 D{k}(t,:)=A(I1(I2(end)),:);
        %                 D{k}(t+1,:)=A(I1(I2(1)),:);
        %             end
        %             t=t+2;
        %             if t>3
        %                 D{k+1}(1:2,:)=D{k}(t-2:t-1,:);
        %                 k=k+1;
        %                 t=3;
        %             end
    end
    if t>3
        D{k}(t,:)=D{k}(1,:);
        if t==4
            D{k}(t+1,:)=D{k}(1,:);
        end
        if i~=length(B2)
            D{k+1}(1:2,:)=D{k}(t-2:t-1,:);
        end
        k=k+1;
        t=3;
    end
    
end