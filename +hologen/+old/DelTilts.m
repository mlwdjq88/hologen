function [ wf1 ] = DelTilts( wf0,mask)
% ����б����������
    [height,width]=size(wf0);
    xx=ones(height,1)*linspace(-1,1,width);
    yy=linspace(-1,1,height)'*ones(1,width);
    wf0s=wf0;
    wf0s(mask~=1)=NaN;
%     wf0s(abs(wf0-mean(wf0(~isnan(wf0))))>3*std(wf0(~isnan(wf0))))=NaN;
%     wf0s(mask==0)=NaN;
    coef=PlaneFit(xx,yy,wf0s);
    wf1 = wf0 - coef(1)*xx- coef(2)*yy - coef(3);
end

function [ coef ] = PlaneFit(xoi,yoi,zoi)
% ƽ����ϣ�z=A*x+B*y+C��zoi�����������Ч��nan��xoi���Ǿ��������Կ�
    vp=~isnan(zoi);
    xi=xoi(vp);
    yi=yoi(vp);
    zi=zoi(vp);

    %ֻ������б��ʱ������ö�������ϣ�Ч���������ã�����ֻ��ƽ�����
    A = [sum(xi.^2) sum(xi.*yi) sum(xi);
        sum(xi.*yi) sum(yi.^2) sum(yi);
        sum(xi) sum(yi) length(xi)];
    B = [sum(xi.*zi);sum(yi.*zi);sum(zi)];
    coef = A\B;
%     out = zi - coef(1)*xi - coef(2)*yi - coef(3);
end