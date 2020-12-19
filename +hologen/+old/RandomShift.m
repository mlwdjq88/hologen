% ��������㷨
% 2009-10-21 ������Ч
% ��Դ��Advanced iterative algorithm for phase extraction of randomly phase-shifted interferograms
% ��������Ǹ���ͼ�Ҷ�I�ͳ�ʼλ������delta��delta��ֵ�������������������ȣ��Ҵ�С˳��Ҫ��I��˳���Ӧ��
% �������limit�ǵ�����������
% delta phase�ȶ���һά������I�Ƕ�ά��������1��ά����֡������2��ά��������

function [phase delta k] = RandomShift(I,delta,limit)
    for k = 1:limit
        % step1��ͨ�õ������㷨��delta��Ϊ��֪�����phase
        M = length(delta);
        A = [M sum(cos(delta)) sum(sin(delta));
           sum(cos(delta)) sum(cos(delta).^2) sum(cos(delta).*sin(delta));
           sum(sin(delta)) sum(cos(delta).*sin(delta)) sum(sin(delta).^2)];

        B1 = sum(I,1);

        t = I;
        for i = 1:M
            t(i,:) = I(i,:)*cos(delta(i));
        end
        B2 = sum(t,1);

        t2 = I;
        for i = 1:M
            t2(i,:) = I(i,:)*sin(delta(i));
        end
        B3 = sum(t2,1);

        B = [B1;B2;B3];

        X = A\B;

        phase = atan2(-X(3,:),X(2,:));
            %ע����Ҫ��atan���������ȡֵ��Χ������

        % step2��phase��Ϊ��֪�����delta
        N = length(phase);
        Ap = [N sum(cos(phase)) sum(sin(phase));
           sum(cos(phase)) sum(cos(phase).^2) sum(cos(phase).*sin(phase));
           sum(sin(phase)) sum(cos(phase).*sin(phase)) sum(sin(phase).^2)];

        B1p = sum(I,2);

        tp = I;
        for j = 1:N
            tp(:,j) = I(:,j)*cos(phase(j));
        end
        B2p = sum(tp,2);

        t2p = I;
        for j = 1:N
            t2p(:,j) = I(:,j)*sin(phase(j));
        end
        B3p = sum(t2p,2);

        Bp = [B1p B2p B3p];

        Xp = Ap\Bp';

        delta_p = atan2(-Xp(3,:),Xp(2,:));  
            %ע����Ҫ��atan���������ȡֵ��Χ������

        % step3�������ж�
        res = sum( abs( (delta_p(2:M)-delta_p(1)) - (delta(2:M)-delta(1)) ) );
        res = min([abs(res - pi*floor(res/pi)) abs(res - pi*ceil(res/pi))]);
        if res < 1e-4
            break;
        else
            delta = delta_p;
        end
    end
end