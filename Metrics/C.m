function [ ratio ] = C( A,B )
%% Cָ�����
    B_num = size(B,1);
    count = 0;
    for i=1:B_num
        y = B(i,:);
        if Domain( y,A )
            count = count + 1;
        end
    end
    ratio = count/B_num;
end

%% �ж�X���Ƿ����һ����x֧��y
function [ flag ] = Domain( y,X )
    %��ȡ�ⲿ�����������ӵĸ���
    X_num = size(X,1);
    M = length(y);
    flag = false;
    %�ж��������ⲿ�����������ӵ�֧���ϵ
    for i = 1:X_num
        bb1 = 0;
        bb2 = 0;
        for j = 1:M     %x = X(i,:)
            aa1 = X(i,j);
            aa2 = y(j);
            if aa1 < aa2
                bb1 = bb1 + 1;
            elseif aa2 == aa1
                bb2 = bb2 + 1;
            end
        end
        %�ж�֧���ϵ
        if bb1 == M || (bb2 > 0 && bb1 == M - bb2)        %x֧��y
            flag = true;
            break;
        end
    end
end