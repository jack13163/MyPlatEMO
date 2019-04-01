function [ ratio ] = C( A,B )
%% C指标计算
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

%% 判断X中是否存在一个解x支配y
function [ flag ] = Domain( y,X )
    %获取外部储备集中粒子的个数
    X_num = size(X,1);
    M = length(y);
    flag = false;
    %判断粒子与外部储备集中粒子的支配关系
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
        %判断支配关系
        if bb1 == M || (bb2 > 0 && bb1 == M - bb2)        %x支配y
            flag = true;
            break;
        end
    end
end