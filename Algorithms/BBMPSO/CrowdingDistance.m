function [ crowd_value ] = CrowdingDistance( Population )
%计算种群中元素的拥挤距离值
%Population：种群

    num = size(Population,2);   %种群规模
    objs = Population.objs;
    M = size(objs,2); %目标个数
    crowd_value = zeros(1,num);
    lim_f = zeros(M,2);

    %粒子的数量大于等于2才能计算拥挤距离值，否则拥挤距离值取1
    if num >= 2
        for i = 1:M
            lim_f(i,1) = min(objs(:,i));         %目标函数i的极小值
            lim_f(i,2) = max(objs(:,i));        %目标函数i的极大值
        end
        DD = zeros(num,M);
        %遍历目标函数
        for i = 1:M
            [~,ind] = sort(objs(:,i));        %将目标函数i的值排序，ind：索引；val：值
            %遍历储备集中的粒子
            for j = 1:num
                if j == 1
                    DD(ind(j),i) = inf;     %设置目标函数i的极小值的拥挤距离值为无穷大
                elseif j == num
                    DD(ind(j),i) = inf;     %设置目标函数i的极大值的拥挤距离值为无穷大
                else
                    DD(ind(j),i) = objs(ind(j+1),i) - objs(ind(j-1),i)/2*(lim_f(i,2) - lim_f(i,1));     %计算目标函数i非边界点的拥挤距离值
                end
            end
        end
        for jj = 1:num
            crowd_value(jj) = sum(DD(jj,:));        %计算各个粒子的拥挤距离值
        end
    else
        crowd_value(1) = 1;
    end
end