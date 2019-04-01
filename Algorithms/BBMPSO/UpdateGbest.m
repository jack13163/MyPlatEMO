function [ Gbest ] = UpdateGbest( Archive, N )
%UPDATEGBEST 
    %% Find the non-dominated solutions
    Gbest_crowd_val = CrowdingDistance(Archive);
    num = size(Archive,2);
    %遍历所有的粒子，求粒子的全局引导者
    for i = 1:N                                                    %粒子个数为popsize
        a1 = randi(num,1);                           %从外部储备集中随机选择两个粒子
        a2 = randi(num,1);
        %比较两个粒子的拥挤距离值，选择拥挤距离值大的外部储备集中的粒子作为该粒子的全局引导者
        if Gbest_crowd_val(a1) >= Gbest_crowd_val(a2)
            Gbest(i) = Archive(a1);
        else
            Gbest(i) = Archive(a2);
        end
    end
end

