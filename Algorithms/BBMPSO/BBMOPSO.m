function BBMOPSO( Global )
% <algorithm> <BBMOPSO>
% div --- 10 --- The number of divisions in each objective
% Bare bone Multi-objective particle swarm optimization

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    div = Global.D;                               %每一个目标的刻度个数
    Na = 10;                                         %储备集的容量

    %% Generate random population
    Population = Global.Initialization();
    Pbest      = Population;
    Archive    = UpdateArchive(Population,[],Na,div);%population是个体的集合，其中每个个体中含有四个对象：dec，obj，con，add
    
    %% Optimization
    while Global.NotTermination(Archive) %显示储备集合中的信息
        %确定全局引导者
        Gbest      = UpdateGbest(Archive, Global.N);
        %生成下一代种群
        Population = Global.Variation([Population,Pbest,Gbest],Global.N);
        %更新个体引导者
        Pbest      = UpdatePbest(Pbest,Population);
        %更新储备集
        Archive    = UpdateArchive(Population,Archive,Na,div);
    end
end

