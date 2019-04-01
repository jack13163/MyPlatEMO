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
    div = Global.D;                               %ÿһ��Ŀ��Ŀ̶ȸ���
    Na = 10;                                         %������������

    %% Generate random population
    Population = Global.Initialization();
    Pbest      = Population;
    Archive    = UpdateArchive(Population,[],Na,div);%population�Ǹ���ļ��ϣ�����ÿ�������к����ĸ�����dec��obj��con��add
    
    %% Optimization
    while Global.NotTermination(Archive) %��ʾ���������е���Ϣ
        %ȷ��ȫ��������
        Gbest      = UpdateGbest(Archive, Global.N);
        %������һ����Ⱥ
        Population = Global.Variation([Population,Pbest,Gbest],Global.N);
        %���¸���������
        Pbest      = UpdatePbest(Pbest,Population);
        %���´�����
        Archive    = UpdateArchive(Population,Archive,Na,div);
    end
end

