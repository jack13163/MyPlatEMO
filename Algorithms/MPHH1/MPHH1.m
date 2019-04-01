function MPHH1(Global)
% <algorithm> MPHH1 <A>
% 多种群进化算法
% div --- 3 --- The number of divisions in each objective
% imgrate --- 0.1 --- The number of individual in each imgrate
% dc --- 10 --- 种群中熵值平稳的次数

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    %% Parameter setting
    [div,imgrate,dc] = Global.ParameterSet(3,0.1,10);      %获取用户配置参数
    MP = 2;                                              %子种群的个数（必须为偶数）
    POPNUM = fix(Global.N/MP);          %每个子种群中个体的个数
    IMNUM = fix(POPNUM*imgrate);     %子种群移民个体数量
    Global.N = MP * POPNUM;              %重新定义种群大小
    LastPopDiv = zeros(1,MP);              %前一代种群个体方差
    PopDiv = zeros(1,MP);                     %当前种群个体方差
    CountDiv = zeros(1,MP);                  %当前种群个体方差
    
    %% Generate random population
    Population = Global.Initialization(Global.N);
    FrontNo = zeros(1,Global.N);            %非支配等级（指标1）
    CrowdDis = zeros(1,Global.N);         %拥挤距离（指标2）
    Pbest  = Population;
    Offspring  = Population;    
    Target  = UpdateArchive(Population,[],Global.N,div);
    
    %% 各个种群初始化
    for i=1:MP
        %种群个体
        bird = ((i-1)*POPNUM+1):i*POPNUM;
        %% 不同的子种群使用不同的策略
        if mod(i,2)~=0
           	[~,FrontNo(bird),CrowdDis(bird)] = EnvironmentalSelection(Population(bird),POPNUM);     
        else
            Pbest(bird)      = Population(bird);
            Archive    = UpdateArchive(Population(bird),[],Global.N,div);
        end
    end
    
    %% Optimization
    while Global.NotTermination(Target) 
        for i=1:MP
            %子种群个体
            bird = ((i-1)*POPNUM+1):i*POPNUM;
            %% 不同的子种群使用不同的策略
            if mod(i,2)~=0
                MatingPool = TournamentSelection(2,POPNUM,FrontNo(bird),-CrowdDis(bird));
                Offspring(bird)  = Global.Variation(Population(MatingPool),POPNUM,@EAreal);
                [Population(bird),FrontNo(bird),CrowdDis(bird)] = EnvironmentalSelection([Population(bird),Offspring(bird)],POPNUM);
            else
                REP        = REPSelection(Archive.objs,POPNUM,div);                             %确定全局引导者
                Population(bird) = Global.Variation([Population(bird),Pbest(bird),Archive(REP)],POPNUM,@BBPSO);               %生成下一代种群
                Pbest(bird)      = UpdatePbest(Pbest(bird),Population(bird));                   %更新个体引导者
                %更新储备集
                Archive    = UpdateArchive(Population(bird),Archive,Global.N,div);
            end
            PopDiv(i) = GetPopDiv(Global,Population(bird));
            %% 累计群体熵值
            if PopDiv(i) <= LastPopDiv(i)
                CountDiv(i) = CountDiv(i) + 1;
            else
                CountDiv(i) = 0;
            end
        end
        
        %% 累计熵值达到指定次数
       	if any(CountDiv >= dc)
         	%% 迁移
          	Population = Immigrant(Population,MP,IMNUM,3);
            CountDiv = zeros(1,MP);
      	end
                
        %% 记录种群的熵值
        LastPopDiv = PopDiv;
        
        %% 更新全局储备集
        Target    = UpdateArchive(Population,Target,Global.N,div);                   
    end
end