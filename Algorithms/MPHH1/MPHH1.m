function MPHH1(Global)
% <algorithm> MPHH1 <A>
% ����Ⱥ�����㷨
% div --- 3 --- The number of divisions in each objective
% imgrate --- 0.1 --- The number of individual in each imgrate
% dc --- 10 --- ��Ⱥ����ֵƽ�ȵĴ���

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    %% Parameter setting
    [div,imgrate,dc] = Global.ParameterSet(3,0.1,10);      %��ȡ�û����ò���
    MP = 2;                                              %����Ⱥ�ĸ���������Ϊż����
    POPNUM = fix(Global.N/MP);          %ÿ������Ⱥ�и���ĸ���
    IMNUM = fix(POPNUM*imgrate);     %����Ⱥ�����������
    Global.N = MP * POPNUM;              %���¶�����Ⱥ��С
    LastPopDiv = zeros(1,MP);              %ǰһ����Ⱥ���巽��
    PopDiv = zeros(1,MP);                     %��ǰ��Ⱥ���巽��
    CountDiv = zeros(1,MP);                  %��ǰ��Ⱥ���巽��
    
    %% Generate random population
    Population = Global.Initialization(Global.N);
    FrontNo = zeros(1,Global.N);            %��֧��ȼ���ָ��1��
    CrowdDis = zeros(1,Global.N);         %ӵ�����루ָ��2��
    Pbest  = Population;
    Offspring  = Population;    
    Target  = UpdateArchive(Population,[],Global.N,div);
    
    %% ������Ⱥ��ʼ��
    for i=1:MP
        %��Ⱥ����
        bird = ((i-1)*POPNUM+1):i*POPNUM;
        %% ��ͬ������Ⱥʹ�ò�ͬ�Ĳ���
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
            %����Ⱥ����
            bird = ((i-1)*POPNUM+1):i*POPNUM;
            %% ��ͬ������Ⱥʹ�ò�ͬ�Ĳ���
            if mod(i,2)~=0
                MatingPool = TournamentSelection(2,POPNUM,FrontNo(bird),-CrowdDis(bird));
                Offspring(bird)  = Global.Variation(Population(MatingPool),POPNUM,@EAreal);
                [Population(bird),FrontNo(bird),CrowdDis(bird)] = EnvironmentalSelection([Population(bird),Offspring(bird)],POPNUM);
            else
                REP        = REPSelection(Archive.objs,POPNUM,div);                             %ȷ��ȫ��������
                Population(bird) = Global.Variation([Population(bird),Pbest(bird),Archive(REP)],POPNUM,@BBPSO);               %������һ����Ⱥ
                Pbest(bird)      = UpdatePbest(Pbest(bird),Population(bird));                   %���¸���������
                %���´�����
                Archive    = UpdateArchive(Population(bird),Archive,Global.N,div);
            end
            PopDiv(i) = GetPopDiv(Global,Population(bird));
            %% �ۼ�Ⱥ����ֵ
            if PopDiv(i) <= LastPopDiv(i)
                CountDiv(i) = CountDiv(i) + 1;
            else
                CountDiv(i) = 0;
            end
        end
        
        %% �ۼ���ֵ�ﵽָ������
       	if any(CountDiv >= dc)
         	%% Ǩ��
          	Population = Immigrant(Population,MP,IMNUM,3);
            CountDiv = zeros(1,MP);
      	end
                
        %% ��¼��Ⱥ����ֵ
        LastPopDiv = PopDiv;
        
        %% ����ȫ�ִ�����
        Target    = UpdateArchive(Population,Target,Global.N,div);                   
    end
end