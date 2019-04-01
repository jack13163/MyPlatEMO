function MPHH2(Global)
% <algorithm> MPHH2 <A>
% ����Ⱥ�Ŵ��㷨
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
    [imgrate,dc] = Global.ParameterSet(0.1,10);      %��ȡ�û����ò���    
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
    
    %% ������Ⱥ��ʼ��
    for i = 1:MP
        %% ��ͬ������Ⱥʹ�ò�ͬ�Ĳ���
        bird = ((i-1)*POPNUM+1):i*POPNUM;
        [~,FrontNo(bird),CrowdDis(bird)] = EnvironmentalSelection(Population(bird),POPNUM);
    end
    
    %% Optimization
    while Global.NotTermination(Population) %��ʾ���������е���Ϣ
        for i=1:MP
            %% ��ͬ������Ⱥʹ�ò�ͬ�Ĳ���
            bird = ((i-1)*POPNUM+1):i*POPNUM;
            MatingPool = TournamentSelection(2,POPNUM,FrontNo(bird),-CrowdDis(bird));
            if i==1
             	Offspring  = Global.Variation(Population(MatingPool),POPNUM,@EAreal,{1,20,1,20});
            elseif i==2
            	Offspring  = Global.Variation(Population(MatingPool),POPNUM,@DE,{1,0.5,1,20});
            end
           	[Population(bird),FrontNo(bird),CrowdDis(bird)] = EnvironmentalSelection([Population(bird),Offspring],POPNUM);
            %% ������Ⱥ�и������ֵ
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
          	Population = immigrant(Population,MP,IMNUM,4);
            CountDiv = zeros(1,MP);
      	end
                
        %% ��¼��Ⱥ����ֵ
        LastPopDiv = PopDiv;
    end
end