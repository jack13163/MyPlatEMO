function varargout = DMOPs3_OS(Operation,Global,input)
% <problem> <DMOPs>
% 五目标原油调度(实例1)

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
persistent PF;

    switch Operation
        case 'init'
            Global.M        = 5;
            Global.M        = 5;
            Global.D        = 75;
            Global.lower    = zeros(1,Global.D);
            Global.upper    = ones(1,Global.D);
            %加载近似PF
            problem = func2str(Global.problem);
            pf_filepath = sprintf('Data/%s/%s.mat',problem,problem);
            if ~exist(pf_filepath,'file')
                fprintf('pf数据文件不存在！\n');
            else
                pf = load(pf_filepath);
                [PF] = pf.PF.objs;                
            end
            
            PopDec    = rand(input,Global.D);
            varargout = {PopDec};
        case 'value'
            pop     = input;
            [N,w] = size(pop);
            [eff] = fat_3(pop);
           	PopObj(:,1) = eff(1:N,w+1);         %管道混合成本
          	PopObj(:,2) = eff(1:N,w+2);         %罐底混合成本
          	PopObj(:,3) = eff(1:N,w+3);      	%蒸馏塔的供油罐切换次数
           	PopObj(:,4) = eff(1:N,w+4);         %供油罐个数
          	PopObj(:,5) = eff(1:N,w+5);         %能耗成本
            
            PopCon = [];
            
            varargout = {eff(:,1:w),PopObj,PopCon};  %input：决策变量；PopObj：目标函数值；PopCon：约束函数值
        case 'PF'
            RefPoint  = PF;
            varargout = {RefPoint};
    end
end
