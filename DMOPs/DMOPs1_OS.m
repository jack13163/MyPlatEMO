function varargout = DMOPs1_OS(Operation,Global,input)
% <problem> <DMOPs>
% ��Ŀ��ԭ�͵���(ʵ��1)

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
            Global.M        = 4;
            Global.M        = 4;
            Global.D        = 50;
            Global.lower    = zeros(1,Global.D);
            Global.upper    = ones(1,Global.D);
            %���ؽ���PF
            problem = func2str(Global.problem);
            pf_filepath = sprintf('Data/%s/%s.mat',problem,problem);
            if ~exist(pf_filepath,'file')
                fprintf('pf�����ļ������ڣ�\n');
            else
                pf = load(pf_filepath);
                [PF] = pf.PF.objs;
            end
            
            PopDec    = rand(input,Global.D);
            varargout = {PopDec};
        case 'value'
            pop     = input;
            [N,w] = size(pop);
            [eff] = fat_1(pop);
           	PopObj(:,1) = eff(1:N,w+1);
          	PopObj(:,2) = eff(1:N,w+2);
          	PopObj(:,3) = eff(1:N,w+3); 
           	PopObj(:,4) = eff(1:N,w+4); 
            
            PopCon = [];
            
            varargout = {eff(:,1:w),PopObj,PopCon};  %input�����߱�����PopObj��Ŀ�꺯��ֵ��PopCon��Լ������ֵ
        case 'PF'
            RefPoint  = PF;
            varargout = {RefPoint};
    end
end