%% ɸѡ���
function [  ] = GetFiltedResult(  )
    a=load('F:\��Ŀ���Ż�����\ʵ�顪������\PlatEMO v1.5 (2017-12)\Data\DMOPs4_OS\DMOPs4_OS.mat');
        
    %����6�����Ž�
    %filted_optimization_result = UpdateArchive(a.PF,[],6,2);
    filted_optimization_result = a.PF([25]);
    
    %��ʾ��ϸ����
    title = {'������',	'���͹�','��ʼʱ��','����ʱ��','ԭ������','�ٶ�','ԭ�����'};
    for i=1:length(filted_optimization_result)
        objs = [filted_optimization_result(i).objs,0,0;40,35,12,9,168.32,0,0];
        filted_optimization_result(i).objs
        [~,schduleplan] = fat_4(filted_optimization_result(i).decs);
        title_str = sprintf('optimal scheduling %d',i);
        figure;
        printgante(schduleplan{1},title_str);
        filename = sprintf('F:/��Ŀ���Ż�����/ʵ�顪������/PlatEMO v1.5 (2017-12)/Data/DMOPs4_OS/���ŵ���_%d.xls',i);
        if exist(filename,'file')
            delete(filename);
        end
        save2excel([schduleplan{1};objs],title,filename);
    end
    
    %�����״�ͼ
    figure;    
    data=filted_optimization_result.objs;
    limit=[15,50;
        10,50;
        10,15;
        5,11;
        140,180];%�����᷶Χ
    clf;    
    label = {'Pipeline mixing cost','Bottom mixing cost','Number of oil tank switching','Number of charge tanks','Energy costs'};
    for i=1:size(data,1)
        Draw_radar(data(i,:),limit, label);
        if i~=size(data,1)
            figure;
        end
    end
end

function Archive = UpdateArchive( Population,OldArchive,N,div)
% Update the archive

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Find the non-dominated solutions
    Archive = CombineArchive(Population,OldArchive);
    
    %% Grid-based retention
    if length(Archive) > N
        Del = Delete(Archive.objs,length(Archive)-N,div);
        Archive(Del) = [];
    end
end

function Del = Delete(PopObj,K,div)   
    N = size(PopObj,1);

    %% Calculate the grid location of each solution
    fmax = max(PopObj,[],1);%�����Ŀ������ֵ
    fmin = min(PopObj,[],1);%�����Ŀ�����Сֵ
    d    = (fmax-fmin)/div;%ÿһ���̶ȵĳ���
    GLoc = floor((PopObj-repmat(fmin,N,1))./repmat(d,N,1));
    GLoc(GLoc>=div)   = div - 1;
    GLoc(isnan(GLoc)) = 0;

    %% Calculate the crowding degree of each grid
    [~,~,Site] = unique(GLoc,'rows');
    CrowdG     = hist(Site,1:max(Site));

    %% Delete K solutions
    Del = false(1,N);
    while sum(Del) < K
        % Select the most crowded grid
        maxGrid = find(CrowdG==max(CrowdG));
        Temp    = randi(length(maxGrid));
        Grid    = maxGrid(Temp);
        % And delete one solution randomly from the grid
        InGrid  = find(Site==Grid);
        Temp    = randi([1,length(InGrid)]);
        p       = InGrid(Temp);
        Del(p)  = true;
        Site(p) = NaN;
        CrowdG(Grid) = CrowdG(Grid) - 1;
    end
end