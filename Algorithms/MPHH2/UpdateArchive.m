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
    fmax = max(PopObj,[],1);%求各个目标的最大值
    fmin = min(PopObj,[],1);%求各个目标的最小值
    d    = (fmax-fmin)/div;%每一个刻度的长度
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