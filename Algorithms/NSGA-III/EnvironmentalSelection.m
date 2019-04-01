function Population = EnvironmentalSelection(Population,N,Z,Zmin)
% The environmental selection of NSGA-III

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    if isempty(Zmin)
        Zmin = ones(1,size(Z,2));
    end

    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(Population.objs,Population.cons,N);
    Next = FrontNo < MaxFNo;
    
    %% Select the solutions in the last front
    Last   = find(FrontNo==MaxFNo);%最后一层解所在的位置
    Choose = LastSelection(Population(Next).objs,Population(Last).objs,N-sum(Next),Z,Zmin);%决定选择最后一层中哪些解
    Next(Last(Choose)) = true;
    % Population for next generation
    Population = Population(Next);
end

function Choose = LastSelection(PopObj1,PopObj2,K,Z,Zmin)
% Select part of the solutions in the last front

    PopObj = [PopObj1;PopObj2] - repmat(Zmin,size(PopObj1,1)+size(PopObj2,1),1);
    [N,M]  = size(PopObj);
    N1     = size(PopObj1,1);
    N2     = size(PopObj2,1);
    NZ     = size(Z,1);

    %% Normalization
    % Detect the extreme points 检测极值点
    Extreme = zeros(1,M);
    w       = zeros(M)+1e-6+eye(M);
    for i = 1 : M
        [~,Extreme(i)] = min(max(PopObj./repmat(w(i,:),N,1),[],2));
    end
    % Calculate the intercepts of the hyperplane constructed by the extreme
    % points and the axes 计算由极限点构造的超平面的截距
    Hyperplane = PopObj(Extreme,:)\ones(M,1);
    a = 1./Hyperplane;
    if any(isnan(a))
        a = max(PopObj,[],1)';
    end
    % Normalization 标准化
    PopObj = PopObj./repmat(a',N,1);
    
    %% Associate each solution with one reference point
    % Calculate the distance of each solution to each reference vector 
    % 将每个解与一个参考点关联，计算每个解到参考向量的距离
    Cosine   = 1 - pdist2(PopObj,Z,'cosine');
    Distance = repmat(sqrt(sum(PopObj.^2,2)),1,NZ).*sqrt(1-Cosine.^2);%若种群1和种群2规模总和为121，参考点的个数为91，那么Distance为一个121*91大小的二维数组，每一行代表个体到各个参考点的距离
    % Associate each solution with its nearest reference point
    % 将每个解与其最近的一个参考点关联
    [d,pi] = min(Distance',[],1);%计算与各个个体距离最近的参考点及其距离值

    %% Calculate the number of associated solutions except for the last front of each reference point
    % 计算除了最后一层以外的所有关联过解的参考点的个数
    rho = hist(pi(1:N1),1:NZ);
    
    %% Environmental selection
    Choose  = false(1,N2);
    Zchoose = true(1,NZ);
    % Select K solutions one by one
    while sum(Choose) < K
        % Select the least crowded reference point
        % 选择最不拥挤的参考点
        Temp = find(Zchoose);
        Jmin = find(rho(Temp)==min(rho(Temp))); %最不拥挤的参考点的集合
        j    = Temp(Jmin(randi(length(Jmin)))); %选定的不拥挤的参考点
        I    = find(Choose==0 & pi(N1+1:end)==j); %找到一个未被选过的，同时选择与参考点j对应的个体I
        % Then select one solution associated with this reference point
        if ~isempty(I)
            % 之前未被选中过
            if rho(j) == 0
                [~,s] = min(d(N1+I));%选择距离参考点最近的个体
            else
                s = randi(length(I));
            end
            Choose(I(s)) = true;%设置标志，不再选择最后一层的该参考点
            rho(j) = rho(j) + 1;
        else
            Zchoose(j) = false;%选择后将不再选择
        end
    end
end