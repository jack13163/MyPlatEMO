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
    Last   = find(FrontNo==MaxFNo);%���һ������ڵ�λ��
    Choose = LastSelection(Population(Next).objs,Population(Last).objs,N-sum(Next),Z,Zmin);%����ѡ�����һ������Щ��
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
    % Detect the extreme points ��⼫ֵ��
    Extreme = zeros(1,M);
    w       = zeros(M)+1e-6+eye(M);
    for i = 1 : M
        [~,Extreme(i)] = min(max(PopObj./repmat(w(i,:),N,1),[],2));
    end
    % Calculate the intercepts of the hyperplane constructed by the extreme
    % points and the axes �����ɼ��޵㹹��ĳ�ƽ��Ľؾ�
    Hyperplane = PopObj(Extreme,:)\ones(M,1);
    a = 1./Hyperplane;
    if any(isnan(a))
        a = max(PopObj,[],1)';
    end
    % Normalization ��׼��
    PopObj = PopObj./repmat(a',N,1);
    
    %% Associate each solution with one reference point
    % Calculate the distance of each solution to each reference vector 
    % ��ÿ������һ���ο������������ÿ���⵽�ο������ľ���
    Cosine   = 1 - pdist2(PopObj,Z,'cosine');
    Distance = repmat(sqrt(sum(PopObj.^2,2)),1,NZ).*sqrt(1-Cosine.^2);%����Ⱥ1����Ⱥ2��ģ�ܺ�Ϊ121���ο���ĸ���Ϊ91����ôDistanceΪһ��121*91��С�Ķ�ά���飬ÿһ�д�����嵽�����ο���ľ���
    % Associate each solution with its nearest reference point
    % ��ÿ�������������һ���ο������
    [d,pi] = min(Distance',[],1);%��������������������Ĳο��㼰�����ֵ

    %% Calculate the number of associated solutions except for the last front of each reference point
    % ����������һ����������й�������Ĳο���ĸ���
    rho = hist(pi(1:N1),1:NZ);
    
    %% Environmental selection
    Choose  = false(1,N2);
    Zchoose = true(1,NZ);
    % Select K solutions one by one
    while sum(Choose) < K
        % Select the least crowded reference point
        % ѡ���ӵ���Ĳο���
        Temp = find(Zchoose);
        Jmin = find(rho(Temp)==min(rho(Temp))); %�ӵ���Ĳο���ļ���
        j    = Temp(Jmin(randi(length(Jmin)))); %ѡ���Ĳ�ӵ���Ĳο���
        I    = find(Choose==0 & pi(N1+1:end)==j); %�ҵ�һ��δ��ѡ���ģ�ͬʱѡ����ο���j��Ӧ�ĸ���I
        % Then select one solution associated with this reference point
        if ~isempty(I)
            % ֮ǰδ��ѡ�й�
            if rho(j) == 0
                [~,s] = min(d(N1+I));%ѡ�����ο�������ĸ���
            else
                s = randi(length(I));
            end
            Choose(I(s)) = true;%���ñ�־������ѡ�����һ��ĸòο���
            rho(j) = rho(j) + 1;
        else
            Zchoose(j) = false;%ѡ��󽫲���ѡ��
        end
    end
end