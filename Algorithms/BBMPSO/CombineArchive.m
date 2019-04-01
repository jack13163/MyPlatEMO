function [ Archive ] = CombineArchive( Population,OldArchive )
%   合并储备集
%   Detailed explanation goes here
    %% 若种群为空，则直接返回原储备集
    if isempty(Population)
        Archive = OldArchive;
        return;
    end
        
    %% 原来储备集为空
    if isempty(OldArchive)
        PopulationObj = Population.objs;%目标矩阵
        for i = 1:size(Population,2)
            %找到第一个可行的解，作为储备集的初始解，并更新种群为剩余解
            if PopulationObj(i,1) ~= inf
                Archive = Population(i);
                Population = Population(i+1:end);
                break;
            end
        end
    else
        Archive = OldArchive;
    end
    
	M = size(Population.objs,2);%目标个数
    
    %% 原储备集非空
    for i = 1:size(Population,2)
        %将每个候选者与外部储备集中的元素进行比较，更新储备集
        Archive = up_vac0(Population(i),Archive,M);
    end
end

%% 更新外部储备集
function [ new_AC ] = up_vac0( par_eff,old_AC,M )
%par_eff：候选粒子
%old_AC：外部储备集
%new_AC：新的外部储备集

    %获取外部储备集中粒子的个数
    ss1 = size(old_AC,2);
    old_ACobjs = old_AC.objs;
    par_effobj = par_eff.objs;
    %判断粒子与外部储备集中粒子的支配关系
    for i = 1:ss1
        bb1 = 0;
        bb2 = 0;
        for j = 1:M     %目标函数的个数M
            aa1 = old_ACobjs(i,j);
            aa2 = par_effobj(1,j);
            if aa2 < aa1
                bb1 = bb1 + 1;
            elseif aa2 == aa1
                bb2 = bb2 + 1;
            end
        end
        %判断支配关系
        if bb1 == M         %候选者支配old_AC(i,:)，即候选者支配外部储备集中的当前粒子
            old_ACobjs(i,:) = inf;          %删除外部储备集中被候选者支配的粒子
        elseif bb2 > 0 && bb1 == M - bb2       %候选者支配old_AC(i,:)，即候选者支配外部储备集中的当前粒子
            old_ACobjs(i,:) = inf;
        elseif bb1 + bb2 == 0       %候选者不支配old_AC(i,:)
            par_effobj(1,:) = inf;         %删除候选者粒子，候选者粒子被淘汰，不能进入外部储备集
            break;      %当前候选者支配于外部储备集中的当前粒子，退出
        elseif bb2 ~= 0 && bb1 == 0         %候选者不支配old_AC(i,:)
            par_effobj(1,:) = inf;
            break;
        end
    end
    %将候选者加入到临时外部储备集中，临时外部储备集中被支配的粒子所在的行变为inf，即无穷大
    tmp_AC = [par_eff,old_AC]; 
    tmp_ACobjs = [par_effobj;old_ACobjs];
    new_AC = [];
    %遍历合并后的临时外部储备集
    for i = 1:ss1+1
        %判断当前粒子是否被其他粒子所支配，若否，则将其纳入新的外部储备集张
        if tmp_ACobjs(i,1) ~= inf
            new_AC = [new_AC,tmp_AC(i)];
        end
    end
end
